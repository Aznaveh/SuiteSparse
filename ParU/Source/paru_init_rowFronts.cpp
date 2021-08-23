/** =========================================================================  /
 * =======================  paru_init_rowFronts  ============================  /
 * ========================================================================== */
/*!  @brief  Initializing row fronts; fronts will be assembled later.
 *         Initializing Row and column tuple lists:
 *            allocating memory and updating lists and initializing matrix
 *            structre
 *        Assemble each front:
 *            Adding numerical values, allocating data
 *            updating the list
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

ParU_ResultCode paru_init_rowFronts(
    paru_matrix **paruMatInfo_handle,  // in/out
                                       // inputs, not modified
    cholmod_sparse *A,
    int scale,  // scales the matrix if > 0
    // symbolic analysis
    paru_symbolic *LUsym)
{
    mallopt(M_MMAP_MAX, 0);                // disable mmap; it's too slow
    mallopt(M_TRIM_THRESHOLD, -1);         // disable sbrk trimming
    mallopt(M_TOP_PAD, 16 * 1024 * 1024);  // increase padding to speedup malloc

    DEBUGLEVEL(0);
    if (!A->packed)
    {
        printf("A is not packed; Wrong format \n");
        return PARU_INVALID;
    }

    if (LUsym == NULL)
    {
        printf("LUsym is NULL\n");
        return PARU_INVALID;
    }

    // initializing paruMat
    paru_matrix *paruMatInfo;
    paruMatInfo = (paru_matrix *)paru_alloc(1, sizeof(paru_matrix));
    if (paruMatInfo == NULL)
    {  // out of memory
        printf("Out of memory: paruMatInfo\n");
        // Nothing to be freed
        return PARU_OUT_OF_MEMORY;
    }
    *paruMatInfo_handle = paruMatInfo;

    Int m, n, nf;
    paruMatInfo->LUsym = LUsym;
    m = paruMatInfo->m = LUsym->m - LUsym->n1;
    n = paruMatInfo->n = LUsym->n - LUsym->n1;
    nf = LUsym->nf;
    paruMatInfo->panel_width = 32;

    work_struct *Work = paruMatInfo->Work =
        (work_struct *)paru_alloc(1, sizeof(work_struct));

    if (Work == NULL)
    {  // out of memory
        printf("Out of memory: Work\n");
        paru_freemat(&paruMatInfo);
        return PARU_OUT_OF_MEMORY;
    }
    // Memory allocations for paruMatInfo
    Int *rowMark = Work->rowMark = (Int *)paru_alloc(m + nf + 1, sizeof(Int));
    Int *elRow = Work->elRow = (Int *)paru_alloc(m + nf, sizeof(Int));
    Int *elCol = Work->elCol = (Int *)paru_alloc(m + nf, sizeof(Int));
    Int *rowSize = Work->rowSize = (Int *)paru_alloc(m, sizeof(Int));
    Int *row_degree_bound = paruMatInfo->row_degree_bound =
        (Int *)paru_alloc(m, sizeof(Int));
    tupleList *RowList = paruMatInfo->RowList =
        (tupleList *)paru_alloc(1, m * sizeof(tupleList));
    paruMatInfo->lacList = (Int *)paru_alloc(m + nf, sizeof(Int));
    paruMatInfo->frowCount = (Int *)paru_alloc(1, nf * sizeof(Int));
    paruMatInfo->fcolCount = (Int *)paru_alloc(1, nf * sizeof(Int));
    paruMatInfo->frowList = (Int **)paru_calloc(1, nf * sizeof(Int *));
    paruMatInfo->fcolList = (Int **)paru_calloc(1, nf * sizeof(Int *));
    paruMatInfo->partial_Us =  // Initialize with NULL
        (paru_fac *)paru_calloc(1, nf * sizeof(paru_fac));
    paruMatInfo->partial_LUs =  // Initialize with NULL
        (paru_fac *)paru_calloc(1, nf * sizeof(paru_fac));
    paruMatInfo->time_stamp = (Int *)paru_alloc(1, nf * sizeof(Int));

    std::vector<Int> **heapList = paruMatInfo->heapList =
        (std::vector<Int> **)paru_calloc(
            1, (m + nf + 1) * sizeof(std::vector<Int> *));
    paru_Element **elementList;
    elementList = paruMatInfo->elementList =  // Initialize with NULL
        (paru_Element **)paru_calloc(1, (m + nf + 1) * sizeof(paru_Element));

    if (rowMark == NULL || elRow == NULL || elCol == NULL || rowSize == NULL ||
        paruMatInfo->lacList == NULL || RowList == NULL ||
        row_degree_bound == NULL || elementList == NULL ||
        paruMatInfo->frowCount == NULL || paruMatInfo->fcolCount == NULL ||
        paruMatInfo->frowList == NULL || paruMatInfo->fcolList == NULL ||
        paruMatInfo->partial_Us == NULL || paruMatInfo->partial_LUs == NULL ||
        paruMatInfo->time_stamp == NULL || heapList == NULL)
    {
        paru_freemat(&paruMatInfo);
        return PARU_OUT_OF_MEMORY;
    }

    // Initializations
    PRLEVEL(1, ("%% $RowList =%p\n", RowList));
    paru_memset(rowSize, -1, m * sizeof(Int));
    PRLEVEL(1, ("%% rowSize pointer=%p size=%ld \n", rowSize, m * sizeof(Int)));

    PRLEVEL(1, ("%% rowMark pointer=%p size=%ld \n", rowMark,
                (m + nf) * sizeof(Int)));

    paru_memset(elRow, -1, (m + nf) * sizeof(Int));
    PRLEVEL(1, ("%% elRow=%p\n", elRow));

    paru_memset(elCol, -1, (m + nf) * sizeof(Int));
    PRLEVEL(1, ("%% elCol=%p\n", elCol));

    PRLEVEL(1, ("%% Work =%p\n ", Work));

#ifdef COUNT_FLOPS
    // flop count info init
    paruMatInfo->flp_cnt_dgemm = paruMatInfo->flp_cnt_trsm =
        paruMatInfo->flp_cnt_dger = 0.0;
    paruMatInfo->flp_cnt_real_dgemm = 0.0;
#endif

    PRLEVEL(1, ("%% m=%ld, n=%ld\n", m, n));
    // RowList, ColList and elementList are place holders
    // pointers to pointers that are allocated
    if (m == 0 || n == 0)
    {
        printf("%%The dimension of matrix is zero: %ld x %ld \n", m, n);
        paru_free(1, sizeof(paru_matrix), paruMatInfo);
        return PARU_INVALID;
    }

    Int snz = LUsym->snz;
    double *Sx = LUsym->Sx;
    Int *Sp = LUsym->Sp;
    Int *Sj = LUsym->Sj;

    //~~~~~~~~~~ scaling the A matrix
    scale = 1;
    if (scale)
    {
        double *max_row = (double *)paru_calloc(m, sizeof(double));
        if (max_row == NULL)
        {  // out of memory
            paru_freemat(&paruMatInfo);
            printf("of memory: max_row\n");
            return PARU_OUT_OF_MEMORY;
        }

        for (Int row = 0; row < m; row++)
        {
            double max = fabs(Sx[Sp[row]]);
            for (Int p = Sp[row] + 1; p < Sp[row + 1]; p++)
            {
                // el_colrowNum[j++] = Sx[p]/scale[p];
                max = MAX(fabs(Sx[p]), max);
                PRLEVEL(1, ("Sj[%ld] =%ld Sx[%ld]=%lf\n", p, Sj[p], p, Sx[p]));
                // for Matlab
                PRLEVEL(0, ("%ld,%ld, %.16lf;\n", row + 1, Sj[p] + 1, Sx[p]));
            }
            PRLEVEL(1, ("max in row %ld is %lf\n", row, max));
            max_row[row] = max;
        }

        // TODO Do the division when puting stuff into element rows
        paruMatInfo->scale_row = max_row;

#ifndef NDEBUG
        Int p = -1;
        PRLEVEL(p, ("%% scale =[ "));
        for (Int row = 0; row < m; row++) PRLEVEL(p, ("%lf ", max_row[row]));
        PRLEVEL(p, ("]\n"));
#endif
    }
    else
        paruMatInfo->scale_row = NULL;

    /// ------------------------------------------------------------------------
    // create S = A (p,q)', or S=A(p,q) if S is considered to be in row-form
    // -------------------------------------------------------------------------
#ifndef NDEBUG
    Int p = 1;
    PRLEVEL(p, ("\n%% Insid init row fronts\n"));
    PRLEVEL(p, ("%% Sp =\n%%"));
    for (Int i = 0; i <= m; i++) PRLEVEL(p, ("%ld ", Sp[i]));
    PRLEVEL(p, ("\n"));

    PRLEVEL(p, ("Sj =\n"));
    for (Int k = 0; k < snz; k++) PRLEVEL(p, ("%ld ", Sj[k]));
    PRLEVEL(p, ("\n"));

#endif

    // constants for initialzing lists
    Int slackRow = 2;

    PRLEVEL(0, ("InMatrix=[\n"));  // MATLAB matrix,

    // Activating comments after this parts will break the matlab input matrix
    // allocating row tuples, elements and updating column tuples

    ParU_ResultCode info;
    info = PARU_SUCCESS;
    //#pragma omp parallel shared (info)
    {
        //#pragma omp for
        for (Int row = 0; row < m; row++)
        {
            Int e = LUsym->row2atree[row];
            Int nrows = 1,
                ncols = Sp[row + 1] -
                        Sp[row];  // nrows and ncols of current front/row

            PRLEVEL(1, ("%% element %ld = %ld x %ld\n", e, nrows, ncols));

            row_degree_bound[row] = ncols;  // Initialzing row degree

            paru_Element *curEl = elementList[e] =
                paru_create_element(nrows, ncols, 0);
            if (curEl == NULL)
            {  // out of memory
                paru_freemat(&paruMatInfo);
                printf("Out of memory: curEl\n");
                info = PARU_OUT_OF_MEMORY;
                //#pragma omp cancel for
                return PARU_OUT_OF_MEMORY;
            }

            rowMark[e] = 0;

            // My new is calling paru_alloc now; so there is no need
            std::vector<Int> *curHeap;
            curHeap = paruMatInfo->heapList[e] = new std::vector<Int>;
            PRLEVEL(1, ("%%Heap allocated %p id=%ld \n", curHeap, e));
            curHeap->push_back(e);

#ifndef NDEBUG  // Printing the pointers info
            Int p = 1;
            PRLEVEL(p, ("%% curEl = %p ", curEl));
            Int size = sizeof(paru_Element) +
                       sizeof(Int) * (2 * (nrows + ncols)) +
                       sizeof(double) * nrows * ncols;
            PRLEVEL(p, ("size= %ld", size));
            PRLEVEL(p, ("\n"));
#endif

            // Allocating Rowlist and updating its tuples
            RowList[row].list =
                (paru_Tuple *)paru_alloc(slackRow * nrows, sizeof(paru_Tuple));
            if (RowList[row].list == NULL)
            {  // out of memory
                paru_freemat(&paruMatInfo);
                printf("Out of memory: RowList[row].list \n");
                info = PARU_OUT_OF_MEMORY;
                //#pragma omp cancel for
                return PARU_OUT_OF_MEMORY;
            }
            RowList[row].numTuple = 0;
            RowList[row].len = slackRow;

            paru_Tuple rowTuple;
            rowTuple.e = e;
            rowTuple.f = 0;
            if (paru_add_rowTuple(RowList, row, rowTuple))
            {
                paru_freemat(&paruMatInfo);
                printf("Out of memory: add_rowTuple \n");
                info = PARU_OUT_OF_MEMORY;
                //#pragma omp cancel for
                return PARU_OUT_OF_MEMORY;
            }

            // Allocating elements
            Int *el_colrowIndex = colIndex_pointer(curEl);
            double *el_colrowNum = numeric_pointer(curEl);

            PRLEVEL(1, ("el_colrowIndex =%p, el_colrowNum = %p \n",
                        el_colrowIndex, el_colrowNum));

            Int j = 0;  // Index inside an element
            // TODO choosing p as a variable can shadow p in debug mode
            double *s = paruMatInfo->scale_row;
            double r_scale;
            if (s != NULL) r_scale = s[row];
            for (Int p = Sp[row]; p < Sp[row + 1]; p++)
            {
                // PRLEVEL(0, ("scale = %lf\t", r_scale));
                el_colrowIndex[j] = Sj[p];
                // TODO: adding the scale here
                // el_colrowNum[j++] = Sx[p]/scale[p];
                el_colrowNum[j++] = (scale == 0) ? Sx[p] : (Sx[p] / r_scale);
                PRLEVEL(1, ("Sj[%ld] =%ld Sx[%ld]=%lf scaled=%lf\n", p, Sj[p],
                            p, Sx[p], (Sx[p] / s[row])));
                // for Matlab
                PRLEVEL(0, ("%ld,%ld, %.16lf;\n", row + 1, Sj[p] + 1,
                            (scale == 0) ? Sx[p] : (Sx[p] / s[row])));
                // Sx[p]));
            }
            el_colrowIndex[j++] = row;  // initializing element row index
            paruMatInfo->lacList[e] = lac_el(elementList, e);
        }
    }

    PRLEVEL(0, ("];\n"));
    PRLEVEL(0, ("I = InMatrix(:,1);\n"));
    PRLEVEL(0, ("J = InMatrix(:,2);\n"));
    PRLEVEL(0, ("X = InMatrix(:,3);\n"));
    PRLEVEL(0, ("S = sparse(I,J,X);\n"));

    // Free here or if not wil be freed in paru_mem anyway
    paru_free(snz, sizeof(double), Sx);
    paru_free(snz, sizeof(Int), Sj);
    LUsym->Sx = NULL;
    LUsym->Sj = NULL;
    return info;
}
