////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_init_rowFronts  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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

ParU_Res paru_init_rowFronts(
    paru_matrix **paruMatInfo_handle,  // in/out
                                       // inputs, not modified
    cholmod_sparse *A,
    // symbolic analysis
    ParU_symbolic *Sym)
{
    //mallopt(M_TRIM_THRESHOLD, -1);         // disable sbrk trimming
    //mallopt(M_TOP_PAD, 16 * 1024 * 1024);  // increase padding to speedup malloc

    DEBUGLEVEL(-1);
    if (!A->packed)
    {
        printf("Paru: A is not packed; Wrong format \n");
        return PARU_INVALID;
    }

    if (Sym == NULL)
    {
        printf("Paru: Sym is NULL\n");
        return PARU_INVALID;
    }

    // initializing paruMat
    paru_matrix *paruMatInfo;
    paruMatInfo = (paru_matrix *)paru_alloc(1, sizeof(paru_matrix));
    if (paruMatInfo == NULL)
    {  // out of memory
        printf("Paru: out of memory, paruMatInfo\n");
        // Nothing to be freed
        return PARU_OUT_OF_MEMORY;
    }
    *paruMatInfo_handle = paruMatInfo;

    Int m, n, nf;
    paruMatInfo->Sym = Sym;
    m = paruMatInfo->m = Sym->m - Sym->n1;
    n = paruMatInfo->n = Sym->n - Sym->n1;
    nf = Sym->nf;
    paruMatInfo->panel_width = 32;
    paruMatInfo->res = PARU_SUCCESS;

    work_struct *Work = paruMatInfo->Work =
        (work_struct *)paru_alloc(1, sizeof(work_struct));

    if (Work == NULL)
    {  // out of memory
        printf("Paru: out of memory: Work\n");
        paru_freemat(&paruMatInfo);
        return PARU_OUT_OF_MEMORY;
    }

    Work->rowMark = Work->elRow = NULL;
    Work->elCol = Work->rowSize = NULL;
    paruMatInfo->row_degree_bound = NULL;
    paruMatInfo->RowList = NULL;
    paruMatInfo->lacList = NULL;
    paruMatInfo->frowCount = NULL;
    paruMatInfo->fcolCount = NULL;
    paruMatInfo->frowList = NULL;
    paruMatInfo->fcolList = NULL;
    paruMatInfo->partial_Us = NULL;
    paruMatInfo->partial_LUs = NULL;
    paruMatInfo->heapList = NULL;
    paruMatInfo->elementList = NULL;
    paruMatInfo->time_stamp = NULL;
    paruMatInfo->Diag_map = NULL;
    paruMatInfo->inv_Diag_map = NULL;

    if (nf == 0)
    {  // nothing to be done
        return PARU_SUCCESS;
    }
    // Memory allocations for paruMatInfo
    Int *rowMark = Work->rowMark = (Int *)paru_alloc(m + nf + 1, sizeof(Int));
    Int *elRow = Work->elRow = (Int *)paru_alloc(m + nf, sizeof(Int));
    Int *elCol = Work->elCol = (Int *)paru_alloc(m + nf, sizeof(Int));
    Int *rowSize = Work->rowSize = (Int *)paru_alloc(m, sizeof(Int));
    Int *row_degree_bound = paruMatInfo->row_degree_bound =
        (Int *)paru_alloc(m, sizeof(Int));
    ParU_TupleList *RowList = paruMatInfo->RowList =
        (ParU_TupleList *)paru_alloc(1, m * sizeof(ParU_TupleList));
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
    ParU_Element **elementList;
    elementList = paruMatInfo->elementList =  // Initialize with NULL
        (ParU_Element **)paru_calloc(1, (m + nf + 1) * sizeof(ParU_Element));
    Int *Diag_map = paruMatInfo->Diag_map = NULL;
    Int *inv_Diag_map = paruMatInfo->inv_Diag_map = NULL;
    if (Sym->strategy == PARU_STRATEGY_SYMMETRIC)
    {
        Diag_map = paruMatInfo->Diag_map =
            (Int *)paru_alloc(Sym->n, sizeof(Int));
        inv_Diag_map = paruMatInfo->inv_Diag_map =
            (Int *)paru_alloc(Sym->n, sizeof(Int));
#ifndef NDEBUG
        paru_memset(Diag_map, 0, Sym->n * sizeof(Int));
        paru_memset(inv_Diag_map, 0, Sym->n * sizeof(Int));
#endif
    }

    if (rowMark == NULL || elRow == NULL || elCol == NULL || rowSize == NULL ||
        paruMatInfo->lacList == NULL || RowList == NULL ||
        row_degree_bound == NULL || elementList == NULL ||
        paruMatInfo->frowCount == NULL || paruMatInfo->fcolCount == NULL ||
        paruMatInfo->frowList == NULL || paruMatInfo->fcolList == NULL ||
        paruMatInfo->partial_Us == NULL || paruMatInfo->partial_LUs == NULL ||
        paruMatInfo->time_stamp == NULL || heapList == NULL ||
        (Sym->strategy == PARU_STRATEGY_SYMMETRIC &&
         (Diag_map == NULL || inv_Diag_map == NULL)))
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
    paruMatInfo->flp_cnt_dgemm = 0.0;
    paruMatInfo->flp_cnt_trsm = 0.0;
    paruMatInfo->flp_cnt_dger = 0.0;
    paruMatInfo->flp_cnt_real_dgemm = 0.0;
#endif

    PRLEVEL(1, ("%% m=%ld, n=%ld\n", m, n));
    // RowList, ColList and elementList are place holders
    // pointers to pointers that are allocated
    if (m == 0 || n == 0)
    {
        printf("Paru: the dimension of matrix is zero: %ld x %ld \n", m, n);
        paru_free(1, sizeof(paru_matrix), paruMatInfo);
        return PARU_INVALID;
    }

    Int snz = Sym->snz;
    double *Sx = Sym->Sx;
    Int *Sp = Sym->Sp;
    Int *Sj = Sym->Sj;

    /// ------------------------------------------------------------------------
    // create S = A (p,q)', or S=A(p,q) if S is considered to be in row-form
    // -------------------------------------------------------------------------
#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("\n%% Insid init row fronts\n"));
    PRLEVEL(PR, ("%% Sp =\n%%"));
    for (Int i = 0; i <= m; i++) PRLEVEL(PR, ("%ld ", Sp[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Sj =\n"));
    for (Int k = 0; k < snz; k++) PRLEVEL(PR, ("%ld ", Sj[k]));
    PRLEVEL(PR, ("\n"));

#endif

    PRLEVEL(0, ("InMatrix=[\n"));  // MATLAB matrix,

    // copying Diag_map
    if (Diag_map)
    {
        #pragma omp taskloop default(none)\
        shared(Sym, Diag_map, inv_Diag_map) grainsize(512)
        for (Int i = 0; i < Sym->n; i++)
        {
            // paru_memcpy(Diag_map, Sym->Diag_map, (Sym->n) * sizeof(Int));
            Diag_map[i] = Sym->Diag_map[i];
            inv_Diag_map[Diag_map[i]] = i;
        }
#ifndef NDEBUG
        PR = -1;
        PRLEVEL(PR, ("init_row Diag_map (%ld) =\n", Sym->n));
        for (Int i = 0; i < MIN(64, Sym->n); i++)
            PRLEVEL(PR, ("%ld ", Diag_map[i]));
        PRLEVEL(PR, ("\n"));
        PRLEVEL(PR, ("inv_Diag_map =\n"));
        for (Int i = 0; i < MIN(64, Sym->n); i++)
            PRLEVEL(PR, ("%ld ", inv_Diag_map[i]));
        PRLEVEL(PR, ("\n"));
        for (Int i = 0; i < Sym->n; i++)
        {
            if (Diag_map[i] == -1)
                PRLEVEL(PR,
                        ("Diag_map[%ld] is not correctly initialized\n", i));

            if (inv_Diag_map[i] == -1)
                PRLEVEL(PR, ("inv_Diag_map[%ld] is not correctly initialized\n",
                             i));

            ASSERT(Diag_map[i] != -1);
            // ASSERT(inv_Diag_map[i] != -1);
        }
        PR = 1;
#endif
    }

    // Activating comments after this parts will break the matlab input matrix
    // allocating row tuples, elements and updating column tuples

    ParU_Res info;
    Int out_of_memory = 0;
    #pragma omp taskloop default(none) \
    shared(out_of_memory, Sym, Sp, row_degree_bound, elementList, m, \
            paruMatInfo, rowMark, RowList, Sj, Sx) grainsize(512)
    for (Int row = 0; row < m; row++)
    {
        Int e = Sym->row2atree[row];
        Int nrows = 1,
            ncols =
                Sp[row + 1] - Sp[row];  // nrows and ncols of current front/row
        //printf("%% element %ld = %ld x %ld\n", e, nrows, ncols);

        row_degree_bound[row] = ncols;  // Initialzing row degree

        ParU_Element *curEl = elementList[e] =
            paru_create_element(nrows, ncols, 0);
        if (curEl == NULL)
        {  // out of memory
            paru_freemat(&paruMatInfo);
            printf("Paru: Out of memory: curEl\n");
            #pragma omp atomic update
            out_of_memory += 1;
        }

        rowMark[e] = 0;

        // My new is calling paru_alloc
        std::vector<Int> *curHeap;
        try
        {
            curHeap = paruMatInfo->heapList[e] = new std::vector<Int>;
        }
        catch (std::bad_alloc const &)
        {  // out of memory
            paru_freemat(&paruMatInfo);
            printf("Paru: Out of memory: curHeap\n");
            #pragma omp atomic update
            out_of_memory += 1;
        }
        //printf("%%Heap allocated %p id=%ld \n", curHeap, e);

        curHeap->push_back(e);

#ifndef NDEBUG  // Printing the pointers info
        //printf ("%% curEl = %p ", curEl);
        //Int size = sizeof(ParU_Element) + sizeof(Int) * (2 * (nrows + ncols)) 
        //           + sizeof(double) * nrows * ncols;
        //printf("size= %ld\n", size);
#endif

        // constants for initialzing lists
        Int slackRow = 2;

        // Allocating Rowlist and updating its tuples
        RowList[row].list =
            (ParU_Tuple *)paru_alloc(slackRow * nrows, sizeof(ParU_Tuple));
        if (RowList[row].list == NULL)
        {  // out of memory
            paru_freemat(&paruMatInfo);
            printf("Paru: out of memory, RowList[row].list \n");
            #pragma omp atomic update
            out_of_memory += 1;
        }
        RowList[row].numTuple = 0;
        RowList[row].len = slackRow;

        ParU_Tuple rowTuple;
        rowTuple.e = e;
        rowTuple.f = 0;
        if (paru_add_rowTuple(RowList, row, rowTuple))
        {
            paru_freemat(&paruMatInfo);
            printf("Paru: out of memory, add_rowTuple \n");
            #pragma omp atomic update
            out_of_memory += 1;
        }

        // Allocating elements
        Int *el_colrowIndex = colIndex_pointer(curEl);
        double *el_colrowNum = numeric_pointer(curEl);

        //printf("el_colrowIndex =%p, el_colrowNum = %p \n",
        // el_colrowIndex, el_colrowNum);

        Int j = 0;  // Index inside an element
        for (Int p = Sp[row]; p < Sp[row + 1]; p++)
        {
            el_colrowIndex[j] = Sj[p];
            el_colrowNum[j++] = Sx[p];
            //printf("Sj[%ld] =%ld Sx[%ld]=%lf \n", p, Sj[p], p, Sx[p]);
            // for Matlab
            //printf("%ld,%ld, %.16lf;\n", row + 1, Sj[p] + 1, Sx[p]);
        }
        el_colrowIndex[j++] = row;  // initializing element row index
        paruMatInfo->lacList[e] = lac_el(elementList, e);
    }
    if (out_of_memory)
        info = PARU_OUT_OF_MEMORY;
    else 
        info = PARU_SUCCESS;

    PRLEVEL(0, ("];\n"));
    PRLEVEL(0, ("I = InMatrix(:,1);\n"));
    PRLEVEL(0, ("J = InMatrix(:,2);\n"));
    PRLEVEL(0, ("X = InMatrix(:,3);\n"));
    PRLEVEL(0, ("S = sparse(I,J,X);\n"));

    // Free here or if not wil be freed in paru_mem anyway
    paru_free(snz, sizeof(double), Sx);
    paru_free(snz, sizeof(Int), Sj);
    Sym->Sx = NULL;
    Sym->Sj = NULL;
    return info;
}
