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

ParU_Ret paru_init_rowFronts(
        paru_work *Work,
        ParU_Numeric **Num_handle,  // in/out
                                   // inputs, not modified
                             cholmod_sparse *A,
                             // symbolic analysis
                             ParU_Symbolic *Sym, ParU_Control *Control)
{
    // mallopt(M_TRIM_THRESHOLD, -1);         // disable sbrk trimming
    // mallopt(M_TOP_PAD, 16 * 1024 * 1024);  // increase padding to speedup
    // malloc

    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;

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

    // initializing Numeric
    ParU_Numeric *Num;
    Num = (ParU_Numeric *)paru_alloc(1, sizeof(ParU_Numeric));
    if (Num == NULL)
    {  // out of memory
        printf("Paru: out of memory, Num\n");
        // Nothing to be freed
        return PARU_OUT_OF_MEMORY;
    }
    *Num_handle = Num;

    Int m, nf;
    Work->Sym = Sym;
    m = Num->m = Sym->m - Sym->n1;
    nf = Sym->nf;
    Num->res = PARU_SUCCESS;
    Num->Control = Control;

    Int *row_degree_bound = Num->row_degree_bound = NULL;
    Num->frowCount = NULL;
    Num->fcolCount = NULL;
    Num->frowList = NULL;
    Num->fcolList = NULL;
    Num->partial_Us = NULL;
    Num->partial_LUs = NULL;
    Num->Sx = NULL;
    Num->Sux = NULL;
    Num->Slx = NULL;
    Num->Rs = NULL;

    //Workd DS
    Int *rowMark = Work->rowMark = NULL;
    Int *elRow = Work->elRow = NULL;
    Int *elCol = Work->elCol = NULL;
    Int *rowSize = Work->rowSize = NULL;
    Work->time_stamp = NULL;
    paru_tupleList *RowList = Work->RowList = NULL;
    Int *Diag_map = Work->Diag_map = NULL;
    Int *inv_Diag_map = Work->inv_Diag_map = NULL;
    paru_element **elementList = Work->elementList = NULL;
    Work->lacList = NULL;
    std::vector<Int> **heapList = Work->heapList = NULL;

    if (nf != 0)
    {
        // Memory allocations for Num
        rowMark = Work->rowMark = (Int *)paru_alloc(m + nf + 1, sizeof(Int));
        elRow = Work->elRow = (Int *)paru_alloc(m + nf, sizeof(Int));
        elCol = Work->elCol = (Int *)paru_alloc(m + nf, sizeof(Int));
        rowSize = Work->rowSize = (Int *)paru_alloc(m, sizeof(Int));
        row_degree_bound = Num->row_degree_bound =
            (Int *)paru_alloc(m, sizeof(Int));
        RowList = Work->RowList =
            (paru_tupleList *)paru_alloc(1, m * sizeof(paru_tupleList));
        Work->lacList = (Int *)paru_alloc(m + nf, sizeof(Int));
        Num->frowCount = (Int *)paru_alloc(1, nf * sizeof(Int));
        Num->fcolCount = (Int *)paru_alloc(1, nf * sizeof(Int));
        Num->frowList = (Int **)paru_calloc(1, nf * sizeof(Int *));
        Num->fcolList = (Int **)paru_calloc(1, nf * sizeof(Int *));
        Num->partial_Us =  // Initialize with NULL
            (ParU_Factors *)paru_calloc(1, nf * sizeof(ParU_Factors));
        Num->partial_LUs =  // Initialize with NULL
            (ParU_Factors *)paru_calloc(1, nf * sizeof(ParU_Factors));

        Work->time_stamp = (Int *)paru_alloc(1, nf * sizeof(Int));

        heapList = Work->heapList = (std::vector<Int> **)paru_calloc(
            1, (m + nf + 1) * sizeof(std::vector<Int> *));
        elementList = Work->elementList =  // Initialize with NULL
            (paru_element **)paru_calloc(1,
                                         (m + nf + 1) * sizeof(paru_element));
        if (Sym->strategy == PARU_STRATEGY_SYMMETRIC)
        {
            Diag_map = Work->Diag_map = (Int *)paru_alloc(Sym->n, sizeof(Int));
            inv_Diag_map = Work->inv_Diag_map =
                (Int *)paru_alloc(Sym->n, sizeof(Int));
#ifndef NDEBUG
            paru_memset(Diag_map, 0, Sym->n * sizeof(Int), Control);
            paru_memset(inv_Diag_map, 0, Sym->n * sizeof(Int), Control);
            PR = 2;
#endif
        }
    }

    Int snz = Sym->snz;
    double *Sx = NULL;
    Sx = Num->Sx = (double *)paru_alloc(snz, sizeof(double));
    Int *cSp = NULL;  // copy of Sp, temporary for making Sx
    cSp = (Int *)paru_alloc(m + 1, sizeof(Int));
    double *Sux = NULL;
    Int *cSup = NULL;  // copy of Sup temporary for making Sux
    Int cs1 = Sym->cs1;
    Int rs1 = Sym->rs1;
    Int sunz = 0;
    if (cs1 > 0)
    {
        sunz = Sym->ustons.nnz;
        Sux = (double *)paru_alloc(sunz, sizeof(double));
        cSup = (Int *)paru_alloc(cs1 + 1, sizeof(Int));
    }
    Num->Sux = Sux;
    double *Slx = NULL;
    Int *cSlp = NULL;  // copyf of Slp temporary, for making Slx
    Int slnz = 0;
    if (rs1 > 0)
    {
        slnz = Sym->lstons.nnz;
        Slx = (double *)paru_alloc(slnz, sizeof(double));
        cSlp = (Int *)paru_alloc(rs1 + 1, sizeof(Int));
    }
    Num->Slx = Slx;
    double *Rs = NULL;
    Int scale = Control->scale;  // if 1 the S will be scaled by max_row
    if (scale == 1) Rs = (double *)paru_calloc(Sym->m, sizeof(double));
    Num->Rs = Rs;
    if ((nf != 0 &&
         (rowMark == NULL || elRow == NULL || elCol == NULL ||
          rowSize == NULL || Work->lacList == NULL || RowList == NULL ||
          row_degree_bound == NULL || elementList == NULL ||
          Num->frowCount == NULL || Num->fcolCount == NULL ||
          Num->frowList == NULL || Num->fcolList == NULL ||
          Num->partial_Us == NULL || Num->partial_LUs == NULL ||
          Work->time_stamp == NULL || heapList == NULL ||
          (Sym->strategy == PARU_STRATEGY_SYMMETRIC &&
           (Diag_map == NULL || inv_Diag_map == NULL)))) ||

        // stuff that can be allocated even when nf==0
        Sx == NULL || (scale == 1 && Rs == NULL) ||
        (cs1 > 0 && (Sux == NULL || cSup == NULL)) ||
        (rs1 > 0 && (Slx == NULL || cSlp == NULL)) || cSp == NULL)
    {
        paru_free(m + 1, sizeof(Int), cSp);
        if (cs1 > 0) paru_free((cs1 + 1), sizeof(Int), cSup);
        if (rs1 > 0) paru_free((rs1 + 1), sizeof(Int), cSlp);
        return PARU_OUT_OF_MEMORY;
    }

    // Initializations
    if (nf != 0)
    {
        PRLEVEL(PR, ("%% $RowList =%p\n", RowList));
        paru_memset(rowSize, -1, m * sizeof(Int), Control);
        PRLEVEL(PR, ("%% rowSize pointer=%p size=%ld \n", rowSize,
                     m * sizeof(Int)));

        PRLEVEL(PR, ("%% rowMark pointer=%p size=%ld \n", rowMark,
                     (m + nf) * sizeof(Int)));

        paru_memset(elRow, -1, (m + nf) * sizeof(Int), Control);
        PRLEVEL(PR, ("%% elRow=%p\n", elRow));

        paru_memset(elCol, -1, (m + nf) * sizeof(Int), Control);
        PRLEVEL(PR, ("%% elCol=%p\n", elCol));

        PRLEVEL(PR, ("%% Work =%p\n ", Work));
    }

    //////////////////Initializing numerics Sx, Sux and Slx //////////////////{
    Int *Ap = (Int *)A->p;
    Int *Ai = (Int *)A->i;
    double *Ax = (double *)A->x;
    Int *Sp = Sym->Sp;
    Int *Slp = NULL;
    Int *Sup = NULL;
    paru_memcpy(cSp, Sp, (m + 1) * sizeof(Int), Control);
    if (cs1 > 0)
    {
        Sup = Sym->ustons.Sup;
        paru_memcpy(cSup, Sup, (cs1 + 1) * sizeof(Int), Control);
    }
    if (rs1 > 0)
    {
        Slp = Sym->lstons.Slp;
        paru_memcpy(cSlp, Slp, (rs1 + 1) * sizeof(Int), Control);
    }
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("Init Sup and Slp in the middle\n"));
    if (cs1 > 0)
    {
        PRLEVEL(PR, ("(%ld) Sup =", sunz));
        for (Int k = 0; k <= cs1; k++)
        {
            PRLEVEL(PR, ("%ld ", Sup[k]));
            PRLEVEL(PR + 2, ("c%ld ", cSup[k]));
            if (Sup[k] != cSup[k])
                PRLEVEL(PR, ("Sup[%ld] =%ld, cSup=%ld", k, Sup[k], cSup[k]));
            ASSERT(Sup[k] == cSup[k]);
        }
        PRLEVEL(PR, ("\n"));
    }
    if (rs1 > 0)
    {
        PRLEVEL(PR, ("(%ld) Slp =", slnz));
        for (Int k = 0; k <= rs1; k++)
        {
            PRLEVEL(PR, ("%ld ", Slp[k]));
            PRLEVEL(PR + 2, ("o%ld ", cSlp[k]));
            if (Slp[k] != cSlp[k])
                PRLEVEL(PR,
                        ("\nSup[%ld] =%ld, cSup=%ld\n", k, Slp[k], cSlp[k]));
            ASSERT(Slp[k] == cSlp[k]);
        }
        PRLEVEL(PR, ("\n"));
    }
#endif

    Int n1 = Sym->n1;

    Int *Qinit = Sym->Qfill;
    Int *Pinv = Sym->Pinv;
#ifndef NDEBUG
    PR = -1;
    PRLEVEL(PR, ("Iniit Pinv =\n"));
    for (Int i = 0; i < m; i++) PRLEVEL(PR, ("%ld ", Pinv[i]));
    PRLEVEL(PR, ("\n"));
#endif

    if (Rs)
    {
        for (Int newcol = 0; newcol < Sym->n; newcol++)
        {
            Int oldcol = Qinit[newcol];
            for (Int p = Ap[oldcol]; p < Ap[oldcol + 1]; p++)
            {
                Int oldrow = Ai[p];
                Rs[oldrow] = MAX(Rs[oldrow], fabs(Ax[p]));
            }
        }
    }

    PRLEVEL(PR, ("%% Rs:\n["));
    if (Rs)
    {  // making sure that every row has at most one element more than zero
        for (Int k = 0; k < m; k++)
        {
            PRLEVEL(PR, ("%lf ", Rs[k]));
            if (Rs[k] <= 0)
            {
                printf("Paru: Matrix is singular, row %ld is zero\n", k);
                Num->res = PARU_SINGULAR; 
                return PARU_SINGULAR;
            }
        }
    }
    PRLEVEL(PR, ("]\n"));

    for (Int newcol = 0; newcol < Sym->n; newcol++)
    {
        Int oldcol = Qinit[newcol];
        for (Int p = Ap[oldcol]; p < Ap[oldcol + 1]; p++)
        {
            Int oldrow = Ai[p];
            Int newrow = Pinv[oldrow];
            Int srow = newrow - n1;
            Int scol = newcol - n1;
            if (srow >= 0 && scol >= 0)
            {  // it is insdie S otherwise it is part of singleton
                Sx[cSp[srow]++] = (Rs == NULL) ? Ax[p] : Ax[p] / Rs[oldrow];
            }
            else if (srow < 0 && scol >= 0)
            {  // inside the U singletons
                PRLEVEL(PR, ("Usingleton newcol = %ld newrow=%ld\n", newcol,
                             newrow));
                // let the diagonal entries be first
                Sux[++cSup[newrow]] = (Rs == NULL) ? Ax[p] : Ax[p] / Rs[oldrow];
            }
            else
            {
                if (newrow < cs1)
                {  // inside U singletons CSR
                    // PRLEVEL(PR, ("Inside U singletons\n"));
                    if (newcol == newrow)
                    {  // diagonal entry
                        Sux[Sup[newrow]] =
                            (Rs == NULL) ? Ax[p] : Ax[p] / Rs[oldrow];
                    }
                    else
                    {
                        Sux[++cSup[newrow]] =
                            (Rs == NULL) ? Ax[p] : Ax[p] / Rs[oldrow];
                    }
                }
                else
                {  // inside L singletons CSC
                    // PRLEVEL(PR, ("Inside L singletons\n"));
                    if (newcol == newrow)
                    {  // diagonal entry
                        Slx[Slp[newcol - cs1]] =
                            (Rs == NULL) ? Ax[p] : Ax[p] / Rs[oldrow];
                    }
                    else
                    {
                        Slx[++cSlp[newcol - cs1]] =
                            (Rs == NULL) ? Ax[p] : Ax[p] / Rs[oldrow];
                    }
                }
            }
        }
    }

    if (Sym->cs1 > 0) paru_free((cs1 + 1), sizeof(Int), cSup);
    if (Sym->rs1 > 0) paru_free((rs1 + 1), sizeof(Int), cSlp);
        //////////////////Initializing numerics Sx, Sux and Slx
        /////////////////////}
#ifdef COUNT_FLOPS
    // flop count info init
    Work->flp_cnt_dgemm = 0.0;
    Work->flp_cnt_trsm = 0.0;
    Work->flp_cnt_dger = 0.0;
    Work->flp_cnt_real_dgemm = 0.0;
#endif

    // RowList, ColList and elementList are place holders
    // pointers to pointers that are allocated

    Int *Sj = Sym->Sj;

    /// ------------------------------------------------------------------------
    // create S = A (p,q)', or S=A(p,q), S is considered to be in row-form
    // -------------------------------------------------------------------------
#ifndef NDEBUG
    Int n = Num->n = Sym->n - Sym->n1;
    PRLEVEL(1, ("%% m=%ld, n=%ld\n", m, n));
    PR = 1;
    PRLEVEL(PR, ("\n%% Inside init row fronts\n"));
    PRLEVEL(PR, ("%% Sp =\n%%"));
    for (Int i = 0; i <= m; i++) PRLEVEL(PR, ("%ld ", Sp[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Sj =\n"));
    for (Int k = 0; k < snz; k++) PRLEVEL(PR, ("%ld ", Sj[k]));
    PRLEVEL(PR, ("\n"));
#endif

    PRLEVEL(1, ("InMatrix=[\n"));  // MATLAB matrix,

    // copying Diag_map
    if (Diag_map)
    {
#pragma omp taskloop default(none) shared(Sym, Diag_map, inv_Diag_map) \
    grainsize(512)
        for (Int i = 0; i < Sym->n; i++)
        {
            // paru_memcpy(Diag_map, Sym->Diag_map, (Sym->n) * sizeof(Int));
            Diag_map[i] = Sym->Diag_map[i];
            inv_Diag_map[Diag_map[i]] = i;
        }
#ifndef NDEBUG
        PR = 1;
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

    ParU_Ret info;
    Int out_of_memory = 0;
#pragma omp taskloop grainsize(512)
    for (Int row = 0; row < m; row++)
    {
        Int e = Sym->row2atree[row];
        Int nrows = 1,
            ncols =
                Sp[row + 1] - Sp[row];  // nrows and ncols of current front/row
        // printf("%% element %ld = %ld x %ld\n", e, nrows, ncols);

        row_degree_bound[row] = ncols;  // Initialzing row degree

        paru_element *curEl = elementList[e] =
            paru_create_element(nrows, ncols, 0);
        if (curEl == NULL)
        {  // out of memory
            printf("Paru: Out of memory: curEl\n");
#pragma omp atomic update
            out_of_memory += 1;
        }

        rowMark[e] = 0;

        // My new is calling paru_alloc
        std::vector<Int> *curHeap;
        try
        {
            curHeap = Work->heapList[e] = new std::vector<Int>;
        }
        catch (std::bad_alloc const &)
        {  // out of memory
            printf("Paru: Out of memory: curHeap\n");
#pragma omp atomic update
            out_of_memory += 1;
        }
        // printf("%%Heap allocated %p id=%ld \n", curHeap, e);

        curHeap->push_back(e);

#ifndef NDEBUG  // Printing the pointers info
                // printf ("%% curEl = %p ", curEl);
        // Int size = sizeof(paru_element) + sizeof(Int) * (2 * (nrows + ncols))
        //           + sizeof(double) * nrows * ncols;
        // printf("size= %ld\n", size);
#endif

        // constants for initialzing lists
        Int slackRow = 2;

        // Allocating Rowlist and updating its tuples
        RowList[row].list =
            (paru_tuple *)paru_alloc(slackRow * nrows, sizeof(paru_tuple));
        if (RowList[row].list == NULL)
        {  // out of memory
            printf("Paru: out of memory, RowList[row].list \n");
#pragma omp atomic update
            out_of_memory += 1;
        }
        RowList[row].numTuple = 0;
        RowList[row].len = slackRow;

        paru_tuple rowTuple;
        rowTuple.e = e;
        rowTuple.f = 0;
        if (paru_add_rowTuple(RowList, row, rowTuple))
        {
            printf("Paru: out of memory, add_rowTuple \n");
#pragma omp atomic update
            out_of_memory += 1;
        }

        // Allocating elements
        Int *el_colrowIndex = colIndex_pointer(curEl);
        double *el_colrowNum = numeric_pointer(curEl);

        // printf("el_colrowIndex =%p, el_colrowNum = %p \n",
        // el_colrowIndex, el_colrowNum);

        Int j = 0;  // Index inside an element
        for (Int p = Sp[row]; p < Sp[row + 1]; p++)
        {
            el_colrowIndex[j] = Sj[p];
            el_colrowNum[j++] = Sx[p];
        }
        el_colrowIndex[j++] = row;  // initializing element row index
        Work->lacList[e] = lac_el(elementList, e);
    }
    if (out_of_memory)
        info = PARU_OUT_OF_MEMORY;
    else
        info = PARU_SUCCESS;

    PRLEVEL(1, ("];\n"));
    PRLEVEL(1, ("I = InMatrix(:,1);\n"));
    PRLEVEL(1, ("J = InMatrix(:,2);\n"));
    PRLEVEL(1, ("X = InMatrix(:,3);\n"));
    PRLEVEL(1, ("S = sparse(I,J,X);\n"));

    return info;
}
