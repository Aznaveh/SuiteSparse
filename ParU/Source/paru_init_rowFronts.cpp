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

ParU_Ret paru_init_rowFronts(ParU_Numeric **Num_handle,  // in/out
                                                         // inputs, not modified
                             cholmod_sparse *A,
                             // symbolic analysis
                             ParU_Symbolic *Sym,
                             ParU_Control *Control)
{
    // mallopt(M_TRIM_THRESHOLD, -1);         // disable sbrk trimming
    // mallopt(M_TOP_PAD, 16 * 1024 * 1024);  // increase padding to speedup
    // malloc

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

    Int m, n, nf;
    Num->Sym = Sym;
    m = Num->m = Sym->m - Sym->n1;
    n = Num->n = Sym->n - Sym->n1;
    nf = Sym->nf;
    Num->res = PARU_SUCCESS;
    Num->Control = Control;

    Paru_Work *Work = Num->Work = (Paru_Work *)paru_alloc(1, sizeof(Paru_Work));

    if (Work == NULL)
    {  // out of memory
        printf("Paru: out of memory: Work\n");
        ParU_Freenum(&Num, Control);
        return PARU_OUT_OF_MEMORY;
    }

    Work->rowMark = Work->elRow = NULL;
    Work->elCol = Work->rowSize = NULL;
    Num->row_degree_bound = NULL;
    Num->RowList = NULL;
    Num->lacList = NULL;
    Num->frowCount = NULL;
    Num->fcolCount = NULL;
    Num->frowList = NULL;
    Num->fcolList = NULL;
    Num->partial_Us = NULL;
    Num->partial_LUs = NULL;
    Num->heapList = NULL;
    Num->elementList = NULL;
    Num->time_stamp = NULL;
    Num->Diag_map = NULL;
    Num->inv_Diag_map = NULL;
    Num->Sx = NULL;
    Num->Sux = NULL;
    Num->Slx = NULL;
    Num->Rs= NULL;

    if (nf == 0)
    {  // nothing to be done
        return PARU_SUCCESS;
    }
    // Memory allocations for Num
    Int *rowMark = Work->rowMark = (Int *)paru_alloc(m + nf + 1, sizeof(Int));
    Int *elRow = Work->elRow = (Int *)paru_alloc(m + nf, sizeof(Int));
    Int *elCol = Work->elCol = (Int *)paru_alloc(m + nf, sizeof(Int));
    Int *rowSize = Work->rowSize = (Int *)paru_alloc(m, sizeof(Int));
    Int *row_degree_bound = Num->row_degree_bound =
        (Int *)paru_alloc(m, sizeof(Int));
    ParU_TupleList *RowList = Num->RowList =
        (ParU_TupleList *)paru_alloc(1, m * sizeof(ParU_TupleList));
    Num->lacList = (Int *)paru_alloc(m + nf, sizeof(Int));
    Num->frowCount = (Int *)paru_alloc(1, nf * sizeof(Int));
    Num->fcolCount = (Int *)paru_alloc(1, nf * sizeof(Int));
    Num->frowList = (Int **)paru_calloc(1, nf * sizeof(Int *));
    Num->fcolList = (Int **)paru_calloc(1, nf * sizeof(Int *));
    Num->partial_Us =  // Initialize with NULL
        (ParU_Factors *)paru_calloc(1, nf * sizeof(ParU_Factors));
    Num->partial_LUs =  // Initialize with NULL
        (ParU_Factors *)paru_calloc(1, nf * sizeof(ParU_Factors));
    Num->time_stamp = (Int *)paru_alloc(1, nf * sizeof(Int));

    std::vector<Int> **heapList = Num->heapList =
        (std::vector<Int> **)paru_calloc(
            1, (m + nf + 1) * sizeof(std::vector<Int> *));
    ParU_Element **elementList;
    elementList = Num->elementList =  // Initialize with NULL
        (ParU_Element **)paru_calloc(1, (m + nf + 1) * sizeof(ParU_Element));
    Int *Diag_map = Num->Diag_map = NULL;
    Int *inv_Diag_map = Num->inv_Diag_map = NULL;
    if (Sym->strategy == PARU_STRATEGY_SYMMETRIC)
    {
        Diag_map = Num->Diag_map = (Int *)paru_alloc(Sym->n, sizeof(Int));
        inv_Diag_map = Num->inv_Diag_map =
            (Int *)paru_alloc(Sym->n, sizeof(Int));
#ifndef NDEBUG
        paru_memset(Diag_map, 0, Sym->n * sizeof(Int), Control);
        paru_memset(inv_Diag_map, 0, Sym->n * sizeof(Int), Control);
#endif
    }

    Int snz = Sym->snz;
    double *SSx = Num->Sx = (double *)paru_alloc(snz, sizeof(double));
    Int *cSp = NULL; //copy of Sp, temporary for making Sx
    cSp = (Int *)paru_alloc(m + 1, sizeof(Int));
    double *Sux = NULL; 
    Int *cSup = NULL;  //copy of Sup temporary for making Sux
    Int cs1 = Sym->cs1; Int rs1 = Sym->rs1;
    if (cs1 > 0)
    {
        Int  sunz = Sym->ustons.nnz;
        Sux = (double *)paru_alloc(sunz, sizeof(double));
        cSup = (Int *)paru_alloc(cs1 + 1, sizeof(Int));
    }
    Num->Sux = Sux;
    double *Slx = NULL; 
    Int *cSlp = NULL; //copyf of Slp temporary, for making Slx
    if (rs1 > 0)
    {
        Int  slnz = Sym->lstons.nnz;
        Slx = (double *)paru_alloc(slnz, sizeof(double));
        cSlp = (Int *)paru_alloc(rs1 + 1, sizeof(Int));
    }
    Num->Slx = Slx;
    Int scale = Control->scale; // if 1 the S will be scaled by max_row
    double *Rs = NULL;
    if (scale == 1) Rs = (double *)paru_calloc(Sym->m, sizeof(double));

    if (rowMark == NULL || elRow == NULL || elCol == NULL || rowSize == NULL ||
            Num->lacList == NULL || RowList == NULL || row_degree_bound == NULL ||
            elementList == NULL || Num->frowCount == NULL ||
            Num->fcolCount == NULL || Num->frowList == NULL ||
            Num->fcolList == NULL || Num->partial_Us == NULL ||
            Num->partial_LUs == NULL || Num->time_stamp == NULL ||
            heapList == NULL || SSx == NULL || (scale == 1 && Rs == NULL) || 
            (cs1 > 0 && (Sux == NULL || cSup == NULL)) || 
            (rs1 > 0 && (Slx == NULL|| cSlp == NULL)) || cSup == NULL ||
            (Sym->strategy == PARU_STRATEGY_SYMMETRIC &&
             (Diag_map == NULL || inv_Diag_map == NULL)))
    {
        ParU_Freenum(&Num, Control);
        paru_free(m+1, sizeof(Int), cSp);
        if (cs1 > 0) paru_free((cs1 + 1), sizeof(Int), cSup);
        if (rs1 > 0) paru_free((rs1 + 1), sizeof(Int), cSlp);
        return PARU_OUT_OF_MEMORY;
    }

    // Initializations
    PRLEVEL(1, ("%% $RowList =%p\n", RowList));
    paru_memset(rowSize, -1, m * sizeof(Int), Control);
    PRLEVEL(1, ("%% rowSize pointer=%p size=%ld \n", rowSize, m * sizeof(Int)));

    PRLEVEL(1, ("%% rowMark pointer=%p size=%ld \n", rowMark,
                (m + nf) * sizeof(Int)));

    paru_memset(elRow, -1, (m + nf) * sizeof(Int), Control);
    PRLEVEL(1, ("%% elRow=%p\n", elRow));

    paru_memset(elCol, -1, (m + nf) * sizeof(Int), Control);
    PRLEVEL(1, ("%% elCol=%p\n", elCol));

    PRLEVEL(1, ("%% Work =%p\n ", Work));

    //////////////////Initializing numerics Sx, Sux and Slx //////////////////{
    Int *Sp = Sym->Sp;
    {
        Int *Ap = (Int *)A->p;
        Int *Ai = (Int *)A->i;
        double *Ax = (double *)A->x;

        paru_memcpy(cSp, Sp, (m + 1) * sizeof(Int), Control);
        if (cs1 > 0)
        {
            paru_memcpy(cSup, Sup, (cs1 + 1) * sizeof(Int), Control);
        }
        if (rs1 > 0)
        {
            paru_memcpy(cSlp, Slp, (rs1 + 1) * sizeof(Int), Control);
        }



        //for (Int newcol = n1; newcol < n; newcol++)
        //{
        //    Int oldcol = Qinit[newcol];
        //    for (Int p = Ap[oldcol]; p < Ap[oldcol + 1]; p++)
        //    {
        //        Int oldrow = Ai[p];
        //        Int newrow = Pinv[oldrow];
        //        Int srow = newrow - n1;
        //        Int scol = newcol - n1;
        //        if (srow >= 0)
        //        {  // it is insdie S otherwise it is part of singleton
        //            Sj[cSp[srow]] = scol;
        //            Sx[cSp[srow]++] = (Rs == NULL) ? Ax[p] : Ax[p] / Rs[oldrow];
        //        }
        //        else
        //        {  // inside the U singletons
        //            PRLEVEL(PR, ("Usingleton rest newcol = %ld newrow=%ld\n",
        //                        newcol, newrow));
        //            ASSERT(newrow != newcol);  // not a diagonal entry
        //            Suj[++cSup[newrow]] = newcol;
        //            Sux[cSup[newrow]] = (Rs == NULL) ? Ax[p] : Ax[p] / Rs[oldrow];
        //        }
        //    }
        //}

        paru_free(m+1, sizeof(Int), cSp);
        if (Sym->cs1 > 0) paru_free((cs1 + 1), sizeof(Int), cSup);
        if (Sym->rs1 > 0) paru_free((rs1 + 1), sizeof(Int), cSlp);
    }
    //////////////////Initializing numerics Sx, Sux and Slx //////////////////}
#ifdef COUNT_FLOPS
    // flop count info init
    Num->flp_cnt_dgemm = 0.0;
    Num->flp_cnt_trsm = 0.0;
    Num->flp_cnt_dger = 0.0;
    Num->flp_cnt_real_dgemm = 0.0;
#endif

    PRLEVEL(1, ("%% m=%ld, n=%ld\n", m, n));
    // RowList, ColList and elementList are place holders
    // pointers to pointers that are allocated
    if (m == 0 || n == 0)
    {
        printf("Paru: the dimension of matrix is zero: %ld x %ld \n", m, n);
        paru_free(1, sizeof(ParU_Numeric), Num);
        return PARU_INVALID;
    }

    double *Sx = Sym->Sx;
    Int *Sj = Sym->Sj;

    /// ------------------------------------------------------------------------
    // create S = A (p,q)', or S=A(p,q), S is considered to be in row-form
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
        #pragma omp taskloop default(none) shared(Sym, Diag_map, inv_Diag_map) \
         grainsize(512)
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

        ParU_Element *curEl = elementList[e] =
            paru_create_element(nrows, ncols, 0);
        if (curEl == NULL)
        {  // out of memory
            ParU_Freenum(&Num, Control);
            printf("Paru: Out of memory: curEl\n");
            #pragma omp atomic update
            out_of_memory += 1;
        }

        rowMark[e] = 0;

        // My new is calling paru_alloc
        std::vector<Int> *curHeap;
        try
        {
            curHeap = Num->heapList[e] = new std::vector<Int>;
        }
        catch (std::bad_alloc const &)
        {  // out of memory
            ParU_Freenum(&Num, Control);
            printf("Paru: Out of memory: curHeap\n");
            #pragma omp atomic update
            out_of_memory += 1;
        }
        // printf("%%Heap allocated %p id=%ld \n", curHeap, e);

        curHeap->push_back(e);

#ifndef NDEBUG  // Printing the pointers info
                // printf ("%% curEl = %p ", curEl);
        // Int size = sizeof(ParU_Element) + sizeof(Int) * (2 * (nrows + ncols))
        //           + sizeof(double) * nrows * ncols;
        // printf("size= %ld\n", size);
#endif

        // constants for initialzing lists
        Int slackRow = 2;

        // Allocating Rowlist and updating its tuples
        RowList[row].list =
            (ParU_Tuple *)paru_alloc(slackRow * nrows, sizeof(ParU_Tuple));
        if (RowList[row].list == NULL)
        {  // out of memory
            ParU_Freenum(&Num, Control);
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
            ParU_Freenum(&Num, Control);
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
            // printf("Sj[%ld] =%ld Sx[%ld]=%lf \n", p, Sj[p], p, Sx[p]);
            // for Matlab
            // printf("%ld,%ld, %.16lf;\n", row + 1, Sj[p] + 1, Sx[p]);
        }
        el_colrowIndex[j++] = row;  // initializing element row index
        Num->lacList[e] = lac_el(elementList, e);
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

    return info;
}
