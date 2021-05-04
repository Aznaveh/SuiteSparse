/** =========================================================================  /
 * =======================  paru_mem.cpp ====================================  /
 * ========================================================================== */
/*! @brief  Wrappers for managing memory
 *  allocating and freeing is done through cholmod and these wrappers
 *
 * @author Aznaveh
 *
 */
#include "Parallel_LU.hpp"

void *paru_alloc(size_t n, size_t size, cholmod_common *cc)
{
    DEBUGLEVEL(0);
    //#ifndef NDEBUG
    static Int alloc_count = 0;
    alloc_count += n * size;
    //#endif
    void *p = cholmod_l_malloc(n, size, cc);
    PRLEVEL(1,
            ("%% allocated %ld in %p total= %ld\n", n * size, p, alloc_count));
    return p;
}

void *paru_calloc(size_t n, size_t size, cholmod_common *cc)
{
    DEBUGLEVEL(0);
    //#ifndef NDEBUG
    static Int calloc_count = 0;
    calloc_count += n * size;
    //#endif
    void *p = cholmod_l_calloc(n, size, cc);
    PRLEVEL(
        1, ("%% callocated %ld in %p total= %ld\n", n * size, p, calloc_count));
    return p;
}

void *paru_stack_calloc(size_t n, size_t size, paru_matrix *paruMatInfo,
                        cholmod_common *cc)
{
    DEBUGLEVEL(0);
    PRLEVEL(1,
            ("%% STACK ALLOC n= %ld size= %ld tot=%ld\n", n, size, n * size));
    static Int calloc_count = 0;
    calloc_count += n * size;

    if (paruMatInfo->stack_mem.mem_bank[0] == NULL)
    // first time
    {
        paru_symbolic *LUsym = paruMatInfo->LUsym;
        Int Us_bound_size = LUsym->Us_bound_size;
        Int LUs_bound_size = LUsym->LUs_bound_size;
        Int double_size = LUs_bound_size + Us_bound_size;
        Int row_Int_bound = LUsym->row_Int_bound;
        Int col_Int_bound = LUsym->col_Int_bound;
        Int int_size = row_Int_bound + col_Int_bound;
        size_t upperBoundSize =
            double_size * sizeof(double) + int_size * sizeof(Int);
        PRLEVEL(1, ("%% ALLOC upper = %ld\n", upperBoundSize));

        size_t size0;
        if (double_size < 1000)
            size0 = upperBoundSize;
        else
            // min
            size0 = (upperBoundSize < 64 * n * size) ? upperBoundSize
                                                     : 64 * n * size;

        PRLEVEL(1, ("%% ALLOC size0= %zu\n", size0));
        void *p = cholmod_l_calloc(size0, 1, cc);

        if (p == NULL)
        {
            printf("Memory allocation problem for the first bank alloc!\n");
            return NULL;
        }

        paruMatInfo->stack_mem.mem_bank[0] = p;
        paruMatInfo->stack_mem.avail = p;
        paruMatInfo->stack_mem.remaining = size0;
        paruMatInfo->stack_mem.size_bank[0] = size0;
        paruMatInfo->stack_mem.cur = 0;
    }

    Int remaining = paruMatInfo->stack_mem.remaining;
    remaining -= n * size;
    PRLEVEL(1, ("%% ALLOC remaining= %ld \n", remaining));

    if (remaining < 0)  // alloc memory for the next bank
    {
        Int cur = paruMatInfo->stack_mem.cur;
        PRLEVEL(1, ("%% ALLOC NEW cur=%ld", cur));
        if (cur == 63)
        {
            printf("LAST BANK ALREADY USED!!\n");
            return NULL;
        }
        size_t cur_size = paruMatInfo->stack_mem.size_bank[cur];

        // The only do while loop I have ever used
        do
        {
            cur_size *= 2;
            PRLEVEL(1, ("%% $$$cur_size= %zu \n", cur_size));
        } while (cur_size < n * size);
        void *p = cholmod_l_calloc(cur_size, 1, cc);
        if (p == NULL)
        {
            printf("Memory allocation problem for the %ld bank alloc!\n", cur);
            return NULL;
        }
        paruMatInfo->stack_mem.mem_bank[++cur] = p;
        paruMatInfo->stack_mem.cur = cur;
        paruMatInfo->stack_mem.size_bank[cur] = cur_size;
        paruMatInfo->stack_mem.avail = p;
        paruMatInfo->stack_mem.remaining = cur_size;
        remaining = cur_size - n * size;
    }

    PRLEVEL(1, ("%% ALLOC2 remaining= %ld \n", remaining));
    paruMatInfo->stack_mem.remaining = (size_t)remaining;
    void *avail = paruMatInfo->stack_mem.avail;
    // void *new_avail = (void*) ( (size_t *)avail + n*size + 1);
    void *new_avail = (void *)((char *)avail + n * size);
    paruMatInfo->stack_mem.avail = new_avail;

    PRLEVEL(1, ("%% callocated %ld in %p total= %ld\n", n * size, avail,
                calloc_count));
    return avail;
}
void *paru_realloc(
    Int newsize,     // requested size
    Int size_Entry,  // size of each Entry
    void *oldP,      // pointer to the old allocated space
    Int *size,       // a single number, input: old size, output: new size
    cholmod_common *cc)
{
    DEBUGLEVEL(0);
    //#ifndef NDEBUG
    static Int realloc_count = 0;
    realloc_count += newsize * size_Entry - *size;
    //#endif
    void *p = cholmod_l_realloc(newsize, size_Entry, oldP, (size_t *)size, cc);
    PRLEVEL(1, ("%% reallocated %ld in %p and freed %p total= %ld\n",
                newsize * size_Entry, p, oldP, realloc_count));
    return p;
}

void paru_free(Int n, Int size, void *p, cholmod_common *cc)
{
    DEBUGLEVEL(0);
    static Int free_count = 0;
    free_count += n * size;

    // Valgrind is unhappy about some part here
    //    PRLEVEL (1, ("%% free %ld in %p total= %ld\n",
    //                n*size, p, free_count));

    //#endif
    if (p != NULL)
        cholmod_l_free(n, size, p, cc);
    else
    {
        PRLEVEL(1, ("%% freeing a NULL pointer  \n"));
    }
}

void paru_freesym(paru_symbolic **LUsym_handle,
                  // workspace and parameters
                  cholmod_common *cc)
{
    DEBUGLEVEL(0);
    if (LUsym_handle == NULL || *LUsym_handle == NULL)
        // nothing to do; caller probably ran out of memory
        return;

    paru_symbolic *LUsym;
    LUsym = *LUsym_handle;

    Int m = LUsym->m;
    Int n = LUsym->n;
    Int n1 = LUsym->n1;
    Int nf = LUsym->nf;
    // Int anz = LUsym->anz;
    Int snz = LUsym->snz;
    Int rjsize = LUsym->rjsize;
    PRLEVEL(1, ("%% In free sym: m=%ld n=%ld\n nf=%ld "
                "LUsym->anz=%ld rjsize=%ld\n",
                m, n, nf, LUsym->anz, rjsize));

    paru_free(nf + 1, sizeof(Int), LUsym->Parent, cc);
    paru_free(nf + 1, sizeof(Int), LUsym->Child, cc);
    paru_free(nf + 2, sizeof(Int), LUsym->Childp, cc);
    paru_free(nf + 1, sizeof(Int), LUsym->Super, cc);
    paru_free(n, sizeof(Int), LUsym->Qfill, cc);
    paru_free(m, sizeof(Int), LUsym->Pinv, cc);
    paru_free((m + 1), sizeof(Int), LUsym->Pinit, cc);
    paru_free(nf + 1, sizeof(Int), LUsym->Fm, cc);
    paru_free(nf + 1, sizeof(Int), LUsym->Cm, cc);

    paru_free(rjsize, sizeof(Int), LUsym->Rj, cc);
    paru_free(nf + 1, sizeof(Int), LUsym->Rp, cc);

    paru_free(m + 1 - n1, sizeof(Int), LUsym->Sp, cc);
    paru_free(snz, sizeof(Int), LUsym->Sj, cc);
    paru_free(snz, sizeof(double), LUsym->Sx, cc);
    paru_free(n + 2 - n1, sizeof(Int), LUsym->Sleft, cc);

    paru_free((n + 1), sizeof(Int), LUsym->Chain_start, cc);
    paru_free((n + 1), sizeof(Int), LUsym->Chain_maxrows, cc);
    paru_free((n + 1), sizeof(Int), LUsym->Chain_maxcols, cc);

    Int ms = m - n1;  // submatrix is msxns

    paru_free(ms + nf, sizeof(Int), LUsym->aParent, cc);
    paru_free(ms + nf + 1, sizeof(Int), LUsym->aChild, cc);
    paru_free(ms + nf + 2, sizeof(Int), LUsym->aChildp, cc);
    paru_free(ms, sizeof(Int), LUsym->row2atree, cc);
    paru_free(nf, sizeof(Int), LUsym->super2atree, cc);
    paru_free(ms + nf, sizeof(Int), LUsym->first, cc);

    paru_free(1, sizeof(paru_symbolic), LUsym, cc);

    *LUsym_handle = NULL;
}

void paru_free_el(Int e, paru_Element **elementList, cholmod_common *cc)
/* fee element e from elementList */
{
    DEBUGLEVEL(0);
    paru_Element *el = elementList[e];
    if (el == NULL) return;
    Int nrows = el->nrows, ncols = el->ncols;
    PRLEVEL(1, ("%%Free the element e =%ld\t", e));
    PRLEVEL(1, ("%% nrows =%ld ", nrows));
    PRLEVEL(1, ("%% ncols =%ld\n", ncols));
    Int tot_size = sizeof(paru_Element) + sizeof(Int) * (2 * (nrows + ncols)) +
                   sizeof(double) * nrows * ncols;
    // TODO use GB_dealloc_memory and the free_pool
    // paru_free(1, tot_size, el, cc);
//  GB_free_memory ((void **) &el, tot_size /* FIXME: wrong size */) ;

    size_t size_allocated = el->size_allocated ;
    printf ("freeing element, tot_size %ld allocated %ld\n", tot_size,
        size_allocated) ;
    GB_dealloc_memory ((void **) &el, size_allocated) ;

    elementList[e] = NULL;
}

/*! It uses LUsym, Do not free LUsym before*/
void paru_freemat(paru_matrix **paruMatInfo_handle, cholmod_common *cc)
{
    DEBUGLEVEL(0);
    if (paruMatInfo_handle == NULL || *paruMatInfo_handle == NULL) return;

    paru_matrix *paruMatInfo;
    paruMatInfo = *paruMatInfo_handle;

    Int m = paruMatInfo->m;  // m and n is different than LUsym
    // Int n = paruMatInfo->n;       // Here there are submatrix size

    tupleList *RowList = paruMatInfo->RowList;
    PRLEVEL(1, ("%% RowList =%p\n", RowList));

    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int nf = LUsym->nf;

    for (Int row = 0; row < m; row++)
    {
        Int len = RowList[row].len;
        paru_free(len, sizeof(paru_Tuple), RowList[row].list, cc);
    }
    paru_free(1, m * sizeof(tupleList), RowList, cc);

    paru_Element **elementList;
    elementList = paruMatInfo->elementList;

    PRLEVEL(1, ("%% LUsym = %p\n", LUsym));
    PRLEVEL(1, ("%% freeing initialized elements:\n"));
    for (Int i = 0; i < m; i++)
    {  // freeing all row elements
        if (LUsym == NULL)
        {
            printf("Probably LUsym has been freed before! Wrong usage\n");
            return;
        }
        Int e = LUsym->row2atree[i];  // element number in augmented tree
        PRLEVEL(1, ("%% e =%ld\t", e));
        paru_free_el(e, elementList, cc);
    }

    PRLEVEL(1, ("\n%% freeing CB elements:\n"));
    for (Int i = 0; i < nf; i++)
    {                                   // freeing all other elements
        Int e = LUsym->super2atree[i];  // element number in augmented tree
        paru_free_el(e, elementList, cc);
    }

    // free the answer
    paru_fac *LUs = paruMatInfo->partial_LUs;
    paru_fac *Us = paruMatInfo->partial_Us;

    PRLEVEL(1, ("%% Done LUs\n"));
    paru_free(1, nf * sizeof(Int), paruMatInfo->frowCount, cc);
    paru_free(1, nf * sizeof(Int), paruMatInfo->fcolCount, cc);

    paru_free(1, nf * sizeof(Int *), paruMatInfo->frowList, cc);
    paru_free(1, nf * sizeof(Int *), paruMatInfo->fcolList, cc);

    paru_free(1, nf * sizeof(paru_fac), LUs, cc);
    paru_free(1, nf * sizeof(paru_fac), Us, cc);

#ifndef NDEBUG
    Int Us_bound_size = LUsym->Us_bound_size;
    Int LUs_bound_size = LUsym->LUs_bound_size;
    Int double_size = LUs_bound_size + Us_bound_size;
    Int row_Int_bound = LUsym->row_Int_bound;
    Int col_Int_bound = LUsym->col_Int_bound;
    Int int_size = row_Int_bound + col_Int_bound;
    Int upperBoundSize = double_size * sizeof(double) + int_size * sizeof(Int);
    PRLEVEL(1, ("%% FREE upperBoundSize =%ld \n", upperBoundSize));
#endif

    for (Int i = 0; i < 64; i++)
    {
        if (paruMatInfo->stack_mem.mem_bank[i] == NULL) break;
        paru_free(paruMatInfo->stack_mem.size_bank[i], 1,
                  paruMatInfo->stack_mem.mem_bank[i], cc);
    }

    paru_free(1, nf * sizeof(Int), paruMatInfo->time_stamp, cc);
    // in practice each parent should deal with the memory for the children
#ifndef NDEBUG
    std::vector<Int> **heapList = paruMatInfo->heapList;
    // freeing memory of heaps.
    for (Int eli = 0; eli < m + nf + 1; eli++)
    {
        if (heapList[eli] != nullptr)
        {
            PRLEVEL(1, ("%% %ld has not been freed %p\n", eli, heapList[eli]));
            delete heapList[eli];
            heapList[eli] = nullptr;
        }
        ASSERT(heapList[eli] == nullptr);
    }
#endif
    paru_free(1, (m + nf + 1) * sizeof(std::vector<Int> **),
              paruMatInfo->heapList, cc);

    paru_free(1, (m + nf + 1) * sizeof(paru_Element), elementList, cc);
    work_struct *Work = paruMatInfo->Work;
    paru_free(m, sizeof(Int), Work->rowSize, cc);
    paru_free(m + nf + 1, sizeof(Int), Work->rowMark, cc);
    paru_free(m + nf, sizeof(Int), Work->elRow, cc);
    paru_free(m + nf, sizeof(Int), Work->elCol, cc);

    paru_free(m + nf, sizeof(Int), paruMatInfo->lacList, cc);

    paru_free(m, sizeof(Int), paruMatInfo->scale_row, cc);
    paru_free(m, sizeof(Int), paruMatInfo->row_degree_bound, cc);
    paru_free(1, sizeof(work_struct), paruMatInfo->Work, cc);
    paru_free(1, sizeof(paru_matrix), paruMatInfo, cc);
    *paruMatInfo_handle = NULL;
    GB_free_pool_finalize ( ) ;
}
