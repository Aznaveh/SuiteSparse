/** =========================================================================  /
 * =======================  paru_mem.cpp ====================================  /
 * ========================================================================== */
/*! @brief  Wrappers for managing memory
 *  allocating and freeing is done through SuiteSparse and these wrappers
 *
 * @author Aznaveh
 *
 */
#include "Parallel_LU.hpp"
//  Wrapper around malloc routine
//
//  Uses a pointer to the malloc routine.
void *paru_alloc(size_t n, size_t size)
{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int alloc_count = 0;
#endif
    void *p;
    if (size == 0)
    {
        printf("size must be > 0\n");
        return NULL;
    }
    else if (n >= (Size_max / size) || n >= INT_MAX)
    {
        // object is too big to allocate without causing integer overflow
        printf("problem too large\n");
        p = NULL;
    }
    else
    {
        p = SuiteSparse_malloc(n, size);
        if (p == NULL)
        {
            // out of memory
            printf("out of memory\n");
        }
        else
        {
#ifndef NDEBUG
            PRLEVEL(1, ("%% allocated %ld in %p total= %ld\n", n * size, p,
                        alloc_count));
            alloc_count += n * size;
#endif
        }
    }
    return p;
}

//  Wrapper around calloc routine
//
//  Uses a pointer to the calloc routine.
void *paru_calloc(size_t n, size_t size)
{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int calloc_count = 0;
#endif
    void *p;
    if (size == 0)
    {
        printf("size must be > 0\n");
        return NULL;
    }
    else if (n >= (Size_max / size) || n >= INT_MAX)
    {
        // object is too big to allocate without causing integer overflow
        printf("problem too large\n");
        p = NULL;
    }
    else
    {
        p = SuiteSparse_calloc(n, size);
        if (p == NULL)
        {
            // out of memory
            printf("out of memory\n");
        }
        else
        {
#ifndef NDEBUG
            PRLEVEL(1, ("%% callocated %ld in %p total= %ld\n", n * size, p,
                        calloc_count));
            calloc_count += n * size;
#endif
        }
    }
    return p;
}

//  Wrapper around realloc routine
//
//  Uses a pointer to the realloc routine.
void *paru_realloc(
    size_t newsize,     // requested size
    size_t size_Entry,  // size of each Entry
    void *oldP,         // pointer to the old allocated space
    size_t *size)       // a single number, input: old size, output: new size
{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int realloc_count = 0;
#endif
    void *p = NULL;
    if (size == 0)
    {
        printf("size must be > 0\n");
        return NULL;
    }
    else if (oldP == NULL)
    {  // A new alloc
        p = SuiteSparse_malloc(newsize, size_Entry);
        *size = (p == NULL) ? 0 : newsize * size_Entry;
    }
    else if (newsize == *size)
    {
        PRLEVEL(1, ("%% reallocating nothing %ld, %ld in %p \n", newsize, *size,
                    oldP, ));
    }
    else if (newsize >= (Size_max / size_Entry) || newsize >= INT_MAX)
    {
        // object is too big to allocate without causing integer overflow
        printf("problem too large\n");
    }

    else
    {  // The object exists, and is changing to some other nonzero size.
        PRLEVEL(1, ("realloc : %d to %d, %d\n", *size, newsize, size_Entry));
        int ok = TRUE;
        p = SuiteSparse_realloc(newsize, *size, size_Entry, oldP, &ok);
        if (ok)
        {
#ifndef NDEBUG
            realloc_count += newsize * size_Entry - *size;
#endif
            PRLEVEL(1, ("%% reallocated %ld in %p and freed %p total= %ld\n",
                        newsize * size_Entry, p, oldP, realloc_count));
            *size = newsize;
        }
    }
    return p;
}

//  Wrapper around free routine
//
void paru_free(Int n, Int size, void *p)
{
    DEBUGLEVEL(0);
    static Int free_count = 0;
    free_count += n * size;

    // Valgrind is unhappy about some part here
    //    PRLEVEL (1, ("%% free %ld in %p total= %ld\n",
    //                n*size, p, free_count));

    //#endif
    if (p != NULL)
        SuiteSparse_free(p);
    else
    {
        PRLEVEL(1, ("%% freeing a NULL pointer  \n"));
    }
}

void paru_freesym(paru_symbolic **LUsym_handle)
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

    paru_free(nf + 1, sizeof(Int), LUsym->Parent);
    paru_free(nf + 1, sizeof(Int), LUsym->Child);
    paru_free(nf + 2, sizeof(Int), LUsym->Childp);
    paru_free(nf + 1, sizeof(Int), LUsym->Super);
    paru_free(n, sizeof(Int), LUsym->Qfill);
    paru_free(m, sizeof(Int), LUsym->Pinv);
    paru_free((m + 1), sizeof(Int), LUsym->Pinit);
    paru_free(nf + 1, sizeof(Int), LUsym->Fm);
    paru_free(nf + 1, sizeof(Int), LUsym->Cm);

    paru_free(rjsize, sizeof(Int), LUsym->Rj);
    paru_free(nf + 1, sizeof(Int), LUsym->Rp);

    paru_free(m + 1 - n1, sizeof(Int), LUsym->Sp);
    paru_free(snz, sizeof(Int), LUsym->Sj);
    paru_free(snz, sizeof(double), LUsym->Sx);
    paru_free(n + 2 - n1, sizeof(Int), LUsym->Sleft);

    paru_free((n + 1), sizeof(Int), LUsym->Chain_start);
    paru_free((n + 1), sizeof(Int), LUsym->Chain_maxrows);
    paru_free((n + 1), sizeof(Int), LUsym->Chain_maxcols);

    Int ms = m - n1;  // submatrix is msxns

    paru_free(ms + nf, sizeof(Int), LUsym->aParent);
    paru_free(ms + nf + 1, sizeof(Int), LUsym->aChild);
    paru_free(ms + nf + 2, sizeof(Int), LUsym->aChildp);
    paru_free(ms, sizeof(Int), LUsym->row2atree);
    paru_free(nf, sizeof(Int), LUsym->super2atree);
    paru_free(ms + nf, sizeof(Int), LUsym->first);

    paru_free(1, sizeof(paru_symbolic), LUsym);

    *LUsym_handle = NULL;
}

void paru_free_el(Int e, paru_Element **elementList)
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
    paru_free(1, tot_size, el);
    elementList[e] = NULL;
}

/*! It uses LUsym, Do not free LUsym before*/
void paru_freemat(paru_matrix **paruMatInfo_handle)
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
        paru_free(len, sizeof(paru_Tuple), RowList[row].list);
    }
    paru_free(1, m * sizeof(tupleList), RowList);

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
        paru_free_el(e, elementList);
    }

    PRLEVEL(1, ("\n%% freeing CB elements:\n"));
    for (Int i = 0; i < nf; i++)
    {                                   // freeing all other elements
        Int e = LUsym->super2atree[i];  // element number in augmented tree
        paru_free_el(e, elementList);
    }

    // free the answer
    paru_fac *LUs = paruMatInfo->partial_LUs;
    paru_fac *Us = paruMatInfo->partial_Us;

    for (Int i = 0; i < nf; i++)
    {
        paru_free(paruMatInfo->frowCount[i], sizeof(Int),
                  paruMatInfo->frowList[i]);

        paru_free(paruMatInfo->fcolCount[i], sizeof(Int),
                  paruMatInfo->fcolList[i]);

        PRLEVEL(1, ("%% Freeing Us=%p\n", Us[i].p));
        if (Us[i].p != NULL)
        {
            Int m = Us[i].m;
            Int n = Us[i].n;
            paru_free(m * n, sizeof(double), Us[i].p);
        }
        PRLEVEL(1, ("%% Freeing LUs=%p\n", LUs[i].p));
        if (LUs[i].p != NULL)
        {
            Int m = LUs[i].m;
            Int n = LUs[i].n;
            paru_free(m * n, sizeof(double), LUs[i].p);
        }
    }

    PRLEVEL(1, ("%% Done LUs\n"));
    paru_free(1, nf * sizeof(Int), paruMatInfo->frowCount);
    paru_free(1, nf * sizeof(Int), paruMatInfo->fcolCount);

    paru_free(1, nf * sizeof(Int *), paruMatInfo->frowList);
    paru_free(1, nf * sizeof(Int *), paruMatInfo->fcolList);

    paru_free(1, nf * sizeof(paru_fac), LUs);
    paru_free(1, nf * sizeof(paru_fac), Us);

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

    paru_free(1, nf * sizeof(Int), paruMatInfo->time_stamp);
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
              paruMatInfo->heapList);

    paru_free(1, (m + nf + 1) * sizeof(paru_Element), elementList);
    work_struct *Work = paruMatInfo->Work;
    paru_free(m, sizeof(Int), Work->rowSize);
    paru_free(m + nf + 1, sizeof(Int), Work->rowMark);
    paru_free(m + nf, sizeof(Int), Work->elRow);
    paru_free(m + nf, sizeof(Int), Work->elCol);

    paru_free(m + nf, sizeof(Int), paruMatInfo->lacList);

    paru_free(m, sizeof(Int), paruMatInfo->scale_row);
    paru_free(m, sizeof(Int), paruMatInfo->row_degree_bound);
    paru_free(1, sizeof(work_struct), paruMatInfo->Work);
    paru_free(1, sizeof(paru_matrix), paruMatInfo);
    *paruMatInfo_handle = NULL;
}
