////////////////////////////////////////////////////////////////////////////////
////////////////////////// paru_mem.cpp ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  Wrappers for managing memory
 *  allocating and freeing is done through SuiteSparse and these wrappers
 *
 * @author Aznaveh
 *
 */
#include "paru_internal.hpp"
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
        printf("Paru: size must be > 0\n");
        return NULL;
    }
    else if (n >= (Size_max / size) || n >= INT_MAX)
    {
        // object is too big to allocate without causing integer overflow
        printf("Paru: problem too large\n");
        p = NULL;
    }
    else
    {
        p = SuiteSparse_malloc(n, size);
        if (p == NULL)
        {
            // out of memory
            printf("Paru: out of memory\n");
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
//  TODO consider making it parallel
void *paru_calloc(size_t n, size_t size)
{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int calloc_count = 0;
#endif
    void *p;
    if (size == 0)
    {
        printf("Paru: size must be > 0\n");
        return NULL;
    }
    else if (n >= (Size_max / size) || n >= INT_MAX)
    {
        // object is too big to allocate without causing integer overflow
        printf("Paru: problem too large\n");
        p = NULL;
    }
    else
    {
        p = SuiteSparse_calloc(n, size);
        if (p == NULL)
        {
            // out of memory
            printf("Paru: out of memory\n");
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
    if (*size == 0)
    {
        printf("Paru: size must be > 0\n");
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
                    oldP));
        p = oldP;
    }
    else if (newsize >= (Size_max / size_Entry) || newsize >= INT_MAX)
    {
        // object is too big to allocate without causing integer overflow
        printf("Paru: problem too large\n");
    }

    else
    {  // The object exists, and is changing to some other nonzero size.
        PRLEVEL(1, ("realloc : %ld to %ld, %ld\n", *size, newsize, size_Entry));
        int ok = TRUE;
        p = SuiteSparse_realloc(newsize, *size, size_Entry, oldP, &ok);
        if (ok)
        {
#ifndef NDEBUG
            realloc_count += newsize * size_Entry - *size;
            PRLEVEL(1, ("%% reallocated %ld in %p and freed %p total= %ld\n",
                        newsize * size_Entry, p, oldP, realloc_count));
#endif
            *size = newsize;
        }
    }
    return p;
}

//  Wrapper around free routine
//
void paru_free(size_t n, size_t size, void *p)
{
    DEBUGLEVEL(0);

    // static Int free_count = 0;
    // free_count += n * size;

    // Valgrind is unhappy about some part here
    //    PRLEVEL (1, ("%% free %ld in %p total= %ld\n",
    //                n*size, p, free_count));

    if (p != NULL)
        SuiteSparse_free(p);
    else
    {
        PRLEVEL(1, ("%% freeing a NULL pointer  \n"));
    }
}

//  Global replacement of new and delete
//
void *operator new(size_t size)
{  // no inline, required by [replacement.functions]/3
    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int cpp_count = 0;
    cpp_count += size;
    PRLEVEL(1, ("global op new called, size = %zu tot=%ld\n", size, cpp_count));
#endif

    if (size == 0)
        ++size;  // avoid malloc(0) which may return nullptr on success

    if (void *ptr = paru_alloc(1, size)) return ptr;
    throw std::bad_alloc{};
}
void operator delete(void *ptr) noexcept
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("global op delete called"));
    paru_free(0, 0, ptr);
}

//  freeing symbolic analysis data structure
ParU_Ret ParU_Freesym(ParU_Symbolic **Sym_handle, ParU_Control *Control)
{
    DEBUGLEVEL(0);
    if (Sym_handle == NULL || *Sym_handle == NULL)
        // nothing to do; caller probably ran out of memory
        return PARU_OUT_OF_MEMORY;

    ParU_Symbolic *Sym;
    Sym = *Sym_handle;

    Int m = Sym->m;
    Int n = Sym->n;
    Int n1 = Sym->n1;
    Int nf = Sym->nf;
    Int snz = Sym->snz;
    PRLEVEL(1, ("%% In free sym: m=%ld n=%ld\n nf=%ld "
                "Sym->anz=%ld \n",
                m, n, nf, Sym->anz));

    paru_free(nf + 1, sizeof(Int), Sym->Parent);
    paru_free(nf + 1, sizeof(Int), Sym->Child);
    paru_free(nf + 2, sizeof(Int), Sym->Childp);
    paru_free(nf + 1, sizeof(Int), Sym->Super);
    paru_free(nf, sizeof(Int), Sym->Depth);
    paru_free(n, sizeof(Int), Sym->Qfill);
    paru_free(n, sizeof(Int), Sym->Diag_map);
    paru_free(m, sizeof(Int), Sym->Ps);
    paru_free(m, sizeof(Int), Sym->Pfin);
    paru_free((m + 1), sizeof(Int), Sym->Pinit);
    paru_free(nf + 1, sizeof(Int), Sym->Fm);
    paru_free(nf + 1, sizeof(Int), Sym->Cm);

    // paru_free(Sym->num_roots, sizeof(Int), Sym->roots);

    paru_free(m + 1 - n1, sizeof(Int), Sym->Sp);
    paru_free(snz, sizeof(Int), Sym->Sj);
    paru_free(snz, sizeof(double), Sym->Sx);
    paru_free(n + 2 - n1, sizeof(Int), Sym->Sleft);

    // paru_free((n + 1), sizeof(Int), Sym->Chain_start);
    // paru_free((n + 1), sizeof(Int), Sym->Chain_maxrows);
    // paru_free((n + 1), sizeof(Int), Sym->Chain_maxcols);

    paru_free(nf + 1, sizeof(double), Sym->front_flop_bound);
    paru_free(nf + 1, sizeof(double), Sym->stree_flop_bound);

    Int ms = m - n1;  // submatrix is msxns

    paru_free(ms + nf, sizeof(Int), Sym->aParent);
    paru_free(ms + nf + 1, sizeof(Int), Sym->aChild);
    paru_free(ms + nf + 2, sizeof(Int), Sym->aChildp);
    paru_free(ms, sizeof(Int), Sym->row2atree);
    paru_free(nf, sizeof(Int), Sym->super2atree);
    paru_free(nf + 1, sizeof(Int), Sym->first);
    paru_free(m, sizeof(Int), Sym->scale_row);

    if (n1 > 0)
    {  // freeing singletons
        Int cs1 = Sym->cs1;
        if (cs1 > 0)
        {
            ParU_U_singleton ustons = Sym->ustons;
            paru_free(cs1 + 1, sizeof(Int), ustons.Sup);
            Int nnz = ustons.nnz;
            paru_free(nnz, sizeof(Int), ustons.Suj);
            paru_free(nnz, sizeof(Int), ustons.Sux);
        }

        Int rs1 = Sym->rs1;
        if (rs1 > 0)
        {
            ParU_L_singleton lstons = Sym->lstons;
            paru_free(rs1 + 1, sizeof(Int), lstons.Slp);
            Int nnz = lstons.nnz;
            paru_free(nnz, sizeof(Int), lstons.Sli);
            paru_free(nnz, sizeof(double), lstons.Slx);
        }
    }
    Int ntasks = Sym->ntasks;
    paru_free(ntasks + 1, sizeof(Int), Sym->task_map);
    paru_free(ntasks, sizeof(Int), Sym->task_parent);
    paru_free(ntasks, sizeof(Int), Sym->task_num_child);
    paru_free(ntasks, sizeof(Int), Sym->task_depth);

    paru_free(1, sizeof(ParU_Symbolic), Sym);

    *Sym_handle = NULL;
    return PARU_SUCCESS;
}

// free element e from elementList
void paru_free_el(Int e, ParU_Element **elementList)
{
    DEBUGLEVEL(0);
    ParU_Element *el = elementList[e];
    if (el == NULL) return;
    Int nrows = el->nrows, ncols = el->ncols;
    PRLEVEL(1, ("%%Free the element e =%ld\t", e));
    PRLEVEL(1, ("%% nrows =%ld ", nrows));
    PRLEVEL(1, ("%% ncols =%ld\n", ncols));
    Int tot_size = sizeof(ParU_Element) + sizeof(Int) * (2 * (nrows + ncols)) +
                   sizeof(double) * nrows * ncols;
    paru_free(1, tot_size, el);
    elementList[e] = NULL;
}

// It uses Sym, Do not free Sym before
ParU_Ret ParU_Freenum (ParU_Numeric **Num_handle, ParU_Control *Control)
{
    DEBUGLEVEL(0);
    if (Num_handle == NULL || *Num_handle == NULL) return PARU_INVALID;

    ParU_Numeric *Num;
    Num = *Num_handle;

    Int m = Num->m;  // m and n is different than Sym
    Int n = Num->n;  // Here there are submatrix size

    ParU_TupleList *RowList = Num->RowList;
    PRLEVEL(1, ("%% RowList =%p\n", RowList));

    ParU_Symbolic *Sym = Num->Sym;
    Int nf = Sym->nf;

    for (Int row = 0; row < m; row++)
    {
        Int len = RowList[row].len;
        paru_free(len, sizeof(ParU_Tuple), RowList[row].list);
    }
    paru_free(1, m * sizeof(ParU_TupleList), RowList);

    ParU_Element **elementList;
    elementList = Num->elementList;

    PRLEVEL(1, ("%% Sym = %p\n", Sym));
    PRLEVEL(1, ("%% freeing initialized elements:\n"));
    for (Int i = 0; i < m; i++)
    {  // freeing all row elements
        if (Sym == NULL)
        {
            printf("Paru: probably Sym has been freed before! Wrong usage\n");
            return PARU_INVALID;
        }
        Int e = Sym->row2atree[i];  // element number in augmented tree
        PRLEVEL(1, ("%% e =%ld\t", e));
        paru_free_el(e, elementList);
    }

    PRLEVEL(1, ("\n%% freeing CB elements:\n"));
    for (Int i = 0; i < nf; i++)
    {                                 // freeing all other elements
        Int e = Sym->super2atree[i];  // element number in augmented tree
        paru_free_el(e, elementList);
    }

    // free the answer
    ParU_Factors *LUs = Num->partial_LUs;
    ParU_Factors *Us = Num->partial_Us;

    for (Int i = 0; i < nf; i++)
    {
        paru_free(Num->frowCount[i], sizeof(Int), Num->frowList[i]);

        paru_free(Num->fcolCount[i], sizeof(Int), Num->fcolList[i]);

        PRLEVEL(1, ("%% Freeing Us=%p\n", Us[i].p));
        if (Us[i].p != NULL)
        {
            Int mm = Us[i].m;
            Int nn = Us[i].n;
            paru_free(mm * nn, sizeof(double), Us[i].p);
        }
        PRLEVEL(1, ("%% Freeing LUs=%p\n", LUs[i].p));
        if (LUs[i].p != NULL)
        {
            Int mm = LUs[i].m;
            Int nn = LUs[i].n;
            paru_free(mm * nn, sizeof(double), LUs[i].p);
        }
    }

    PRLEVEL(1, ("%% Done LUs\n"));
    paru_free(1, nf * sizeof(Int), Num->frowCount);
    paru_free(1, nf * sizeof(Int), Num->fcolCount);

    paru_free(1, nf * sizeof(Int *), Num->frowList);
    paru_free(1, nf * sizeof(Int *), Num->fcolList);

    paru_free(1, nf * sizeof(ParU_Factors), LUs);
    paru_free(1, nf * sizeof(ParU_Factors), Us);

    if (Num->Diag_map)
    {
        paru_free(n, sizeof(Int), Num->Diag_map);
        paru_free(n, sizeof(Int), Num->inv_Diag_map);
    }
#ifndef NDEBUG
    Int Us_bound_size = Sym->Us_bound_size;
    Int LUs_bound_size = Sym->LUs_bound_size;
    Int double_size = LUs_bound_size + Us_bound_size;
    Int row_Int_bound = Sym->row_Int_bound;
    Int col_Int_bound = Sym->col_Int_bound;
    Int int_size = row_Int_bound + col_Int_bound;
    Int upperBoundSize = double_size * sizeof(double) + int_size * sizeof(Int);
    PRLEVEL(1, ("%% FREE upperBoundSize =%ld \n", upperBoundSize));
#endif

    paru_free(1, nf * sizeof(Int), Num->time_stamp);
    // in practice each parent should deal with the memory for the children
#ifndef NDEBUG
    std::vector<Int> **heapList = Num->heapList;
    // freeing memory of heaps.
    if (heapList != NULL)
    {
        for (Int eli = 0; eli < m + nf + 1; eli++)
        {
            if (heapList[eli] != nullptr)
            {
                PRLEVEL(1,
                        ("%% %ld has not been freed %p\n", eli, heapList[eli]));
                delete heapList[eli];
                heapList[eli] = nullptr;
            }
            ASSERT(heapList[eli] == nullptr);
        }
    }
#endif
    paru_free(1, (m + nf + 1) * sizeof(std::vector<Int> **), Num->heapList);

    paru_free(1, (m + nf + 1) * sizeof(ParU_Element), elementList);

    Paru_Work *Work = Num->Work;
    if (Work != NULL)
    {
        paru_free(m, sizeof(Int), Work->rowSize);
        paru_free(m + nf + 1, sizeof(Int), Work->rowMark);
        paru_free(m + nf, sizeof(Int), Work->elRow);
        paru_free(m + nf, sizeof(Int), Work->elCol);
        paru_free(1, sizeof(Paru_Work), Num->Work);
    }

    paru_free(m + nf, sizeof(Int), Num->lacList);

    paru_free(m, sizeof(Int), Num->row_degree_bound);
    paru_free(1, sizeof(ParU_Numeric), Num);
    *Num_handle = NULL;
    return PARU_SUCCESS;
}
