////////////////////////////////////////////////////////////////////////////////
////////////////////////// paru_mem.cpp ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*! @brief  Wrappers for managing memory
 *  allocating and freeing is done through SuiteSparse and these wrappers
 *
 * @author Aznaveh
 *
 */

#ifdef PARU_ALLOC_TESTING
// global variables
bool paru_malloc_tracking = false;
Int paru_nmalloc = 0;

/* sample usage:

    To test:

        int info = paru_stuff ( ...) ;

    do this instead:

        paru_set_malloc_tracking (true) ;
        for (nmalloc = 0 ; ; nmalloc++)
        {
            paru_set_nmalloc (nmalloc) ;
            // do stuff
            int info = paru_stuff (...) ;
            if (info != PARU_OUT_OF_MEMORY) break ;
        }
        paru_set_malloc_tracking (false) ;

    or:

        int info ;
        BRUTAL_ALLOC_TEST (info, paru_sym (...)) ;
        if (info == PARU_SUCCESS)
        {
            BRUTAL_ALLOC_TEST (info, paru_num (...)) ;
        }
        if (info == PARU_SUCCESS)
        {
            BRUTAL_ALLOC_TEST (info, paru_solve (...)) ;
        }
        TEST_ASSERT (info == PARU_SUCCESS)
        free sym, num, soln

    // put in test coverage *.h file:
    #ifdef PARU_ALLOC_TESTING
    #define BRUTAL_ALLOC_TEST(info,method)              \
    {                                                   \
        paru_set_malloc_tracking (true) ;               \
        for (Int nmalloc = 0 ; ; Int nmalloc++)         \
        {                                               \
            paru_set_nmalloc (nmalloc)                  \
            // do stuff                                 \
            info = method ;                             \
            if (info != PARU_OUT_OF_MEMORY) break ;     \
            if (nmalloc > 1000000) { test failure }
        }                                               \
        paru_set_malloc_tracking (false) ;              \
    }
    #else
    #define BRUTAL_ALLOC_TEST(info,method)              \
    {                                                   \
        info = method ;                                 \
    }
    #endif
*/

bool paru_get_malloc_tracking(void)
{
    bool track;
#pragma omp critical paru_malloc_testing
    {
        track = paru_malloc_tracking;
    }
    return (track);
}

void paru_set_malloc_tracking(bool track)
{
#pragma omp critical paru_malloc_testing
    {
        paru_malloc_tracking = track;
    }
}

void paru_set_nmalloc(bool nmalloc)
{
#pragma omp critical paru_malloc_testing
    {
        paru_nmalloc = nmalloc;
    }
}

Int paru_decr_nmalloc(void)
{
    Int nmalloc = 0;
#pragma omp critical paru_malloc_testing
    {
        if (paru_nmalloc > 0)
        {
            nmalloc = paru_nmalloc--;
        }
    }
    return (nmalloc);
}

Int paru_get_nmalloc(void)
{
    Int nmalloc = 0;
#pragma omp critical paru_malloc_testing
    {
        nmalloc = paru_nmalloc;
    }
    return (nmalloc);
}

#endif

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
    void *p = NULL;
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
#ifdef PARU_ALLOC_TESTING
        // brutal memory testing only
        if (paru_get_malloc_tracking())
        {
            Int nmalloc = paru_decr_nmalloc();
            if (nmalloc > 0)
            {
                p = SuiteSparse_malloc(n, size);
            }
        }
        else
        {
            p = SuiteSparse_malloc(n, size);
        }
#else
        // in production
        p = SuiteSparse_malloc(n, size);
#endif

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
#ifdef PARU_ALLOC_TESTING
        // brutal memory testing only
        if (paru_get_malloc_tracking())
        {
            Int nmalloc = paru_decr_nmalloc();
            if (nmalloc > 0)
            {
                p = SuiteSparse_calloc(n, size);
            }
        }
        else
        {
            p = SuiteSparse_calloc(n, size);
        }
#else
        // in production
        p = SuiteSparse_calloc(n, size);
#endif

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
        p = paru_alloc(newsize, size_Entry);
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

#ifdef PARU_ALLOC_TESTING
        // brutal memory testing only
        if (paru_get_malloc_tracking())
        {
            Int nmalloc = paru_decr_nmalloc();
            if (nmalloc > 0)
            {
                p = SuiteSparse_realloc(newsize, *size, size_Entry, oldP, &ok);
            }
            else
            {
                // pretend to fail
                ok = FALSE;
            }
        }
        else
        {
            p = SuiteSparse_realloc(newsize, *size, size_Entry, oldP, &ok);
        }
#else
        // in production
        p = SuiteSparse_realloc(newsize, *size, size_Entry, oldP, &ok);
#endif

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
        return PARU_SUCCESS;

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
    paru_free(m, sizeof(Int), Sym->Pinv);

    if (n1 > 0)
    {  // freeing singletons
        Int cs1 = Sym->cs1;
        if (cs1 > 0)
        {
            ParU_U_singleton ustons = Sym->ustons;
            paru_free(cs1 + 1, sizeof(Int), ustons.Sup);
            Int nnz = ustons.nnz;
            paru_free(nnz, sizeof(Int), ustons.Suj);
        }

        Int rs1 = Sym->rs1;
        if (rs1 > 0)
        {
            ParU_L_singleton lstons = Sym->lstons;
            paru_free(rs1 + 1, sizeof(Int), lstons.Slp);
            Int nnz = lstons.nnz;
            paru_free(nnz, sizeof(Int), lstons.Sli);
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
void paru_free_el(Int e, paru_element **elementList)
{
    DEBUGLEVEL(0);
    paru_element *el = elementList[e];
    if (el == NULL) return;
#ifndef NDEBUG
    Int nrows = el->nrows, ncols = el->ncols;
    PRLEVEL(1, ("%%Free the element e =%ld\t", e));
    PRLEVEL(1, ("%% nrows =%ld ", nrows));
    PRLEVEL(1, ("%% ncols =%ld\n", ncols));
    Int tot_size = sizeof(paru_element) + sizeof(Int) * (2 * (nrows + ncols)) +
                   sizeof(double) * nrows * ncols;
    paru_free(1, tot_size, el);
#else
    paru_free(1, 0, el);
#endif
    elementList[e] = NULL;
}

ParU_Ret paru_free_work(ParU_Symbolic *Sym, paru_work *Work)
{
    Int m = Sym->m - Sym->n1;
    Int nf = Sym->nf;
    Int n = Sym->n - Sym->n1;
    paru_free(m, sizeof(Int), Work->rowSize);
    paru_free(m + nf + 1, sizeof(Int), Work->rowMark);
    paru_free(m + nf, sizeof(Int), Work->elRow);
    paru_free(m + nf, sizeof(Int), Work->elCol);

    paru_free(1, nf * sizeof(Int), Work->time_stamp);

    paru_tupleList *RowList = Work->RowList;
    PRLEVEL(1, ("%% RowList =%p\n", RowList));

    for (Int row = 0; row < m; row++)
    {
        Int len = RowList[row].len;
        paru_free(len, sizeof(paru_tuple), RowList[row].list);
    }
    paru_free(1, m * sizeof(paru_tupleList), RowList);

    if (Work->Diag_map)
    {
        paru_free(n, sizeof(Int), Work->Diag_map);
        paru_free(n, sizeof(Int), Work->inv_Diag_map);
    }

    paru_element **elementList;
    elementList = Work->elementList;

    PRLEVEL(1, ("%% Sym = %p\n", Sym));
    PRLEVEL(1, ("%% freeing initialized elements:\n"));
    for (Int i = 0; i < m; i++)
    {                               // freeing all row elements
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

    paru_free(1, (m + nf + 1) * sizeof(paru_element), elementList);

    paru_free(m + nf, sizeof(Int), Work->lacList);

    // in practice each parent should deal with the memory for the children
#ifndef NDEBUG
    std::vector<Int> **heapList = Work->heapList;
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
    paru_free(1, (m + nf + 1) * sizeof(std::vector<Int> **), Work->heapList);
    paru_free(m, sizeof(Int), Work->row_degree_bound);

    return PARU_SUCCESS;
}

ParU_Ret ParU_Freenum(ParU_Numeric **Num_handle,
                     ParU_Control *Control)
{
    DEBUGLEVEL(0);
    if (Num_handle == NULL || *Num_handle == NULL)
    {
        // nothing to do
        return PARU_SUCCESS;
    }

    ParU_Numeric *Num;
    Num = *Num_handle;

    Int m = Num->m;  
    Int nf = Num->nf;

    // freeing the numerical input
    paru_free(Num->snz, sizeof(double), Num->Sx);  
    if (Num->sunz > 0)                            
    {
        paru_free(Num->sunz, sizeof(double), Num->Sux);
    }
    if (Num->slnz > 0)                              
    {
        paru_free(Num->slnz, sizeof(double), Num->Slx); 
    }
 
    paru_free(Num->sym_m, sizeof(Int), Num->Rs);  // HERE: sym_m

    // free the factors
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

    paru_free(1, sizeof(ParU_Numeric), Num);
    *Num_handle = NULL;
    return PARU_SUCCESS;
}
