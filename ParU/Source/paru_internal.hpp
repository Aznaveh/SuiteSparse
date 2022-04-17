////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_internal.hpp //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

#ifndef PARU_INTERNAL_H
#define PARU_INTERNAL_H
//!
//  internal libraries that are not visible to the user
//  @author Aznaveh
//
#include "ParU.hpp"
#include "paru_omp.hpp"

// -----------------------------------------------------------------------------
// debugging and printing macros
// -----------------------------------------------------------------------------
// force debugging off
#ifndef NDEBUG
#define NDEBUG
#endif

#ifndef NPR
#define NPR
#endif

#ifndef NTIME
#define NTIME
#endif

#define DLONG
#ifdef Int  // Int defined in amd as int without DLONG
#undef Int  // and I need inside of amd because of UMFPACK
#endif

#ifdef PARU_ALLOC_TESTING 
#undef PARU_ALLOC_TESTING 
#endif

extern "C"
{
#include "cholmod_blas.h"
#include "umf_internal.h"
}

// for printing information uncomment this; to activate assertions uncomment
//#undef NPR  //<<1>>
// uncomment the following line to turn on debugging mode
//#undef NDEBUG  //<<2>>
// uncomment the following line to turn on OpenMP timing
//#undef NTIME   //<<3>>
#define PARU_ALLOC_TESTING // for coverage test allocations

// uncomment if you want to count hardware flops
//#define COUNT_FLOPS

// defined somewhere else
#ifdef ASSERT
#undef ASSERT
#endif
#ifndef NDEBUG
#include <assert.h>
#define ASSERT(e) assert(e)
#else
#define ASSERT(e)
#endif

#ifndef NPR
static int print_level = 0;
#define PRLEVEL(level, param)                   \
    {                                           \
        if (print_level >= level) printf param; \
    }
#define DEBUGLEVEL(level)    \
    {                        \
        print_level = level; \
    }
#define PARU_DEFINE_PRLEVEL int PR = 1
#else
#define PRLEVEL(level, param)
#define DEBUGLEVEL(level)
#define PARU_DEFINE_PRLEVEL
#endif

// These libraries are included in Suitesparse_config
//#include <stdlib.h>
//#include <math.h>
//#include <float.h>
//#include <stdio.h>
//#include <cstring>
//#include <malloc.h> // mallopt used in paru_init_rowFronts.cpp

// To be able to use set
#include <algorithm>
#include <set>
#include <vector>


// -----------------------------------------------------------------------------
// basic macros
// -----------------------------------------------------------------------------

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define EMPTY (-1)
// already defined in amd
// define TRUE 1
// define FALSE 0
#define IMPLIES(p, q) (!(p) || (q))

// NULL should already be defined, but ensure it is here.
#ifndef NULL
#define NULL ((void *)0)
#endif

#define Size_max ((size_t)(-1))  // the largest value of size_t

// internal data structures
struct heaps_info
{
    Int sum_size, biggest_Child_size, biggest_Child_id;
};

// =============================================================================
//                   ParU_Tuple, Row data structure
// =============================================================================
struct paru_tuple
{
    // The (e,f) tuples for element lists
    Int e,  //  element number
        f;  //   offest
};

struct paru_tupleList
{                      // List of tuples
    Int numTuple,      //  number of Tuples in this element
        len;           //  length of allocated space for current list
    paru_tuple *list;  // list of tuples regarding to this element
};

// =============================================================================
//                  An element, contribution block
// =============================================================================
struct paru_element
{
    Int

        nrowsleft,  // number of rows remaining
        ncolsleft,  // number of columns remaining
        nrows,      // number of rows
        ncols,      // number of columns
        rValid,     // validity of relative row index
        cValid;     // validity of relative column index

    Int lac;  // least active column which is active
    // 0 <= lac <= ncols

    Int nzr_pc;  // number of zero rows in pivotal column of current front

    size_t size_allocated;
    // followed in memory by:
    //   Int
    //   col [0..ncols-1],  column indices of this element
    //   row [0..nrows-1] ; row indices of this element
    //
    //   relColInd [0..ncols-1];    relative indices of this element for
    //   current front
    //   relRowInd [0..nrows-1],    relative indices of this element for
    //   current front
    //   double ncols*nrows; numeric values
};

struct paru_work
{
    // gather scatter space for rows
    Int *rowSize;  // Initalized data structure, size of rows
    // Int rowMark;      // Work->rowSize[x] < rowMark[eli] for each front
    Int *rowMark;  // size = m+nf

    // gather scatter space for elements
    Int *elRow;  // Initalized data structure, size m+nf
    Int *elCol;  // Initalized data structure, size m+nf

    // only used for statistics when debugging is enabled:
    Int actual_alloc_LUs;      // actual memory allocated for LUs
    Int actual_alloc_Us;       // actual memory allocated for Us
    Int actual_alloc_row_int;  // actual memory allocated for rows
    Int actual_alloc_col_int;  // actual memory allocated for cols

#ifdef COUNT_FLOPS
    double flp_cnt_dgemm;
    double flp_cnt_trsm;
    double flp_cnt_dger;
    double flp_cnt_real_dgemm;
#endif

    paru_tupleList *RowList;  // size n of dynamic list
    Int *time_stamp;          // for relative index update; not initialized

    Int *Diag_map;  // size n,
    // Both of these are NULL if the stratey is not symmetric
    // copy of Diag_map from Sym;
    // this copy can be updated during the factorization
    Int *inv_Diag_map;  // size n,
    // inverse of Diag_map from Sym;
    // It helps editing the Diag_map

    Int *row_degree_bound;  // row degree size number of rows

    ParU_Symbolic *Sym;  // point to the symbolic that user sends

    paru_element **elementList;  // pointers to all elements, size = m+nf+1

    Int *lacList;  // size m+nf least active column of each element
                   //    el_colIndex[el->lac]  == lacList [e]
                   //    number of element
                   //
    // each active front owns and manage a heap list. The heap is based on the
    // least numbered column. The active front Takes the pointer of the biggest
    // child and release its other children after concatenating their list to
    // its own. The list of heaps are initialized by nullptr
    std::vector<Int> **heapList;  // size m+nf+1, initialized with nullptr

    Int naft;  // number of actvie frontal tasks
    Int resq;  // number of remainig ready tasks in the queue
};

//------------------------------------------------------------------------------
// inline internal functions

inline Int *colIndex_pointer(paru_element *curEl) { return (Int *)(curEl + 1); }
// Never ever use these functions prior to initializing ncols and nrows
inline Int *rowIndex_pointer(paru_element *curEl)
{
    return (Int *)(curEl + 1) + curEl->ncols;
}

inline Int *relColInd(paru_element *curEl)
//{    return (Int*)(curEl+1) + curEl->ncols + curEl->nrows + 1;}
{
    return (Int *)(curEl + 1) + curEl->ncols + curEl->nrows;
}

inline Int *relRowInd(paru_element *curEl)
//{    return (Int*)(curEl+1) + 2*curEl->ncols + curEl->nrows + 2;}
{
    return (Int *)(curEl + 1) + 2 * curEl->ncols + curEl->nrows;
}

inline double *numeric_pointer(paru_element *curEl)
// sizeof Int and double are same, but I keep it like this for clarity
//{ return (double*)((Int*)(curEl+1) + 2*curEl->ncols + 2*curEl->nrows + 2);}
{
    return (double *)((Int *)(curEl + 1) + 2 * curEl->ncols + 2 * curEl->nrows);
}

inline Int flip(Int colInd) { return -colInd - 2; }

inline Int lac_el(paru_element **elementList, Int eli)
{  // return least numbered column of the element i (eli)
    if (elementList[eli] == NULL)
        return LONG_MAX;
    else
    {
        Int *el_colIndex = (Int *)(elementList[eli] + 1);
        Int lac_ind = elementList[eli]->lac;
        return el_colIndex[lac_ind];
    }
}

//------------------------------------------------------------------------------
// internal routines
//
/* Wrappers for managing memory */
void *paru_alloc(size_t n, size_t size);
void *paru_calloc(size_t n, size_t size);
void *paru_realloc(size_t newsize, size_t size_Entry, void *oldP, size_t *size);

void paru_free(size_t n, size_t size, void *p);
void paru_free_el(Int e, paru_element **elementList);

void *operator new(std::size_t sz);
void operator delete(void *ptr) noexcept;

void paru_memset(void *ptr, Int value, size_t num, ParU_Control *Control);
void paru_memcpy(void *destination, const void *source, size_t num,
                 ParU_Control *Control);

#ifdef PARU_ALLOC_TESTING
bool paru_get_malloc_tracking (void) ;
void paru_set_malloc_tracking (bool track) ;
void paru_set_nmalloc (bool nmalloc) ;
Int paru_decr_nmalloc (void) ;
Int paru_get_nmalloc (void) ;
#endif

/* add tuple functions defintions */
Int paru_add_rowTuple(paru_tupleList *RowList, Int row, paru_tuple T);

Int paru_dgetrf(double *F, Int *frowList, Int m, Int n, BLAS_INT *ipiv);

Int paru_factorize_full_summed(Int f, Int start_fac,
                               std::vector<Int> &panel_row,
                               std::set<Int> &stl_colSet,
                               std::vector<Int> &pivotal_elements,
                               paru_work *Work, ParU_Numeric *Num);

paru_element *paru_create_element(Int nrows, Int ncols, Int init);

void paru_assemble_row_2U(Int e, Int f, Int sR, Int dR,
                          std::vector<Int> &colHash, paru_work *Work,
                          ParU_Numeric *Num);

Int paru_trsm(Int f, double *pF, double *uPart, Int fp, Int rowCount,
              Int colCount, paru_work *Work, ParU_Numeric *Num);
Int paru_dgemm(Int f, double *pF, double *uPart, double *el, Int fp,
               Int rowCount, Int colCount, paru_work *Work, ParU_Numeric *Num);

void paru_print_element(Int e, paru_work *Work, ParU_Numeric *Num);
void paru_print_paru_tupleList(paru_tupleList *listSet, Int index);
void paru_init_rel(Int f, paru_work *Work);

void paru_update_rel_ind_col(Int e, Int f, std::vector<Int> &colHash,
                             paru_work *Work, ParU_Numeric *Num);

void paru_update_rowDeg(Int panel_num, Int row_end, Int f, Int start_fac,
                        std::set<Int> &stl_colSet,
                        std::vector<Int> &pivotal_elements, paru_work *Work,
                        ParU_Numeric *Num);

Int paru_cumsum(Int n, Int *X, ParU_Control *Control);

Int bin_srch_col(Int *srt_lst, Int l, Int r, Int num);
Int bin_srch(Int *srt_lst, Int l, Int r, Int num);

ParU_Ret paru_init_rowFronts(paru_work *Work, ParU_Numeric **Num_handle,
                             cholmod_sparse *A, ParU_Symbolic *Sym,
                             ParU_Control *Control);
ParU_Ret paru_front(Int f, paru_work *Work, ParU_Numeric *Num);

ParU_Ret paru_pivotal(std::vector<Int> &pivotal_elements,
                      std::vector<Int> &panel_row, Int &zero_piv_rows, Int f,
                      heaps_info &hi, paru_work *Work, ParU_Numeric *Num);

int paru_intersection(Int e, paru_element **elementList,
                      std::set<Int> &stl_colSet);

ParU_Ret paru_prior_assemble(Int f, Int start_fac,
                             std::vector<Int> &pivotal_elements,
                             std::vector<Int> &colHash, heaps_info &hi,
                             paru_work *Work, ParU_Numeric *Num);

void paru_assemble_all(Int e, Int f, std::vector<Int> &colHash, paru_work *Work,
                       ParU_Numeric *Num);

void paru_assemble_cols(Int e, Int f, std::vector<Int> &colHash,
                        paru_work *Work, ParU_Numeric *Num);

void paru_assemble_rows(Int e, Int f, std::vector<Int> &colHash,
                        paru_work *Work, ParU_Numeric *Num);

void paru_assemble_el_with0rows(Int e, Int f, std::vector<Int> &colHash,
                                paru_work *Work, ParU_Numeric *Num);

void paru_full_summed(Int e, Int f, paru_work *Work, ParU_Numeric *Num);

// heap related
ParU_Ret paru_make_heap(Int f, Int start_fac,
                        std::vector<Int> &pivotal_elements, heaps_info &hi,
                        std::vector<Int> &colHash, paru_work *Work,
                        ParU_Numeric *Num);

ParU_Ret paru_make_heap_empty_el(Int f, std::vector<Int> &pivotal_elements,
                                 heaps_info &hi, paru_work *Work,
                                 ParU_Numeric *Num);

// hash related
void paru_insert_hash(Int key, Int value, std::vector<Int> &colHash);
Int paru_find_hash(Int key, std::vector<Int> &colHash, Int *fcolList);

// permutation stuff for the solver
ParU_Ret paru_perm(ParU_Symbolic *Sym, ParU_Numeric *Num);
Int paru_apply_inv_perm(const Int *P, const double *b, double *x, Int m);
Int paru_apply_inv_perm(const Int *P, const double *b, double *x, Int m, Int n);
Int paru_apply_perm_scale(const Int *P, const double *s, const double *b,
                          double *x, Int m);
Int paru_apply_perm_scale(const Int *P, const double *s, const double *b,
                          double *x, Int m, Int n);
Int paru_gaxpy(cholmod_sparse *A, const double *x, double *y, double alpha);
double paru_spm_1norm(cholmod_sparse *A);
double paru_vec_1norm(const double *x, Int n);

void paru_Diag_update(Int pivcol, Int pivrow, paru_work *Work);
void paru_tasked_dgemm(Int f, BLAS_INT m, BLAS_INT n, BLAS_INT k, double *A,
                       BLAS_INT lda, double *B, BLAS_INT ldb, double beta,
                       double *C, BLAS_INT ldc, paru_work *Work,
                       ParU_Numeric *Num);

void paru_tasked_trsm(Int f, int m, int n, double alpha, double *a, int lda,
                      double *b, int ldb, paru_work *Work, ParU_Numeric *Num);

ParU_Ret paru_free_work(ParU_Symbolic *Sym, paru_work *Work);
// lsolve and usolve
Int paru_lsolve(double *x, ParU_Symbolic *Sym, ParU_Numeric *Num,
                ParU_Control *Control);
Int paru_lsolve(double *X, Int n, ParU_Symbolic *Sym, ParU_Numeric *Num,
                ParU_Control *Control);
Int paru_usolve(double *x, ParU_Symbolic *Sym, ParU_Numeric *Num,
                ParU_Control *Control);
Int paru_usolve(double *X, Int n, ParU_Symbolic *Sym, ParU_Numeric *Num,
                ParU_Control *Control);

// not user-callable: for testing only
ParU_Ret paru_backward(double *x1, double &resid, double &norm,
                       cholmod_sparse *A, ParU_Symbolic *Sym, ParU_Numeric *Num,
                       ParU_Control *Control);

void paru_write(int scale, char *id, paru_work *Work, ParU_Numeric *Num);
#endif
