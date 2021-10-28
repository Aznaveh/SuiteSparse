////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_internal.hpp //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#ifndef PARU_INTERNAL_H
#define PARU_INTERNAL_H
//!
//  internal libraries that are not visible to the user
//  @author Aznaveh
//
#include "ParU.hpp"

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

#ifdef Int  // defined in amd
#undef Int
#endif

#define DLONG
extern "C"
{
#include "cholmod_blas.h"
#include "umf_internal.h"
}

// for printing information uncomment this; to activate assertions uncomment
//#undef NPR  //<<1>>

// from spqr.hpp
// Aznaveh For MATLAB OUTPUT UNCOMMENT HERE
// uncomment the following line to turn on debugging
//#undef NDEBUG  //<<2>>

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
#else
#define PRLEVEL(level, param)
#define DEBUGLEVEL(level)
#endif

// -----------------------------------------------------------------------------
// basic macros
// -----------------------------------------------------------------------------

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define EMPTY (-1)
// defined in amd #define TRUE 1
// defined in amd #define FALSE 0
#define IMPLIES(p, q) (!(p) || (q))

// NULL should already be defined, but ensure it is here.
#ifndef NULL
#define NULL ((void *)0)
#endif

#define Size_max ((size_t)(-1))  // the largest value of size_t
//------------------------------------------------------------------------------
// inline internal functions

inline Int *colIndex_pointer(paru_Element *curEl) { return (Int *)(curEl + 1); }
// Never ever use these functions prior to initializing ncols and nrows
inline Int *rowIndex_pointer(paru_Element *curEl)
{
    return (Int *)(curEl + 1) + curEl->ncols;
}

inline Int *relColInd(paru_Element *curEl)
//{    return (Int*)(curEl+1) + curEl->ncols + curEl->nrows + 1;}
{
    return (Int *)(curEl + 1) + curEl->ncols + curEl->nrows;
}

inline Int *relRowInd(paru_Element *curEl)
//{    return (Int*)(curEl+1) + 2*curEl->ncols + curEl->nrows + 2;}
{
    return (Int *)(curEl + 1) + 2 * curEl->ncols + curEl->nrows;
}

inline double *numeric_pointer(paru_Element *curEl)
// sizeof Int and double are same, but I keep it like this for clarity
//{ return (double*)((Int*)(curEl+1) + 2*curEl->ncols + 2*curEl->nrows + 2);}
{
    return (double *)((Int *)(curEl + 1) + 2 * curEl->ncols + 2 * curEl->nrows);
}

inline Int flip(Int colInd) { return -colInd - 2; }

inline Int lac_el(paru_Element **elementList, Int eli)
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
void paru_free_el(Int e, paru_Element **elementList);

void *operator new(std::size_t sz);
void operator delete(void *ptr) noexcept;

void paru_memset(void *ptr, Int value, size_t num);
void paru_memcpy(void *destination, const void *source, size_t num);

/* add tuple functions defintions */
Int paru_add_rowTuple(tupleList *RowList, Int row, paru_Tuple T);

Int paru_dgetrf(double *F, Int *frowList, Int m, Int n, BLAS_INT *ipiv);

Int paru_factorize_full_summed(Int f, Int start_fac,
                               std::vector<Int> &panel_row,
                               std::set<Int> &stl_colSet,
                               std::vector<Int> &pivotal_elements,
                               paru_matrix *paruMatInfo);

paru_Element *paru_create_element(Int nrows, Int ncols, Int init);

void paru_assemble_row_2U(Int e, Int f, Int sR, Int dR,
                          std::vector<Int> &colHash, paru_matrix *paruMatInfo);

Int paru_trsm(Int f, double *pF, double *uPart, Int fp, Int rowCount,
              Int colCount);
Int paru_dgemm(Int f, double *pF, double *uPart, double *el, Int fp,
               Int rowCount, Int colCount);

void paru_print_element(paru_matrix *paruMatInfo, Int e);
void paru_print_tupleList(tupleList *listSet, Int index);
void paru_init_rel(Int f, paru_matrix *paruMatInfo);

void paru_update_rel_ind_col(Int e, Int f, std::vector<Int> &colHash,
                             paru_matrix *paruMatInfo);

void paru_update_rowDeg(Int panel_num, Int row_end, Int f, Int start_fac,
                        std::set<Int> &stl_colSet,
                        std::vector<Int> &pivotal_elements,
                        paru_matrix *paruMatInfo);

Int paru_cumsum(Int n, Int *X);

Int bin_srch_col(Int *srt_lst, Int l, Int r, Int num);
Int bin_srch(Int *srt_lst, Int l, Int r, Int num);

ParU_ResultCode paru_init_rowFronts(paru_matrix **paruMatInfo_handle,
                                    cholmod_sparse *A, paru_symbolic *LUsym);
ParU_ResultCode paru_front(Int f, paru_matrix *paruMatInfo);

ParU_ResultCode paru_pivotal(std::vector<Int> &pivotal_elements,
                             std::vector<Int> &panel_row, Int &zero_piv_rows,
                             Int f, heaps_info &hi, paru_matrix *paruMatInfo);

int paru_intersection(Int e, paru_Element **elementList,
                      std::set<Int> &stl_colSet);

ParU_ResultCode paru_prior_assemble(Int f, Int start_fac,
                                    std::vector<Int> &pivotal_elements,
                                    std::vector<Int> &colHash, heaps_info &hi,
                                    paru_matrix *paruMatInfo);

void paru_assemble_all(Int e, Int f, std::vector<Int> &colHash,
                       paru_matrix *paruMatInfo);

void paru_assemble_cols(Int e, Int f, std::vector<Int> &colHash,
                        paru_matrix *paruMatInfo);

void paru_assemble_rows(Int e, Int f, std::vector<Int> &colHash,
                        paru_matrix *paruMatInfo);

void paru_assemble_el_with0rows(Int e, Int f, std::vector<Int> &colHash,
                                paru_matrix *paruMatInfo);

void paru_full_summed(Int e, Int f, paru_matrix *paruMatInfo);

// heap related
ParU_ResultCode paru_make_heap(Int f, Int start_fac,
                               std::vector<Int> &pivotal_elements,
                               heaps_info &hi, std::vector<Int> &colHash,
                               paru_matrix *paruMatInfo);

ParU_ResultCode paru_make_heap_empty_el(Int f,
                                        std::vector<Int> &pivotal_elements,
                                        heaps_info &hi,
                                        paru_matrix *paruMatInfo);

// hash related
void paru_insert_hash(Int key, Int value, std::vector<Int> &colHash);
Int paru_find_hash(Int key, std::vector<Int> &colHash, Int *fcolList);

// permutation stuff for the solver
ParU_ResultCode paru_perm(paru_matrix *paruMatInfo);
Int paru_apply_inv_perm(const Int *P, const double *b, double *x, Int m);
Int paru_apply_perm_scale(const Int *P, const double *s, const double *b,
                          double *x, Int m);

// lsolve and usolve
Int paru_lsolve(paru_matrix *paruMatInfo, double *x);
Int paru_usolve(paru_matrix *paruMatInfo, double *x);

Int paru_gaxpy(cholmod_sparse *A, const double *x, double *y, double alpha);
double paru_spm_1norm(cholmod_sparse *A);
double paru_vec_1norm(const double *x, Int n);
void paru_Diag_update(Int pivcol, Int pivrow, paru_matrix *paruMatInfo);
void paru_tasked_dgemm(Int f, char *transa, char *transb, BLAS_INT *m,
                       BLAS_INT *n, BLAS_INT *k, double *alpha, double *A,
                       BLAS_INT *lda, double *B, BLAS_INT *ldb, double *beta,
                       double *C, BLAS_INT *ldc);
void paru_tasked_trsm(Int f, char *side, char *uplo, char *transa, char *diag,
                      int *m, int *n, double *alpha, double *a, int *lda,
                      double *b, int *ldb);
#endif
