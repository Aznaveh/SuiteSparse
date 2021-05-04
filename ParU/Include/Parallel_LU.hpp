// ============================================================================/
// ======================= Parallel_LU.hpp ====================================/
// ============================================================================/

#ifndef PARU_H
#define PARU_H

// These libraries are included probably in Suitesparse config
//#include <stdlib.h>
//#include <math.h>
//#include <float.h>
//#include <stdio.h>
//#include <cstring>

// To be able to use set
#include <algorithm>
#include <set>
#include <vector>

#include <malloc.h>

#include <omp.h>

extern "C"
{
// #include "umfpack.h"
#include "cholmod.h"
#include "cholmod_blas.h"
#include "umf_internal.h"
}

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

// for printing information uncomment this; to activate assertions uncomment
//#undef NPR    //<<1>>

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

#ifdef Int  // defined in amd
#undef Int
#endif
#define Int SuiteSparse_long

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

// =============================================================================
// === paru_symbolic ===========================================================
// =============================================================================

// The contents of this object do not change during numeric factorization.  The
// Symbolic object depends only on the pattern of the input matrix, and not its
// values.
// This makes parallelism easier to manage, since all threads can
// have access to this object without synchronization.
//
typedef struct
{ /* paru_symbolic*/

    // -------------------------------------------------------------------------
    // row-form of the input matrix and its permutations
    // -------------------------------------------------------------------------

    // During symbolic analysis, the nonzero pattern of S = A(P,Q) is
    // constructed, where A is the user's input matrix.  Its numerical values
    // are also constructed, but they do not become part of the Symbolic
    // object.  The matrix S is stored in row-oriented form.  The rows of S are
    // sorted according to their leftmost column index (via Pinv).  Column
    // indices in each row of S are in strictly ascending order, even though
    // the input matrix A need not be sorted.

    Int m, n, anz;  // S is m-by-n with anz entries

    Int snz;     // nnz in submatrix
    Int *Sp;     // size m+1-n1, row pointers of S
    Int *Sj;     // size snz = Sp [n], column indices of S
    double *Sx;  // size snz = Sp [n], numeric values of S

    Int *Qfill;  // size n, fill-reducing column permutation.
    // Qfill [k] = j if column k of A is column j of S.

    Int *Pinv;  // size m, inverse row permutation that places
    // S=A(P,Q) in increasing order of leftmost column
    // index.  Pinv [i] = k if row i of A is row k of S.

    Int *Pinit;  // size m, row permutation.
    // UMFPACK computes it and I compute Pinv out of it.
    // I need it in several places so I decided to keep it

    Int *Sleft;  // size n-n1+2.  The list of rows of S whose
    // leftmost column index is j is given by
    // Sleft [j] ... Sleft [j+1]-1.  This can be empty (that is, Sleft
    // [j] can equal Sleft [j+1]).  Sleft [n] is the number of
    // non-empty rows of S, and Sleft [n+1] == m.  That is, Sleft [n]
    // ... Sleft [n+1]-1 gives the empty rows of S, if any.

    Int strategy;  // for this package it is important if the strategy is
    // symmetric or if it is unsymmetric

    // -------------------------------------------------------------------------
    // frontal matrices: pattern and tree
    // -------------------------------------------------------------------------

    // Each frontal matrix is fm-by-fn, with fnpiv pivot columns.  The fn
    // column indices are given by a set of size fnpiv pivot columns, defined
    // by Super, followed by the pattern Rj [ Rp[f] ...  Rp[f+1]-1 ].

    // The row indices of the front are not kept.  If the Householder vectors
    // are not kept, the row indices are not needed.  If the Householder
    // vectors are kept, the row indices are computed dynamically during
    // numerical factorization.

    Int nf;        // number of frontal matrices; nf <= MIN (m,n)
    Int n1;        // number of singletons in the matrix
                   // the matrix S is the one without any singletons
    Int rs1, cs1;  // number of row and column singletons, n1 = rs1+cs1;

    // parent, child, and childp define the row merge tree or etree (A'A)
    Int *Parent;  // size nf+1  Add another node just to make the forest a
    Int *Child;   // size nf+1      tree
    Int *Childp;  // size nf+2

    // The parent of a front f is Parent [f], or EMPTY if f=nf.
    // A list of children of f can be obtained in the list
    // Child [Childp [f] ... Childp [f+1]-1].

    // Node nf in the tree is a placeholder; it does not represent a frontal
    // matrix.  All roots of the frontal "tree" (may be a forest) have the
    // placeholder node nf as their parent.  Thus, the tree of nodes 0:nf is
    // truly a tree, with just one parent (node nf).

    Int *aParent;  // size m+nf
    Int *aChild;   // size m+nf+1
    Int *aChildp;  // size m+nf+2
    Int *first;  // size m+nf first successor of front in augmented postordered
    //  tree; all successors are between first[eli]...eli-1

    // pivot column in the front F.  This refers to a column of S.  The
    // number of expected pivot columns in F is thus
    // Super [f+1] - Super [f].

    // Upper bound number of rows for each front
    Int *Fm;  // size nf+1

    // Upper bound  number of rows in the contribution block of each front
    Int *Cm;  // size nf+1

    Int *Super;  // size nf+1.  Super [f] gives the first
    // pivot column in the front F.  This refers to a column of S.  The
    // number of expected pivot columns in F is thus
    // Super [f+1] - Super [f].

    Int *Rp;  // size nf+1
    Int *Rj;  // size rjsize; compressed supernodal form of R

    Int rjsize;  // size of Rj

    Int *row2atree;    // Mapping from rows to augmented tree size m
    Int *super2atree;  // Mapping from super nodes to augmented tree size nf

    Int *Chain_start;  // size = n_col +1;  actual size = nfr+1
    // The kth frontal matrix chain consists of frontal
    // matrices Chain_start [k] through Chain_start [k+1]-1.
    // Thus, Chain_start [0] is always 0 and
    // Chain_start[nchains] is the total number of frontal
    // matrices, nfr. For two adjacent fornts f and f+1
    // within a single chian, f+1 is always the parent of f
    // (that is, Front_parent [f] = f+1).

    Int *Chain_maxrows;  // size = n_col +1;  actual size = nfr+1
    Int *Chain_maxcols;  // The kth frontal matrix chain requires a single
    // working array of dimension Chain_maxrows [k] by
    // Chain_maxcols [k], for the unifrontal technique that
    // factorizes the frontal matrix chain. Since the
    // symbolic factorization only provides

    Int Us_bound_size;   // Upper bound on size of all Us, sum all fp*fn
    Int LUs_bound_size;  // Upper bound on size of all LUs, sum all fp*fm
    Int row_Int_bound;   // Upper bound on size of all ints for rows
    Int col_Int_bound;   // Upper bound on size of all ints for cols
} paru_symbolic;

// =============================================================================
//      paru_Tuple, Row and Column data structure
// =============================================================================
typedef struct
{ /* paru_Tuple */

    /* The (e,f) tuples for element lists */
    Int e, /*  element number */
        f; /*   offest */
} paru_Tuple;

/* -------------------------------------------------------------------------- */
/* An element */
/* -------------------------------------------------------------------------- */

typedef struct
{ /* paru_Element */

    Int

        nrowsleft, /* number of rows remaining */
        ncolsleft, /* number of columns remaining */
        nrows,     /* number of rows */
        ncols,     /* number of columns */
        rValid,    /* validity of relative row index */
        cValid;    /* validity of relative column index */

    Int lac;  // least active column which is active
              // 0 <= lac <= ncols

    /* followed in memory by:
       Int
       col [0..ncols-1],	column indices of this element
       row [0..nrows-1] ;	row indices of this element

       relColInd [0..ncols-1];	relative indices of this element for
       current front
       relRowInd [0..nrows-1],	relative indices of this element for
       current front
       double ncols*nrows; numeric values
       */

} paru_Element;  // contribution block

typedef struct
{
    Int sum_size, biggest_Child_size, biggest_Child_id;
} heaps_info;
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

typedef struct
{ /*List of tuples */

    /*element of a column or a row*/
    Int numTuple,     /*  number of Tuples in this element */
        len;          /*  length of allocated space for current list*/
    paru_Tuple *list; /* list of tuples regarding to this element */

} tupleList;

typedef struct
{ /*work_struct*/

    // gather scatter space for rows
    Int *rowSize;  // Initalized data structure, size of rows
    // Int rowMark;      // Work->rowSize[x] < rowMark[eli] for each front
    Int *rowMark;  // size = m+nf

    Int *elRow;  // Initalized data structure, size m+nf
    Int *elCol;  // Initalized data structure, size m+nf

} work_struct;

typedef struct
{
    size_t size_bank[64],  // size of each bank
        remaining;         // remaining space of current bank
    Int cur;               // points to current bank
    void *mem_bank[64];
    void *avail;  // pointer to available space in current bank
} stack_mem_info;

typedef struct
{              /*dense factorized part pointer*/
    Int m, n;  /* mxn dense matrix */
    double *p; /* point to factorized parts */
} paru_fac;

typedef struct
{             /*Matrix */
    Int m, n; /* size of the sumbatrix that is factorized */
    paru_symbolic *LUsym;
    tupleList *RowList; /* size n of dynamic list */

    paru_Element **elementList; /* pointers to all elements, size = m+nf+1 */
    work_struct *Work;

    Int *time_stamp; /* for relative index update
                        not initialized */

    // Computed parts of each front
    Int *frowCount;        /* size nf   size(CB) = rowCount[f]x         */
    Int *fcolCount;        /* size nf                        colCount[f]*/
    Int **frowList;        /* size nf   frowList[f] is rows of the matrix S */
    Int **fcolList;        /* size nf   colList[f] is non pivotal cols of the
                              matrix S */
    paru_fac *partial_Us;  /* size nf   size(Us)= fp*colCount[f]    */
    paru_fac *partial_LUs; /* size nf   size(LUs)= rowCount[f]*fp   */

#ifndef NDEBUG
    Int actual_alloc_LUs;     /* actual memory allocated for LUs*/
    Int actual_alloc_Us;      /* actual memory allocated for Us*/
    Int actual_alloc_row_int; /* actual memory allocated for rows*/
    Int actual_alloc_col_int; /* actual memory allocated for cols*/
#endif

    Int *row_degree_bound; /* row degree size number of rows */
    Int panel_width;       /* width of panel for dense factorizaiton*/

    Int *lacList; /* sieze m+nf least active column of each element
                     el_colIndex[el->lac]  == lacList [e]
                     number of element*/

    double *scale_row; /* the array for row scaling */

    // each active front owns and manage a heap list. The heap is based on the
    // least numbered column. The active front Takes the pointer of the biggest
    // child and release its other children after concatenating their list to
    // its own. The list of heaps are initialized by nullptr
    std::vector<Int> **heapList; /* size m+nf+1, initialized with nullptr  */

    // analysis information
    double my_time;
    double umf_time;

    stack_mem_info stack_mem;

#ifdef COUNT_FLOPS
    // flop count info
    double flp_cnt_dgemm;
    double flp_cnt_trsm;
    double flp_cnt_dger;
    double flp_cnt_real_dgemm;
#endif

} paru_matrix;

paru_symbolic *paru_analyze(cholmod_sparse *A, cholmod_common *cc);

paru_matrix *paru_init_rowFronts(cholmod_sparse *A, int scale,
                                 paru_symbolic *LUsym, cholmod_common *cc);

/* Wrappers for managing memory */
void *paru_alloc(size_t n, size_t size, cholmod_common *cc);
void *paru_calloc(size_t n, size_t size, cholmod_common *cc);
void *paru_stack_calloc(size_t n, size_t size, paru_matrix *paruMatInfo,
                        cholmod_common *cc);
void *paru_realloc(Int newsize, Int size_Entry, void *oldP, Int *size,
                   cholmod_common *cc);

void paru_free(Int n, Int size, void *p, cholmod_common *cc);
void paru_free_el(Int e, paru_Element **elementList, cholmod_common *cc);
void paru_freesym(paru_symbolic **LUsym_handle, cholmod_common *cc);
void paru_freemat(paru_matrix **paruMatInfo_handle, cholmod_common *cc);

/* add tuple functions defintions */
Int paru_add_rowTuple(tupleList *RowList, Int row, paru_Tuple T,
                      cholmod_common *cc);

int paru_front(paru_matrix *paruMatInfo, Int f, cholmod_common *cc);

Int paru_dgetrf(double *F, Int *frowList, Int m, Int n, BLAS_INT *ipiv);
Int paru_factorize(Int f, Int start_fac, std::vector<Int> &panel_row,
                   std::set<Int> &stl_colSet,
                   std::vector<Int> &pivotal_elements,
                   paru_matrix *paruMatInfo);

paru_Element *paru_create_element(Int nrows, Int ncols, Int init,
                                  cholmod_common *cc);
void assemble_col(const double *sR, double *dR, Int m, const Int *relRowInd);

void assemble_row_toU(Int e, Int f, Int sR, Int dR, std::vector<Int> &colHash,
                      paru_matrix *paruMatInfo);

Int paru_trsm(double *pF, double *uPart, Int fp, Int rowCount, Int colCount);
Int paru_dgemm(double *pF, double *uPart, double *el, Int fp, Int rowCount,
               Int colCount);

void paru_print_element(paru_matrix *paruMatInfo, Int e);
void paru_print_tupleList(tupleList *listSet, Int index);
void paru_init_rel(Int f, paru_matrix *paruMatInfo);

void paru_update_rel_ind_col(Int e, Int f, std::vector<Int> &colHash,
                             paru_matrix *paruMatInfo);

void paru_write(paru_matrix *paruMatInfo, int scale, char *id,
                cholmod_common *cc);

void paru_update_rowDeg(Int panel_num, Int row_end, Int f, Int start_fac,
                        std::set<Int> &stl_colSet,
                        std::vector<Int> &pivotal_elements,
                        paru_matrix *paruMatInfo);

Int paru_cumsum(Int n, Int *X);

Int bin_srch_col(Int *srt_lst, Int l, Int r, Int num);
Int bin_srch(Int *srt_lst, Int l, Int r, Int num);

void paru_pivotal(std::vector<Int> &pivotal_elements,
                  std::vector<Int> &panel_row, Int f, heaps_info &hi,
                  paru_matrix *paruMatInfo, cholmod_common *cc);

int paru_intersection(Int e, paru_Element **elementList,
                      std::set<Int> &stl_colSet);

void paru_prior_assemble(Int f, Int start_fac,
                         std::vector<Int> &pivotal_elements,
                         std::vector<Int> &colHash, heaps_info &hi,
                         paru_matrix *paruMatInfo, cholmod_common *cc);

void paru_eliminate_all(Int e, Int f, std::vector<Int> &colHash,
                        paru_matrix *paruMatInfo, cholmod_common *cc);

void paru_eliminate_cols(Int e, Int f, std::vector<Int> &colHash,
                         paru_matrix *paruMatInfo, cholmod_common *cc);

void paru_eliminate_rows(Int e, Int f, std::vector<Int> &colHash,
                         paru_matrix *paruMatInfo, cholmod_common *cc);

// heap related
void paru_make_heap(Int f, Int start_fac, std::vector<Int> &pivotal_elements,
                    heaps_info &hi, std::vector<Int> &colHash,
                    paru_matrix *paruMatInfo, cholmod_common *cc);

// hash related
void paru_insert_hash(Int key, Int value, std::vector<Int> &colHash);
Int paru_find_hash(Int key, std::vector<Int> &colHash, Int *fcolList);
#endif
