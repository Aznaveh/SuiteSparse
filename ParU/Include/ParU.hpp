// ============================================================================/
// ======================= ParU.hpp ===========================================/
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

//#include <malloc.h> // mallopt used in paru_init_rowFronts.cpp

#include <omp.h>

// Not using definitions in /CHOLMOD/Include/cholmod_blas.h
// using cbals instead
#define CHOLMOD_BLAS_H 
#ifdef BLAS_INT 
#undef BLAS_INT
#endif
#define BLAS_INT int32_t

#ifdef MKLROOT
#define CHOLMOD_BLAS_H
#include <mkl.h>
#define BLAS_set_num_threads(n) mkl_set_num_threads(n)
#else
#include <cblas.h>
#define BLAS_set_num_threads(n) openblas_set_num_threads(n)
#endif

extern "C"
{
#include "cholmod.h"
#include "umfpack.h"
}

#ifdef Int  // defined in amd
#undef Int
#endif
#define Int int64_t

//  Just like UMFPACK_STRATEGY defined in UMFPACK/Include/umfpack.h
//  However I might use different strategies in symbolic analysis and
//  factorization; they match with UMFPACK_STRATEGY
#define PARU_STRATEGY_AUTO 0         // decided to use sym. or unsym. strategy
#define PARU_STRATEGY_UNSYMMETRIC 1  // COLAMD(A), metis, hmetis(?)
#define PARU_STRATEGY_SYMMETRIC 3    // prefer diagonal

// =============================================================================
// === paru_symbolic ===========================================================
// =============================================================================

// The contents of this object do not change during numeric factorization.  The
// Symbolic object depends only on the pattern of the input matrix, and not its
// values.
// This makes parallelism easier to manage, since all threads can
// have access to this object without synchronization.
//

// typedef struct U_singleton
struct U_singleton
{
    // CSR format for U singletons
    Int nnz;      // nnz in submatrix
    Int *Sup;     // size cs1
    Int *Suj;     // size ?
    double *Sux;  // size ?
};

// typedef struct L_singleton
struct L_singleton
{
    // CSC format for U singletons
    Int nnz;      // nnz in submatrix
    Int *Slp;     // size rs1
    Int *Sli;     // size ?
    double *Slx;  // size ?
};

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

    Int m, n, anz;  // S is m-by-n with anz entries; S is scaled

    Int snz;     // nnz in submatrix
    Int *Sp;     // size m+1-n1, row pointers of S
    Int *Sj;     // size snz = Sp [n], column indices of S
    double *Sx;  // size snz = Sp [n], numeric values of S

    double *scale_row;  // the array for row scaling based on original matrix
                        // size = m

    // Usingletons and Lsingltons
    U_singleton ustons;
    L_singleton lstons;

    Int *Qfill;  // size n, fill-reducing column permutation.
    // Qfill [k] = j if column k of A is column j of S.

    Int *Pinit;  // size m, row permutation.
    // UMFPACK computes it and I compute Pinv out of it.
    // I need it in several places so I decided to keep it

    Int *Diag_map;  // size n,
    // UMFPACK computes it and I use it to find original diags out of it

    Int *Ps;  // size m, row permutation.
    // Permutation from S to LU. needed for lsolve and usolve
    // Look paru_perm for more details

    Int *Pfin;  // size m, row permutation.
    // ParU final permutation. Look paru_perm for more details

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
    // column indices are given by a set of size fnpiv pivot columns

    // The row indices of the front are not kept.  If the Householder vectors
    // are not kept, the row indices are not needed.  If the Householder
    // vectors are kept, the row indices are computed dynamically during
    // numerical factorization.

    Int nf;  // number of frontal matrices; nf <= MIN (m,n)
    Int n1;  // number of singletons in the matrix
    // the matrix S is the one without any singletons
    Int rs1, cs1;  // number of row and column singletons, n1 = rs1+cs1;

    // parent, child, and childp define the row merge tree or etree (A'A)
    Int *Parent;  // size nf+1  Add another node just to make the forest a
    Int *Child;   // size nf+1      tree
    Int *Childp;  // size nf+2

    Int *Depth;  // size nf distance of each node from the root

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
    Int *first;    // size nf+1 first successor of front in the tree;
    // all successors are between first[f]...f-1

    // pivot column in the front F.  This refers to a column of S.  The
    // number of expected pivot columns in F is thus
    // Super [f+1] - Super [f].

    //Int num_roots;  // number of roots
    // it is at least one and can be more in case of forest
    Int *roots;

    // Upper bound number of rows for each front
    Int *Fm;  // size nf+1

    // Upper bound  number of rows in the contribution block of each front
    Int *Cm;  // size nf+1

    Int *Super;  // size nf+1.  Super [f] gives the first
    // pivot column in the front F.  This refers to a column of S.  The
    // number of expected pivot columns in F is thus
    // Super [f+1] - Super [f].

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

    // only used for statistics when debugging is enabled:
    Int Us_bound_size;   // Upper bound on size of all Us, sum all fp*fn
    Int LUs_bound_size;  // Upper bound on size of all LUs, sum all fp*fm
    Int row_Int_bound;   // Upper bound on size of all ints for rows
    Int col_Int_bound;   // Upper bound on size of all ints for cols

    double *front_flop_bound;  // bound on m*n*k for each front size nf+1
    double *stree_flop_bound;  // flop bound for front and descendents size nf+1

    //data structure related to tasks
    Int ntasks; // number of tasks; at most nf
    Int *task_map; // each task does the fronts 
                   // from task_map[i]+1 to task_map[i+1]; task_map[0] is -1
    Int *task_parent; //tree data structure for tasks
    Int *task_num_child; //number of children of each task
    Int *task_depth; //max depth of each task

    // symbolic analysis time
    double my_time;

    Int max_chain; //maximum size of the chains in final tree

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

// internal:

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

    // gather scatter space for elements
    Int *elRow;  // Initalized data structure, size m+nf
    Int *elCol;  // Initalized data structure, size m+nf

} work_struct;

typedef struct
{              /*dense factorized part pointer*/
    Int m, n;  /* mxn dense matrix */
    double *p; /* point to factorized parts */
} paru_fac;

enum ParU_ResultCode
{
    PARU_SUCCESS,
    PARU_OUT_OF_MEMORY,
    PARU_INVALID,
    PARU_SINGULAR
};

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

    // only used for statistics when debugging is enabled:
    Int actual_alloc_LUs;     /* actual memory allocated for LUs*/
    Int actual_alloc_Us;      /* actual memory allocated for Us*/
    Int actual_alloc_row_int; /* actual memory allocated for rows*/
    Int actual_alloc_col_int; /* actual memory allocated for cols*/

    Int *row_degree_bound; /* row degree size number of rows */
    Int panel_width;       /* width of panel for dense factorizaiton*/

    Int *lacList; /* sieze m+nf least active column of each element
                     el_colIndex[el->lac]  == lacList [e]
                     number of element*/

    Int *Diag_map;  // size n,
    // Both of these are NULL if the stratey is not symmetric
    // copy of Diag_map from LUsym;
    // this copy can be updated during the factorization
    Int *inv_Diag_map;  // size n,
    // inverse of Diag_map from LUsym;
    // It helps editing the Diag_map

    // each active front owns and manage a heap list. The heap is based on the
    // least numbered column. The active front Takes the pointer of the biggest
    // child and release its other children after concatenating their list to
    // its own. The list of heaps are initialized by nullptr
    std::vector<Int> **heapList; /* size m+nf+1, initialized with nullptr  */

    // analysis information
    double my_time;  // factorization time
    double umf_time;

    // #ifdef COUNT_FLOPS
    // flop count info
    double flp_cnt_dgemm;
    double flp_cnt_trsm;
    double flp_cnt_dger;
    double flp_cnt_real_dgemm;
    // #endif

    Int naft; //number of actvie frontal tasks
    Int resq; //number of remainig ready tasks in the queue
    Int paru_max_threads; //I want to call omp_get_max_threads just once
                         // or user can give me a value less than that
    ParU_ResultCode res;

} paru_matrix;

//------------------------------------------------------------------------------
// user:

paru_symbolic *paru_analyze(Int scale, cholmod_sparse *A);

/* usage:
   S = paru_analyse (A) ;
   LU = paru_factoriz (A,S) ;
info: an enum: PARU_SUCCESS, PARU_OUT_OF_MEMORY, PARU_INVALID, PARU_SINGULAR,
... info = paru_analyse (&S, A) ; info = paru_factoriz (&LU, A,S) ;
*/

// a routine that does init_row and also factorization
ParU_ResultCode paru_factorize(cholmod_sparse *A, paru_symbolic *LUsym,
                               paru_matrix **paruMatInfo_handle);
ParU_ResultCode paru_solve(paru_matrix *paruMatInfo, double *b);
ParU_ResultCode paru_solve(paru_matrix *paruMatInfo, double *b, Int n);
void paru_write(paru_matrix *paruMatInfo, int scale, char *id);

void paru_freesym(paru_symbolic **LUsym_handle);
void paru_freemat(paru_matrix **paruMatInfo_handle);

ParU_ResultCode paru_residual(cholmod_sparse *A, paru_matrix *paruMatInfo,
                              double *b, double *Results);
ParU_ResultCode paru_residual(cholmod_sparse *A, paru_matrix *paruMatInfo,
                              double *b, double *Results, Int n);


ParU_ResultCode paru_backward(cholmod_sparse *A, paru_matrix *paruMatInfo,
                              double *x1, double *Results);
#endif
