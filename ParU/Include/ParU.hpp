// ============================================================================/
// ======================= ParU.hpp ===========================================/
// ============================================================================/

#ifndef PARU_H
#define PARU_H

// These libraries are included probably in Suitesparse_config
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

#define CHOLMOD_BLAS_H
#ifdef BLAS_INT
#undef BLAS_INT
#endif

#ifdef MKLROOT
// MKL BLAS
#define CHOLMOD_BLAS_H
#include <mkl.h>
#define BLAS_set_num_threads(n) mkl_set_num_threads(n)
#define BLAS_INT MKL_INT
#else
// assume OpenBLAS
#include <cblas.h>
#define BLAS_set_num_threads(n) openblas_set_num_threads(n)
#define BLAS_INT int
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
#define PARU_STRATEGY_AUTO 0         // decided to use sym. or unsym. strategy
#define PARU_STRATEGY_UNSYMMETRIC 1  // COLAMD(A), metis, ...
#define PARU_STRATEGY_SYMMETRIC 3    // prefer diagonal

// =============================================================================
// === ParU_Symbolic ===========================================================
// =============================================================================
//
// The contents of this object do not change during numeric factorization.  The
// ParU_U_singleton and ParU_L_singleton are datastructures for singletons that
// has been borrowed from UMFPACK
//
//              ParU_L_singleton is CSC
//                                    l
//                                    v
// 	  ParU_U_singleton is CSR -> L U U U U U U U U
// 	                             . L U U U U U U U
// 	                             . . L U U U U U U
// 	                             . . . L . . . . .
// 	                             . . . L L . . . .
// 	                             . . . L L x x x x
// 	                             . . . L L x x x x
// 	                             . . . L L x x x x
// 	                             . . . L L x x x x

struct ParU_U_singleton
{
    // CSR format for U singletons
    Int nnz;   // nnz in submatrix
    Int *Sup;  // size cs1
    Int *Suj;  // size is computed
    double *Sux;
};

struct ParU_L_singleton
{
    // CSC format for U singletons
    Int nnz;   // nnz in submatrix
    Int *Slp;  // size rs1
    Int *Sli;  // size is computed
    double *Slx;
};

struct ParU_Symbolic
{
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
    ParU_U_singleton ustons;
    ParU_L_singleton lstons;

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

    Int strategy;  // the strategy USED by umfpack
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

    // Int num_roots;  // number of roots
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

    // data structure related to tasks
    Int ntasks;        // number of tasks; at most nf
    Int *task_map;     // each task does the fronts
                       // from task_map[i]+1 to task_map[i+1]; task_map[0] is -1
    Int *task_parent;  // tree data structure for tasks
    Int *task_num_child;  // number of children of each task
    Int *task_depth;      // max depth of each task

    // symbolic analysis time
    double my_time;

    Int max_chain;  // maximum size of the chains in final tree
};

// =============================================================================
//      ParU_Tuple, Row and Column data structure
// =============================================================================
struct ParU_Tuple
{
    // The (e,f) tuples for element lists
    Int e,  //  element number
        f;  //   offest
};

// =============================================================================
// An element, contribution block
// =============================================================================
struct ParU_Element
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
    //   col [0..ncols-1],	column indices of this element
    //   row [0..nrows-1] ;	row indices of this element
    //
    //   relColInd [0..ncols-1];	relative indices of this element for
    //   current front
    //   relRowInd [0..nrows-1],	relative indices of this element for
    //   current front
    //   double ncols*nrows; numeric values
};

struct ParU_TupleList
{                      // List of tuples
    Int numTuple,      //  number of Tuples in this element
        len;           //  length of allocated space for current list
    ParU_Tuple *list;  // list of tuples regarding to this element
};

struct Paru_Work
{
    // gather scatter space for rows
    Int *rowSize;  // Initalized data structure, size of rows
    // Int rowMark;      // Work->rowSize[x] < rowMark[eli] for each front
    Int *rowMark;  // size = m+nf

    // gather scatter space for elements
    Int *elRow;  // Initalized data structure, size m+nf
    Int *elCol;  // Initalized data structure, size m+nf
};

struct ParU_Factors
{               // dense factorized part pointer
    Int m, n;   //  mxn dense matrix
    double *p;  //  point to factorized parts
};

enum ParU_Ret
{
    PARU_SUCCESS,
    PARU_OUT_OF_MEMORY,
    PARU_INVALID,
    PARU_SINGULAR
};

struct ParU_Control
{
    Int mem_chunk = 1024*1024; //chunk size for memset and memcpy
    //Sybmolic controls
    Int scale = 1; // if 1 matrix will be scaled using max_row
    Int umfpack_ordering = UMFPACK_ORDERING_METIS;
    Int umfpack_strategy = UMFPACK_STRATEGY_AUTO; //symmetric or unsymmetric


    //Numeric controls
    Int panel_width = 32;  // width of panel for dense factorizaiton
    Int paru_strategy = PARU_STRATEGY_AUTO; //the same stratey umfpack used

    double piv_toler = 0.1;    //tolerance for accepting sparse pivots
    double diag_toler = 0.001; //tolerance for accepting symmetric pivots
    Int trivial = 4; // dgemms with sizes less than trivial doesn't call BLAS
    Int worthwhile_dgemm = 512; // dgemms bigger than worthwhile are tasked
    Int worthwhile_trsm = 4096; // trsm bigger than worthwhile are tasked
    Int paru_max_threads = 0;  //It will be initialized with omp_max_threads
                               // if the user did not provide
};

struct ParU_Numeric
{
    Int m, n;  // size of the sumbatrix that is factorized
    ParU_Symbolic *Sym;
    ParU_TupleList *RowList;  // size n of dynamic list
    ParU_Control *Control;   // a copy of controls for internal use
                            // it is freed after factorize

    ParU_Element **elementList;  // pointers to all elements, size = m+nf+1
    Paru_Work *Work;

    Int *time_stamp;  // for relative index update; not initialized

    // Computed parts of each front
    Int *frowCount;  // size nf   size(CB) = rowCount[f]x
    Int *fcolCount;  // size nf                        colCount[f]
    Int **frowList;  // size nf   frowList[f] is rows of the matrix S
    Int **fcolList;  // size nf   colList[f] is non pivotal cols of the
                     //   matrix S
    ParU_Factors *partial_Us;   // size nf   size(Us)= fp*colCount[f]
    ParU_Factors *partial_LUs;  // size nf   size(LUs)= rowCount[f]*fp

    // only used for statistics when debugging is enabled:
    Int actual_alloc_LUs;      // actual memory allocated for LUs
    Int actual_alloc_Us;       // actual memory allocated for Us
    Int actual_alloc_row_int;  // actual memory allocated for rows
    Int actual_alloc_col_int;  // actual memory allocated for cols

    Int max_row_count;      // maximum number of rows/cols for all the fronts
    Int max_col_count;      // it is initalized after factorization
    Int *row_degree_bound;  // row degree size number of rows

    Int *lacList;  // sieze m+nf least active column of each element
                   //    el_colIndex[el->lac]  == lacList [e]
                   //    number of element

    Int *Diag_map;  // size n,
    // Both of these are NULL if the stratey is not symmetric
    // copy of Diag_map from Sym;
    // this copy can be updated during the factorization
    Int *inv_Diag_map;  // size n,
    // inverse of Diag_map from Sym;
    // It helps editing the Diag_map

    // each active front owns and manage a heap list. The heap is based on the
    // least numbered column. The active front Takes the pointer of the biggest
    // child and release its other children after concatenating their list to
    // its own. The list of heaps are initialized by nullptr
    std::vector<Int> **heapList;  // size m+nf+1, initialized with nullptr

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

    Int naft;              // number of actvie frontal tasks
    Int resq;              // number of remainig ready tasks in the queue
    ParU_Ret res;          // returning value of numeric phase
};

//------------------------------------------------------------------------------
ParU_Ret ParU_Analyze(cholmod_sparse *A, ParU_Symbolic **Sym_handle,
                      ParU_Control *Control);
ParU_Ret ParU_Factorize(cholmod_sparse *A, ParU_Symbolic *Sym,
                        ParU_Numeric **Num_handle, ParU_Control *Control);
ParU_Ret ParU_Solve(double *b, ParU_Numeric *Num, ParU_Control *Control);
ParU_Ret ParU_Solve(double *B, Int n, ParU_Numeric *Num, ParU_Control *Control);

ParU_Ret ParU_Freesym(ParU_Symbolic **Sym_handle, ParU_Control *Control);
ParU_Ret ParU_Freenum(ParU_Numeric **Num_handle, ParU_Control *Control);

ParU_Ret ParU_Residual(double *b, double &resid, double &norm,
                       cholmod_sparse *A, ParU_Numeric *Num,
                       ParU_Control *Control);

ParU_Ret ParU_Residual(cholmod_sparse *A, ParU_Numeric *Num, double *b,
                       double *Results, Int n, ParU_Control *Control);

ParU_Ret ParU_Backward(double *x1, double &resid, double &norm,
                       cholmod_sparse *A, ParU_Numeric *Num,
                       ParU_Control *Control);
#endif
