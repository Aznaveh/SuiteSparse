// =============================================================================
// === Parallel_LU.hpp =======================================================
// =============================================================================

#include "spqr.hpp"

#define Int SuiteSparse_long

// =============================================================================
// === paru_symbolic ===========================================================
// =============================================================================

// The contents of this object do not change during numeric factorization.  The
// Symbolic object depends only on the pattern of the input matrix, and not its
// values. 
// This makes parallelism easier to manage, since all threads can
// have access to this object without synchronization.
//
struct paru_symbolic
{

    // -------------------------------------------------------------------------
    // row-form of the input matrix and its permutations
    // -------------------------------------------------------------------------

    // During symbolic analysis, the nonzero pattern of S = A(P,Q) is
    // constructed, where A is the user's input matrix.  Its numerical values
    // are also constructed, but they do not become part of the Symbolic
    // object.  The matrix S is stored in row-oriented form.  The rows of S are
    // sorted according to their leftmost column index (via PLinv).  Column
    // indices in each row of S are in strictly ascending order, even though
    // the input matrix A need not be sorted.

    SuiteSparse_long m, n, anz ; // S is m-by-n with anz entries

    SuiteSparse_long *Sp ;       // size m+1, row pointers of S

    SuiteSparse_long *Sj ;       // size anz = Sp [n], column indices of S

//    SuiteSparse_long *Qfill ;    // size n, fill-reducing column permutation.
                        // Qfill [k] = j if column k of A is column j of S.

//    SuiteSparse_long *PLinv ;    // size m, inverse row permutation that places
                        // S=A(P,Q) in increasing order of leftmost column
                        // index.  PLinv [i] = k if row i of A is row k of S.

    SuiteSparse_long *Sleft ;    // size n+2.  The list of rows of S whose
            // leftmost column index is j is given by
            // Sleft [j] ... Sleft [j+1]-1.  This can be empty (that is, Sleft
            // [j] can equal Sleft [j+1]).  Sleft [n] is the number of
            // non-empty rows of S, and Sleft [n+1] == m.  That is, Sleft [n]
            // ... Sleft [n+1]-1 gives the empty rows of S, if any.

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

    SuiteSparse_long nf ;        // number of frontal matrices; nf <= MIN (m,n)
    SuiteSparse_long maxfn ;     // max # of columns in any front

    // parent, child, and childp define the row merge tree or etree (A'A)
    SuiteSparse_long *Parent ;   // size nf+1
    SuiteSparse_long *Child ;    // size nf+1
    SuiteSparse_long *Childp ;   // size nf+2

    // The parent of a front f is Parent [f], or EMPTY if f=nf.
    // A list of children of f can be obtained in the list
    // Child [Childp [f] ... Childp [f+1]-1].

    // Node nf in the tree is a placeholder; it does not represent a frontal
    // matrix.  All roots of the frontal "tree" (may be a forest) have the
    // placeholder node nf as their parent.  Thus, the tree of nodes 0:nf is
    // truly a tree, with just one parent (node nf).

    SuiteSparse_long *aParent; // size m+nf+1
    SuiteSparse_long *aChild; // size m+nf+1
    SuiteSparse_long *aChildp; // size m+nf+2

    SuiteSparse_long *Super ;    // size nf+1.  Super [f] gives the first
        // pivot column in the front F.  This refers to a column of S.  The
        // number of expected pivot columns in F is thus
        // Super [f+1] - Super [f].

   //Upper bound number of rows for each front
    SuiteSparse_long *Fm ;               // size nf+1

    //Upper bound  number of rows in the contribution block of each front
    SuiteSparse_long *Cm ;               // size nf+1

    SuiteSparse_long *row2atree;               //Mapping from rows to augmented tree size m
    SuiteSparse_long *super2atree;               //Mapping from super nodes to augmented tree size m

} ;

/* Wrappers for managing memory */
void *paralloc(int n, int size, cholmod_common* cc);
void paru_freesym(paru_symbolic** LUsym_handle,cholmod_common *cc);


paru_symbolic *paru_sym_analyse
( cholmod_sparse *A, cholmod_common *cc) ;

// =============================================================================
//      Tuple, Row and Column data structure 
// =============================================================================
typedef struct /* Tuple */
{
    /* The (e,f) tuples for element lists */
    int e,   /*  element number */
        f;  /*   offest */
} Tuple;

typedef struct /* Element */
{
    int 
        nrows,
        ncols;
    /* followed by 
     * int col[0..ncols-1], //column indices of the element
     * int row[0..nrows-1]; row indices of the element
     * double C[0..nrows*ncols-1] * Numerical values
     */
} Element;

