// ============================================================================/  
// ======================= Parallel_LU.hpp ====================================/
// ============================================================================/
#include <stdio.h>
#include "spqr.hpp"

// -----------------------------------------------------------------------------
// debugging and printing macros
// -----------------------------------------------------------------------------
#ifndef NPR
    #define NPR
#endif
//for printing information uncomment this
#undef NPR


#ifndef NPR
    static int print_level = 0 ;
    #define PRLEVEL(level,param) { if (print_level >= level) printf param ; }
    #define DEBUGLEVEL(level) { print_level = level ; }
#else
    #define PRLEVEL(level,param)
    #define DEBUGLEVEL(level)
#endif

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
typedef struct /* paru_symbolic*/
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

    Int m, n, anz ; // S is m-by-n with anz entries

    Int *Sp ;       // size m+1, row pointers of S

    Int *Sj ;       // size anz = Sp [n], column indices of S

    Int *Qfill ;    // size n, fill-reducing column permutation.
    // Qfill [k] = j if column k of A is column j of S.

    Int *PLinv ;    // size m, inverse row permutation that places
    // S=A(P,Q) in increasing order of leftmost column
    // index.  PLinv [i] = k if row i of A is row k of S.

    Int *Sleft ;    // size n+2.  The list of rows of S whose
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

    Int nf ;        // number of frontal matrices; nf <= MIN (m,n)
    Int maxfn ;     // max # of columns in any front

    // parent, child, and childp define the row merge tree or etree (A'A)
    Int *Parent ;   // size nf+1
    Int *Child ;    // size nf+1
    Int *Childp ;   // size nf+2

    // The parent of a front f is Parent [f], or EMPTY if f=nf.
    // A list of children of f can be obtained in the list
    // Child [Childp [f] ... Childp [f+1]-1].

    // Node nf in the tree is a placeholder; it does not represent a frontal
    // matrix.  All roots of the frontal "tree" (may be a forest) have the
    // placeholder node nf as their parent.  Thus, the tree of nodes 0:nf is
    // truly a tree, with just one parent (node nf).

    Int *aParent; // size m+nf
    Int *aChild; // size m+nf+1
    Int *aChildp; // size m+nf+2

    // pivot column in the front F.  This refers to a column of S.  The
    // number of expected pivot columns in F is thus
    // Super [f+1] - Super [f].

    //Upper bound number of rows for each front
    Int *Fm ;               // size nf+1

    //Upper bound  number of rows in the contribution block of each front
    Int *Cm ;               // size nf+1

    Int *Super ;    // size nf+1.  Super [f] gives the first
        // pivot column in the front F.  This refers to a column of S.  The
        // number of expected pivot columns in F is thus
        // Super [f+1] - Super [f].

    Int *Rp ;       // size nf+1
    Int *Rj ;       // size rjsize; compressed supernodal form of R

    Int rjsize ;    // size of Rj
 
    Int *row2atree;       //Mapping from rows to augmented tree size m
    Int *super2atree;     //Mapping from super nodes to augmented tree size nf

} paru_symbolic;

// =============================================================================
//      Tuple, Row and Column data structure 
// =============================================================================
typedef struct /* Tuple */
{
    /* The (e,f) tuples for element lists */
    int e,   /*  element number */
        f;  /*   offest */
} Tuple;

/* -------------------------------------------------------------------------- */
/* An element */
/* -------------------------------------------------------------------------- */

typedef struct	/* Element */
{
    Int

        nrowsleft,	/* number of rows remaining */
        ncolsleft,	/* number of columns remaining */
        nrows,		/* number of rows */
        ncols;		/* number of columns */

    /* followed in memory by:
       Int
       col [0..ncols-1],	column indices of this element
       row [0..nrows-1] ;	row indices of this element
       Entry			(suitably aligned, see macro below)
       C [0...nrows-1, 0...ncols-1] ;
       size of C is nrows*ncols Entry's
       */

} Element ; // contribution block


typedef struct  /*List of tuples */
{
    /*element of a column or a row*/
    Int
        numTuple,   /*  number of Tuples in this element */
        len;    /*  length of allocated space for current list*/
    Tuple *list;    /* list of tuples regarding to this element */

}   tupleList;



typedef struct  /*Matrix */
{
    Int m, n;
    paru_symbolic *LUsym;
    tupleList *RowList;
    tupleList *ColList;  
    Element **elementList; //pointers to all elements, size = m+nf+1 

}   paru_matrix;


paru_symbolic *paru_sym_analyse
( cholmod_sparse *A, cholmod_common *cc) ;

paru_matrix *paru_init_rowFronts 
(cholmod_sparse *A, paru_symbolic *LUsym, cholmod_common *cc);

/* Wrappers for managing memory */
void *paru_alloc(int n, int size, cholmod_common *cc);
void *paru_calloc(int n, int size, cholmod_common *cc);
void paru_free(int n, int size, void *p,  cholmod_common *cc);
void paru_freesym(paru_symbolic **LUsym_handle,cholmod_common *cc);
void paru_freemat(paru_matrix **paruMatInfo_handle, cholmod_common *cc);

/* add tuple functions defintions */
int paru_add_rowTuple (tupleList *RowList, Int row, Tuple T, 
        cholmod_common *cc);
int paru_add_colTuple (tupleList *ColList, Int col, 
        Tuple T, cholmod_common *cc);
void paru_assemble(paru_matrix *paruMatInfo, Int f);
