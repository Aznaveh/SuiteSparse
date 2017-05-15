/* Wrappers for managing memory */
//#include "Parallel_LU.hpp"
#include "spqr.hpp"
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

void *paralloc(int n, int size, cholmod_common* cc)
{
//    return malloc(n*size);
     return cholmod_l_malloc(n,size,cc);
    //return cholmod_l_calloc(n,size,cc);
}

void paru_freesym(paru_symbolic** LUsym_handle,
            // workspace and parameters
    cholmod_common *cc
)
{
    //uncomplete

    if (LUsym_handle == NULL || *LUsym_handle == NULL)
    {
        // nothing to do; caller probably ran out of memory
        return ;
    }

    paru_symbolic *LUsym ;
    LUsym = *LUsym_handle  ;

    Long m, n, anz, nf; 

    m = LUsym->m ;
    n = LUsym->n ;
    nf = LUsym->nf ; anz = LUsym->anz ; 
    //PR((">><<>><<%d\n",LUsym->Super));
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Super, cc) ;
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Parent, cc) ;
    cholmod_l_free (nf+2,   sizeof (SuiteSparse_long), LUsym->Childp, cc) ;
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Child, cc) ;
    cholmod_l_free (n+2,    sizeof (SuiteSparse_long), LUsym->Sleft, cc) ;
    cholmod_l_free (m+1,    sizeof (SuiteSparse_long), LUsym->Sp, cc) ;
    cholmod_l_free (anz,    sizeof (SuiteSparse_long), LUsym->Sj, cc) ;
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Fm, cc) ;
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Cm, cc) ;

    cholmod_l_free (m+nf+1,   sizeof (SuiteSparse_long), LUsym->aParent, cc) ;
    cholmod_l_free (m+nf+1,   sizeof (SuiteSparse_long), LUsym->aChild, cc) ;
    cholmod_l_free (m+nf+2,   sizeof (SuiteSparse_long), LUsym->aChildp, cc) ;
    cholmod_l_free (m,   sizeof (SuiteSparse_long), LUsym->row2atree, cc) ;
    cholmod_l_free (nf,   sizeof (SuiteSparse_long), LUsym->super2atree, cc) ;


   /* 
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Rp, cc) ;
    cholmod_l_free (rjsize, sizeof (SuiteSparse_long), LUsym->Rj, cc) ;
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Post, cc) ;
    cholmod_l_free (m,      sizeof (SuiteSparse_long), LUsym->PLinv, cc) ;
    cholmod_l_free (n,      sizeof (SuiteSparse_long), LUsym->ColCount, cc) ;
    */
//    parfree((paru_symbolic *)LUsym->Sp);
//    parfree((paru_symbolic *)LUsym->Sj);
//    parfree((paru_symbolic *)LUsym->Sleft);
//    parfree((paru_symbolic *)LUsym->Parent);
//    parfree((paru_symbolic *)LUsym->Child);
//    parfree((paru_symbolic *)LUsym->Childp);
//    parfree((paru_symbolic *)LUsym->aParent);
//    parfree((paru_symbolic *)LUsym->aChild);
//    parfree((paru_symbolic *)LUsym->aChildp);
//    parfree((paru_symbolic *)LUsym->Super);
//    parfree((paru_symbolic *)LUsym->Fm);
//    parfree((paru_symbolic *)LUsym->Cm);
//    parfree((paru_symbolic *)LUsym->row2atree);
//    parfree((paru_symbolic *)LUsym->super2atree);
 
    cholmod_l_free (1, sizeof (paru_symbolic), LUsym, cc) ;

    *LUsym_handle = NULL ;
}
