/* Wrappers for managing memory */
#include "Parallel_LU.hpp"
void *paralloc(int n, int size, cholmod_common* cc)
{
     return cholmod_l_malloc(n,size,cc);
}

void paru_freesym(paru_symbolic** LUsym_handle,
            // workspace and parameters
    cholmod_common *cc
)
{

    DEBUGLEVEL(1);
//    if (LUsym_handle == NULL || *LUsym_handle == NULL)
//    {
//        // nothing to do; caller probably ran out of memory
//        return ;
//    }
//
//    PRLEVEL (1, ("In free sym\n"));
//    
//    
//    paru_symbolic *LUsym ;
//    LUsym = *LUsym_handle  ;
//
//    Long m, n, anz, nf; 
//
//    m = LUsym->m ;
//    n = LUsym->n ;
//    nf = LUsym->nf ; anz = LUsym->anz ; 
//    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Super, cc) ;
//    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Parent, cc) ;
//    cholmod_l_free (nf+2,   sizeof (SuiteSparse_long), LUsym->Childp, cc) ;
//    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Child, cc) ;
//    cholmod_l_free (n+2,    sizeof (SuiteSparse_long), LUsym->Sleft, cc) ;
//    cholmod_l_free (m+1,    sizeof (SuiteSparse_long), LUsym->Sp, cc) ;
//    cholmod_l_free (anz,    sizeof (SuiteSparse_long), LUsym->Sj, cc) ;
//    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Fm, cc) ;
//    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Cm, cc) ;
//
//    cholmod_l_free (m+nf+1,   sizeof (SuiteSparse_long), LUsym->aParent, cc) ;
//    cholmod_l_free (m+nf+1,   sizeof (SuiteSparse_long), LUsym->aChild, cc) ;
//    cholmod_l_free (m+nf+2,   sizeof (SuiteSparse_long), LUsym->aChildp, cc) ;
//    cholmod_l_free (m,   sizeof (SuiteSparse_long), LUsym->row2atree, cc) ;
//    cholmod_l_free (nf,   sizeof (SuiteSparse_long), LUsym->super2atree, cc) ;
//
//    cholmod_l_free (1, sizeof (paru_symbolic), LUsym, cc) ;
//
    *LUsym_handle = NULL ;
}
