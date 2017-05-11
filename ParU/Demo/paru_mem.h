/* Wrappers for managing memory */
#include "../Include/Parallel_LU.hpp"
void *paralloc(int n, int size, cholmod_common* cc)
{
//    return malloc(n*size);
     return cholmod_l_malloc(n,size,cc);
    //return cholmod_l_calloc(n,size,cc);
}
//void *parfree(void *p)
//{
//    if (p) free(p);
//    return(NULL);
//}
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
