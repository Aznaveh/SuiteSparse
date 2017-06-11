/* Wrappers for managing memory */
#include "Parallel_LU.hpp"

void *paru_alloc(int n, int size, cholmod_common *cc){
     return cholmod_l_malloc(n,size,cc);        }

void *paru_calloc(int n, int size, cholmod_common *cc){
     return cholmod_l_calloc(n,size,cc);        }


void paru_free(int n, int size, void *p,  cholmod_common *cc){
    cholmod_l_free (n,   size, p, cc); }

void paru_freesym(paru_symbolic **LUsym_handle,
            // workspace and parameters
    cholmod_common *cc
){
    if (LUsym_handle == NULL || *LUsym_handle == NULL)
    {
        // nothing to do; caller probably ran out of memory
        return ;
    }

    paru_symbolic *LUsym;
    LUsym = *LUsym_handle;

    Long m, n, anz, nf; 

    m = LUsym->m;
    n = LUsym->n;
    nf = LUsym->nf; anz = LUsym->anz; 
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Super, cc);
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Parent, cc);
    cholmod_l_free (nf+2,   sizeof (SuiteSparse_long), LUsym->Childp, cc);
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Child, cc);
    cholmod_l_free (n+2,    sizeof (SuiteSparse_long), LUsym->Sleft, cc);
    cholmod_l_free (m+1,    sizeof (SuiteSparse_long), LUsym->Sp, cc);
    cholmod_l_free (anz,    sizeof (SuiteSparse_long), LUsym->Sj, cc);
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Fm, cc);
    cholmod_l_free (nf+1,   sizeof (SuiteSparse_long), LUsym->Cm, cc);

    cholmod_l_free (m+nf+1,   sizeof (SuiteSparse_long), LUsym->aParent, cc);
    cholmod_l_free (m+nf+1,   sizeof (SuiteSparse_long), LUsym->aChild, cc);
    cholmod_l_free (m+nf+2,   sizeof (SuiteSparse_long), LUsym->aChildp, cc);
    cholmod_l_free (m,   sizeof (SuiteSparse_long), LUsym->row2atree, cc);
    cholmod_l_free (nf,   sizeof (SuiteSparse_long), LUsym->super2atree, cc);

    cholmod_l_free (1, sizeof (paru_symbolic), LUsym, cc);

    *LUsym_handle = NULL;
}
/*! It uses LUsym, Do not free LUsym before*/
void paru_freemat(paru_matrix **paruMatInfo_handle,
    cholmod_common *cc
){
    if (paruMatInfo_handle == NULL || *paruMatInfo_handle == NULL ){
        return;
    }
    paru_matrix *paruMatInfo;
    paruMatInfo = *paruMatInfo_handle;

    Int m,n;  
    m = paruMatInfo->m;       n = paruMatInfo->n; 
    Int slackRow = paruMatInfo->slackRow;
    Int slackCol = paruMatInfo->slackCol;

    /*! TODO: I am changing tuple list allocation
     *  I have to change the free mechanism here too*/
    tupleList *RowList,*ColList;
    RowList= paruMatInfo->RowList; 
    cholmod_l_free (1, m*sizeof(tupleList), RowList, cc);
    ColList= paruMatInfo->ColList;
    cholmod_l_free (1, n*sizeof(tupleList), ColList, cc);

    paru_symbolic *LUsym = paruMatInfo-> LUsym;
    Element **elementList; 
    elementList = paruMatInfo->elementList;

    for(Int i = 0; i < m ; i++){        // freeing all elements
        Int e = LUsym->row2atree[i];    //element number in augmented tree
        Element *curEl = elementList[e];
        if (curEl == NULL) continue;
        Int nrows = curEl->nrows,
           ncols = curEl->ncols;
        cholmod_l_free (1, sizeof(Element)+sizeof(Int)*(nrows+ncols)+
                    sizeof(double)*nrows*ncols, curEl, cc);
    }
    Int nf = LUsym->nf;
    cholmod_l_free (1, (m+nf+1)*sizeof(Element), elementList, cc);
    *paruMatInfo_handle = NULL;
} 

