/** Wrappers for managing memory 
 * allocating and freeing is done through cholmod and these wrappers
 * */
#include "Parallel_LU.hpp"

void *paru_alloc (Int n, Int size, cholmod_common *cc){
    return cholmod_l_malloc (n,size,cc);
}

void *paru_calloc(Int n, Int size, cholmod_common *cc){
    return cholmod_l_calloc (n,size,cc);       
}

void *paru_realloc(
        Int newsize,    // requested size
        Int size_Entry, // size of each Entry
        void *oldP,     // pointer to the old allocated space
        Int *size,  // a single number, input: old size, output: new size
        cholmod_common *cc){
    return cholmod_l_realloc (newsize, size_Entry, oldP, (size_t*)size, cc);
}



void paru_free (Int n, Int size, void *p,  cholmod_common *cc){
    cholmod_l_free (n,   size, p, cc);
}

void paru_freesym (paru_symbolic **LUsym_handle,
        // workspace and parameters
        cholmod_common *cc
        ){
    DEBUGLEVEL (0);
    if (LUsym_handle == NULL || *LUsym_handle == NULL){
        // nothing to do; caller probably ran out of memory
        return;
    }

    paru_symbolic *LUsym;
    LUsym = *LUsym_handle;

    Int m, n, anz, nf, rjsize; 

    m = LUsym->m;
    n = LUsym->n;
    nf = LUsym->nf; 
    anz = LUsym->anz; 
    rjsize = LUsym->rjsize;
    PRLEVEL (1, ("In free sym: m=%ld n=%ld\n nf=%ld\
                anz=%ld rjsize=%ld\n", m, n, nf, anz, rjsize ));

    cholmod_l_free (nf+1, sizeof (Int), LUsym->Parent, cc);
    cholmod_l_free (nf+1, sizeof (Int), LUsym->Child, cc);
    cholmod_l_free (nf+2, sizeof (Int), LUsym->Childp, cc);
    cholmod_l_free (nf+1, sizeof (Int), LUsym->Super, cc);
    cholmod_l_free (n, sizeof (Int), LUsym->Qfill , cc);
    cholmod_l_free (m, sizeof (Int), LUsym->PLinv, cc);
    cholmod_l_free (nf+1, sizeof (Int), LUsym->Fm, cc);
    cholmod_l_free (nf+1, sizeof (Int), LUsym->Cm, cc);

    cholmod_l_free (rjsize, sizeof (Int), LUsym->Rj, cc);
    cholmod_l_free (nf+1,   sizeof (Int), LUsym->Rp, cc);

    cholmod_l_free (m+1, sizeof (Int), LUsym->Sp, cc);
    cholmod_l_free (anz, sizeof (Int), LUsym->Sj, cc);
    cholmod_l_free (n+2, sizeof (Int), LUsym->Sleft, cc);


    cholmod_l_free (m+nf, sizeof (Int), LUsym->aParent, cc);
    cholmod_l_free (m+nf+1, sizeof (Int), LUsym->aChild, cc);
    cholmod_l_free (m+nf+2, sizeof (Int), LUsym->aChildp, cc);
    cholmod_l_free (m, sizeof (Int), LUsym->row2atree, cc);
    cholmod_l_free (nf, sizeof (Int), LUsym->super2atree, cc);

    cholmod_l_free (1, sizeof (paru_symbolic), LUsym, cc);

    *LUsym_handle = NULL;
}

/*! It uses LUsym, Do not free LUsym before*/
void paru_freemat (paru_matrix **paruMatInfo_handle,
        cholmod_common *cc){

    DEBUGLEVEL(0); 
    if (paruMatInfo_handle == NULL || *paruMatInfo_handle == NULL ){
        return;
    }
    paru_matrix *paruMatInfo;
    paruMatInfo = *paruMatInfo_handle;

    Int m,n;  
    m = paruMatInfo->m;       
    n = paruMatInfo->n; 

    tupleList *RowList = paruMatInfo->RowList;
    PRLEVEL (2, ("RowList =%p\n", RowList));
    tupleList *ColList = paruMatInfo->ColList;
    PRLEVEL (2, ("ColList =%p\n", ColList));

    // free tuple lists 
    for (Int col = 0; col < n; col++) {
        Int len = ColList [col].len;
        cholmod_l_free (len , sizeof (Tuple), ColList[col].list, cc);
    }
    cholmod_l_free (1, n*sizeof(tupleList), ColList, cc);

    for (Int row = 0; row < m; row++) {
        Int len = RowList [row].len;
        cholmod_l_free (len , sizeof (Tuple), RowList[row].list, cc);
    }
    cholmod_l_free (1, m*sizeof(tupleList), RowList, cc);


    paru_symbolic *LUsym = paruMatInfo-> LUsym;
    Element **elementList; 
    elementList = paruMatInfo->elementList;

    /*! TODO: This code just work for row initialized elements	
     * this code must be more general
     *
     * I should free other elements later 
     * */
    for(Int i = 0; i < m ; i++){        // freeing all row elements
        if(LUsym == NULL){
            printf ("Probably LUsym has been freed before! Wrong usage\n");
            return;
        }
        PRLEVEL (1, ("LUsym = %p\n",LUsym));
        Int e = LUsym->row2atree[i];    //element number in augmented tree
        PRLEVEL (1, ("e =%ld\t", e));
        Element *curEl = elementList[e];
        if (curEl == NULL) continue;
        Int nrows = curEl->nrows,
            ncols = curEl->ncols;
        PRLEVEL (1, ("nrows =%ld ", nrows));
        PRLEVEL (1, ("ncols =%ld\n", ncols));
        cholmod_l_free (1, sizeof(Element)+sizeof(Int)*(nrows+ncols)+
                sizeof(double)*nrows*ncols, curEl, cc);
    }

    Int nf = LUsym->nf;
    cholmod_l_free (1, (m+nf+1)*sizeof(Element), elementList, cc);
    work_struct *Work = paruMatInfo->Work;
    cholmod_l_free (m, sizeof(Int), Work->all_Zero, cc);
    cholmod_l_free (m, sizeof(Int), Work->scratch, cc);
    
    cholmod_l_free (1, sizeof(work_struct), paruMatInfo->Work, cc);
    cholmod_l_free (1, sizeof(paru_matrix), paruMatInfo, cc);
    *paruMatInfo_handle = NULL;
} 
