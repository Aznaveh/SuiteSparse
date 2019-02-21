/** Wrappers for managing memory 
 * allocating and freeing is done through cholmod and these wrappers
 * */
#include "Parallel_LU.hpp"

void *paru_alloc (Int n, Int size, cholmod_common *cc){
    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int alloc_count =0;
    alloc_count += n*size ;
#endif
    void *p = cholmod_l_malloc (n,size,cc);
    PRLEVEL (1, ("%% allocated %ld in %p total= %ld\n", 
                n*size, p, alloc_count ));
    return p;
}

void *paru_calloc(Int n, Int size, cholmod_common *cc){
    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int calloc_count =0;
    calloc_count += n*size ;
#endif
    void *p= cholmod_l_calloc (n,size,cc);       
    PRLEVEL (1, ("%% callocated %ld in %p total= %ld\n", 
                n*size, p, calloc_count ));
    return p;
}

void *paru_realloc(
        Int newsize,    // requested size
        Int size_Entry, // size of each Entry
        void *oldP,     // pointer to the old allocated space
        Int *size,  // a single number, input: old size, output: new size
        cholmod_common *cc){

    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int realloc_count =0;
    realloc_count += newsize*size_Entry - *size;
#endif
    void *p= cholmod_l_realloc (newsize, size_Entry, oldP, (size_t*)size, cc);
    PRLEVEL (1, ("%% reallocated %ld in %p and freed %p total= %ld\n", 
                newsize*size_Entry, p, oldP, realloc_count ));
    return p;
}


void paru_free (Int n, Int size, void *p,  cholmod_common *cc){
    DEBUGLEVEL(0);
#ifndef NDEBUG
    static Int free_count =0;
    free_count+= n*size ;
#endif
    if(p != NULL)
        cholmod_l_free (n,   size, p, cc);
    else
        printf("freeing a NULL pointer\n");
    PRLEVEL (1, ("%% free %ld in %p total= %ld\n", 
                n*size, p, free_count));

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
    PRLEVEL (1, ("%% In free sym: m=%ld n=%ld\n nf=%ld\
                anz=%ld rjsize=%ld\n", m, n, nf, anz, rjsize ));

    paru_free (nf+1, sizeof (Int), LUsym->Parent, cc);
    paru_free (nf+1, sizeof (Int), LUsym->Child, cc);
    paru_free (nf+2, sizeof (Int), LUsym->Childp, cc);
    paru_free (nf+1, sizeof (Int), LUsym->Super, cc);
    paru_free (n, sizeof (Int), LUsym->Qfill , cc);
    paru_free (m, sizeof (Int), LUsym->PLinv, cc);
    paru_free (nf+1, sizeof (Int), LUsym->Fm, cc);
    paru_free (nf+1, sizeof (Int), LUsym->Cm, cc);

    paru_free (rjsize, sizeof (Int), LUsym->Rj, cc);
    paru_free (nf+1,   sizeof (Int), LUsym->Rp, cc);

    paru_free (m+1, sizeof (Int), LUsym->Sp, cc);
    paru_free (anz, sizeof (Int), LUsym->Sj, cc);
    paru_free (n+2, sizeof (Int), LUsym->Sleft, cc);

    paru_free (m+nf, sizeof (Int), LUsym->aParent, cc);
    paru_free (m+nf+1, sizeof (Int), LUsym->aChild, cc);
    paru_free (m+nf+2, sizeof (Int), LUsym->aChildp, cc);
    paru_free (m, sizeof (Int), LUsym->row2atree, cc);
    paru_free (nf, sizeof (Int), LUsym->super2atree, cc);

    paru_free (1, sizeof (paru_symbolic), LUsym, cc);

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
    PRLEVEL (1, ("%% RowList =%p\n", RowList));
    tupleList *ColList = paruMatInfo->ColList;
    PRLEVEL (1, ("%% ColList =%p\n", ColList));


    paru_symbolic *LUsym = paruMatInfo-> LUsym;
    Int nf = LUsym->nf;

    // free tuple lists 
    for (Int col = 0; col < n; col++) {
        Int len = ColList [col].len;
        // ASSERT (len < m);  // it is a wrong assertion but there is a good
        //  point
        if (len > m+nf )                        
            printf ("%% too much space used for %ld\n",col);
        paru_free (len , sizeof (Tuple), ColList[col].list, cc);
    }
    paru_free (1, n*sizeof(tupleList), ColList, cc);

    for (Int row = 0; row < m; row++) {
        Int len = RowList [row].len;
        paru_free (len , sizeof (Tuple), RowList[row].list, cc);
    }
    paru_free (1, m*sizeof(tupleList), RowList, cc);


    Element **elementList; 
    elementList = paruMatInfo->elementList;


    PRLEVEL (1, ("%% LUsym = %p\n",LUsym));
    PRLEVEL (1, ("%% freeing initialized elements:\n"));
    for(Int i = 0; i < m ; i++){        // freeing all row elements
        if(LUsym == NULL){
            printf ("Probably LUsym has been freed before! Wrong usage\n");
            return;
        }
        Int e = LUsym->row2atree[i];    //element number in augmented tree
        PRLEVEL (1, ("%% e =%ld\t", e));
        Element *curEl = elementList[e];
        if (curEl == NULL) continue;
        Int nrows = curEl->nrows,
            ncols = curEl->ncols;
        PRLEVEL (1, ("%% nrows =%ld ", nrows));
        PRLEVEL (1, ("%% ncols =%ld\n", ncols));
        paru_free (1, sizeof(Element)+sizeof(Int)*(2*(nrows+ncols)+2)+
                sizeof(double)*nrows*ncols, curEl, cc);
    }


    PRLEVEL (1, ("%% freeing CB elements:\n"));
    for(Int i = 0; i < nf ; i++){        // freeing all other elements
        Int e = LUsym->super2atree[i];    //element number in augmented tree
        PRLEVEL (1, ("%% e =%ld\t", e));
        Element *curEl = elementList[e];
        if (curEl == NULL) continue; /* CB not used */
        Int nrows = curEl->nrows,
            ncols = curEl->ncols;
        Int tot_size = sizeof(Element)+sizeof(Int)*(2*(nrows+ncols))+
            sizeof(double)*nrows*ncols;
        PRLEVEL (1, ("%% nrows =%ld ", nrows));
        PRLEVEL (1, ("%% ncols =%ld tot_size=%ld\n", ncols, tot_size));
       paru_free (1, tot_size, curEl, cc);
    }

    //free the answer
    paru_fac *LUs =  paruMatInfo->partial_LUs;
    paru_fac *Us =  paruMatInfo->partial_Us;
    for(Int i = 0; i < nf ; i++){  
        PRLEVEL (1, ("%% Freeing Us=%p\n", Us[i].p));
        if(Us[i].p != NULL){
            Int m=Us[i].m; Int n=Us[i].n;
            paru_free (m*n, sizeof (double), Us[i].p, cc);
        }
        PRLEVEL (1, ("%% Freeing LUs=%p\n", LUs[i].p));
        if(LUs[i].p != NULL){
            Int m=LUs[i].m; Int n=LUs[i].n;
            paru_free (m*n, sizeof (double), LUs[i].p, cc);
        }
    }

    paru_free(1, nf*sizeof(paru_fac),LUs, cc);
    paru_free(1, nf*sizeof(paru_fac),Us, cc);
    paru_free(1, nf*sizeof(Int),paruMatInfo->time_stamp, cc);


    paru_free (1, (m+nf+1)*sizeof(Element), elementList, cc);
    work_struct *Work = paruMatInfo->Work;
    paru_free (m, sizeof(Int), Work->rowSize, cc);
    paru_free (2*m+n, sizeof(Int), Work->scratch, cc);
    paru_free (n, sizeof(Int), Work->colSize, cc);
    paru_free (m+nf, sizeof(Int), Work->elRow, cc);
    paru_free (m+nf, sizeof(Int), Work->elCol, cc);


    paru_free (1, sizeof(work_struct), paruMatInfo->Work, cc);
    paru_free (1, sizeof(paru_matrix), paruMatInfo, cc);
    *paruMatInfo_handle = NULL;
} 
