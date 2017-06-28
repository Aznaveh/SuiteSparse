/* =========================================================================   /
 * ============================== paru_sym_analyse =========================   /
 * =========================================================================   /
 * Computing augmented tree and using QRsym data structure and add it to LUsym
 *
 * Example: ./Matrix/problem.mtx
 * original post ordered etree:          augmented tree: (Row elements)
 *
 *   4___                                     16_____
 *   |   \                                    |   \  (14-15) 
 *   0    3__                                 3    13_____
 *        |  \                              (0-2)  |  \   \
 *        1   2                                    7   10 (11-12)
 *                                                 |    \
 *                                               (4-6) (8-9) 
 *                                               
 *      Note: the original etree use a dummy parent for all the tree(forest)
 *              the augmented tree does not
 * */
#include "Parallel_LU.hpp"
paru_symbolic *paru_sym_analyse
(
 // inputs, not modified
 cholmod_sparse *A,
 // workspace and parameters
 cholmod_common *cc ){   

    DEBUGLEVEL(0);
    paru_symbolic *LUsym;

    LUsym = (paru_symbolic*) paru_alloc (1, sizeof(paru_symbolic), cc);
    // ... check for LUsym NULL ...
    if(LUsym == NULL){
        //out of memory
        return NULL;
    }

    spqr_symbolic *QRsym;
    cc->SPQR_grain = 1;
    cc->useGPU = -1;
    QRsym = spqr_analyze (A, SPQR_ORDERING_CHOLMOD, FALSE,FALSE , FALSE, cc);

    Int m, n, anz,nf, rjsize ; 
    m = LUsym->m = QRsym->m;
    n = LUsym->n = QRsym->n;
    anz = LUsym->anz = QRsym->anz;
    nf =  LUsym->nf = QRsym->nf;
    rjsize =  LUsym->rjsize = QRsym->rjsize;


    PRLEVEL (0, ("A  is  %ld x %ld \n",m, n ));
    PRLEVEL (0, ("nf=%ld\n",nf ));
    PRLEVEL (0, ("anz = %ld  rjsize=%ld\n", anz, rjsize));


    LUsym->maxfn = QRsym->maxfn;

    Int *Parent, *Child, *Childp, 
        *Super, *Qfill, *PLinv;
    //brain transplant
    Parent = LUsym->Parent = QRsym->Parent;
    QRsym->Parent = NULL;
    Child =  LUsym->Child =  QRsym->Child;  
    QRsym->Child = NULL;
    Childp = LUsym->Childp = QRsym->Childp; 
    QRsym->Childp = NULL;
    Super =  LUsym->Super =  QRsym->Super;  
    QRsym->Super = NULL;
    Qfill =  LUsym->Qfill =  QRsym->Qfill;  
    QRsym->Qfill = NULL;
    PLinv =  LUsym->PLinv =  QRsym->PLinv;  
    QRsym->PLinv = NULL;
    LUsym->Fm = QRsym->Fm; 
    QRsym->Fm = NULL;
    LUsym->Cm = QRsym->Cm; 
    QRsym->Cm = NULL;
    LUsym->Rj = QRsym->Rj; 
    QRsym->Rj = NULL;
    LUsym->Rp = QRsym->Rp; 
    QRsym->Rp = NULL;

    //Staircase structure
    Int *Sp, *Sj, *Sleft;
    Sp =  LUsym->Sp = QRsym->Sp;     
    QRsym->Sp = NULL;
    Sj =  LUsym->Sj = QRsym->Sj;     
    QRsym->Sj = NULL;
    Sleft = LUsym->Sleft = QRsym->Sleft;  
    QRsym->Sleft = NULL;

    spqr_freesym (&QRsym, cc); //No longer needed

    /*Computing augmented tree */
    Int *aParent; //augmented tree size m+nf
    Int *aChild;  // size m+nf+1
    Int *aChildp; // size m+nf+2
    Int *rM, *snM; // row map and supernode map

    //initializing with NULL to avoid freeing not allocated memory
    aParent = LUsym->aParent = NULL; 
    aChildp = LUsym->aChildp = NULL;
    aChild = LUsym->aChild = NULL;
    rM = LUsym->row2atree = NULL;
    snM = LUsym->super2atree = NULL;

    /*! Check if there exist empty row or column in square part*/
   for (Int row = 0; row < m; row++){
        PRLEVEL (1,("Sprow[%ld]=%ld\n", row, Sp[row]));
        if (Sp [row] == Sp[row+1] ){
            printf("Empty Row\n");
            paru_freesym (&LUsym , cc);
            return NULL;
        }
    }
    for (Int col = 0; col < n; col++){
        Int *Ap =(Int*) A->p;
        if (Ap [col] == Ap [col+1]){
            printf("Empty Column\n");
            paru_freesym (&LUsym , cc);
            return NULL;
        }
    }


    //#ifndef NDEBUG
    /* print fronts*/
    for (Int f = 0; f < nf; f++){
        Int fm, fn, fp;
        fm = LUsym->Fm[f];
        fn = LUsym->Rp[f+1]-LUsym->Rp[f];
        fp = Super[f+1]-Super[f];

        PRLEVEL (0,("Front=%ld Parent=%ld\n", f,  Parent[f]));
        PRLEVEL (1,("\nlist of children: "));
        for (Int i = Childp[f]; i <= Childp[f+1]-1; i++) 
            PRLEVEL (1,("%ld ",Child[i]));
        PRLEVEL (1,("\n"));
    }
    //#endif /* end of NDEBUG */

    aParent = (Int*) paru_alloc (m+nf, sizeof(Int),cc);
    aChild =  (Int*) paru_alloc (m+nf+1, sizeof(Int),cc);
    aChildp = (Int*) paru_alloc (m+nf+2, sizeof(Int),cc);


    rM =  (Int*) paru_alloc(m  ,sizeof(Int), cc);
    snM = (Int*) paru_alloc(nf ,sizeof(Int), cc);

    if(aParent == NULL || aChild == NULL || aChildp == NULL ||
            rM == NULL  || snM == NULL ){
        printf ("Out of memory");
        paru_freesym (&LUsym , cc);
        return NULL;
    }

    //initialization
    /*! TODO: Not all of them need initialization this is for debug for now
     * aParent needs initialization*/
    for (Int f = 0; f < nf; f++) snM[f] = -1;
    for (Int i = 0; i < m; i++) rM[i] = -1;
    for (Int i = 0; i < m+nf; i++) aParent[i] = -1;
    for (Int i = 0; i < m+nf+1; i++) aChild[i] = -1;
    for (Int i = 0; i < m+nf+2; i++) aChildp[i] = -1;

    aChildp[0] = 0;
    Int offset = 0; //number of rows visited in each iteration orig front+ rows
    Int lastChildFlag = 0;
    Int childpointer = 0;

    for (Int f = 0; f < nf; f++) {
        PRLEVEL (1,("Front %ld\n", f)) ;
        PRLEVEL (1,("pivot columns [ %ld to %ld ] n: %ld \n",
                    Super [f], Super [f+1]-1, n)) ;
        ASSERT(Super[f+1] <= n);
        Int numRow =Sleft[Super[f+1]]-Sleft[Super[f]] ;

        Int numoforiginalChild=0;
        if (lastChildFlag){  // the current node is the parent
            PRLEVEL (1,("Childs of %ld: ",f)) ;
            numoforiginalChild= Childp[f+1] - Childp[f];

            for (Int i = Childp[f]; i <= Childp[f+1]-1; i++) {
                PRLEVEL (1,("%ld ",Child[i]));
                Int c= Child [i];
                ASSERT(snM[c] < m+nf+1);
                aParent[ snM[c]]=offset+numRow;
                PRLEVEL (1, ("aParent[%ld] =%ld\n", 
                            aParent[snM [c]], offset+numRow));
                ASSERT(childpointer < m+nf+1);
                aChild[childpointer++] = snM[c];
            }
        }

        PRLEVEL (1,("numRow=%ld ",numRow));
        PRLEVEL (1,("#offset=%ld\n",offset));
        for(Int i = offset ; i < offset+numRow ; i++){
            ASSERT (aChildp [i+1] == -1);
            aChildp[i+1] = aChildp[i];
            PRLEVEL (1, ("@i=%ld aCp[%ld]=%ld aCp[%ld]=%ld",
                        i, i+1, aChildp[i+1], i, aChildp[i]));
        }

        for (Int i = Sleft[Super[f]]; i < Sleft[Super[f+1]]; i++){ 
            // number of rows
            ASSERT(i < m);

            rM[i] = i+f;
            ASSERT(i+f < m+nf+1);
            aParent[i+f] = offset+numRow;
            ASSERT(childpointer < m+nf+1);
            aChild[childpointer++] = i+f;
        }

        offset += numRow;
        snM[f] = offset++;
        ASSERT(offset < m+nf+1);
        ASSERT (aChildp [offset] == -1);
        aChildp[offset] = aChildp[offset-1]+numRow+numoforiginalChild;
        PRLEVEL (1, ("\n f=%ld numoforiginalChild=%ld\n", f, numoforiginalChild));

        if( Parent[f] == f+1){  //last child due to staircase
            PRLEVEL (1, ("last Child =%ld\n", f));
            lastChildFlag = 1;  
        }else
            lastChildFlag = 0;  
    }

    LUsym->aParent = aParent;
    LUsym->aChildp = aChildp;
    LUsym->aChild = aChild;
    LUsym->row2atree = rM;
    LUsym->super2atree = snM;

#ifndef NDEBUG
    PRLEVEL (1,("super node mapping (snM): ")); 
    for (Int f = 0; f < nf; f++){
        ASSERT (snM [f] != -1);
        PRLEVEL (1,("%ld ",snM[f]));     
    }
    PRLEVEL (1,("\n"));

    PRLEVEL (1,("row mapping (rM): "));  
    for (Int i = 0; i < m; i++){
        ASSERT (rM [i] != -1);
        PRLEVEL (1,("%ld ",rM[i]));        
    }
    PRLEVEL (1,("\n"));

    PRLEVEL (1,("aParent: ")); 
    for (Int i = 0; i < m+nf; i++) PRLEVEL (1,("%ld ",aParent[i]));   
    PRLEVEL (1,("\n"));

    PRLEVEL (1,("aChildp: "));
    for (Int i = 0; i < m+nf+2; i++) PRLEVEL (1,("%ld ",aChildp[i]));   
    PRLEVEL (1,("\n"));

    PRLEVEL (1,("aChild: ")); 
    for (Int i = 0; i < m+nf+1; i++) PRLEVEL (1,("%ld ",aChild[i]));    
    PRLEVEL (1,("\n"));

    for(Int i=0; i< m+nf; i++){
        PRLEVEL (1,("anode:%ld",i));
        for(Int c=aChildp[i]; c< aChildp[i+1]; c++)
            PRLEVEL (1,(" %ld,",aChild[c]));  
        PRLEVEL (1,("\n"));
    }
#endif
    return (LUsym) ;
}
