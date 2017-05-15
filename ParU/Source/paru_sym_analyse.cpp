// =============================================================================
// ===paru_sym_analyse============================================================
// =============================================================================

#include "Parallel_LU.hpp"


// =============================================================================

paru_symbolic * paru_sym_analyse
(
 // inputs, not modified
 cholmod_sparse *A,
 // workspace and parameters
 cholmod_common *cc
 )
{   
    paru_symbolic *LUsym;

    LUsym=(paru_symbolic *)paralloc(1,sizeof(paru_symbolic),cc);
    // ... check for LUsym NULL ...
    if(LUsym == NULL){
        //out of memory
        return NULL;
    }
        

    spqr_symbolic *QRsym;
    cc->SPQR_grain = 1;
    cc->useGPU = -1;
    QRsym = spqr_analyze (A, SPQR_ORDERING_CHOLMOD, FALSE,FALSE , FALSE, cc) ;

    Int m, n, anz,nf ; 
    m=LUsym->m=QRsym->m;
    n=LUsym->n=QRsym->n;
    anz=LUsym->anz=QRsym->anz;
    nf = LUsym->nf = QRsym->nf;

    LUsym->maxfn= LUsym->maxfn;

    Int *Parent, *Child, *Childp, *Rp, *ColCount, *Super;
    //brain transplant
    Parent= LUsym->Parent= QRsym->Parent; QRsym->Parent=NULL;
    Child=  LUsym->Child=  QRsym->Child;  QRsym->Child=NULL;
    Childp= LUsym->Childp= QRsym->Childp; QRsym->Childp=NULL;
    Super=  LUsym->Super=  QRsym->Super;  QRsym->Super=NULL;

    LUsym->Fm= QRsym->Fm; QRsym->Fm=NULL;
    LUsym->Cm= QRsym->Cm; QRsym->Cm=NULL;

    //Staircase structure
    Int *Sp, *Sj, *Sleft;
    Sp=    LUsym->Sp=    QRsym->Sp;     QRsym->Sp=NULL;
    Sj=    LUsym->Sj=    QRsym->Sj;     QRsym->Sj=NULL;
    Sleft= LUsym->Sleft= QRsym->Sleft;  QRsym->Sleft=NULL;


#ifndef NDEBUG
    /* print fronts*/
    for (Int f = 0; f < nf; f++) {
        Int fm, fn, fp;
        fm = LUsym->Fm[f];
        fn = QRsym->Rp[f+1]-QRsym->Rp[f];
        fp = Super[f+1]-Super[f];

        PR (("Front=%ld #col=%ld #row=%ld #pivotCol=%ld Par=%ld", 
                    f, fn, fm, fp,Parent[f]));
        PR ((" #pivot col= %ld",Super[f+1]-Super[f]));
        PR (("\nlist of children:\t"));
        for (Int i = Childp[f]; i <= Childp[f+1]-1; i++) 
            PR (("%ld ",Child[i]));
        PR (("\n\n"));
    }
#endif


    /*Computing augmented tree */
    Int *aParent; //augmented tree size m+nf+1
    Int *aChild;  // size m+nf+1
    Int *aChildp; // size m+nf+2
    aParent= (Int *) paralloc(m+nf+1, sizeof(Int),cc);
    aChild=  (Int *) paralloc(m+nf+1, sizeof(Int),cc);
    aChildp= (Int *) paralloc(m+nf+2, sizeof(Int),cc);

    Int *rM, *snM; // row map and supernode map
    rM=  (Int *) paralloc(m  ,sizeof(Int), cc);
    snM= (Int *) paralloc(nf ,sizeof(Int), cc);


    //initialization
    for (Int f = 0; f < nf; f++) snM[f]=-1;
    for (Int i = 0; i < m; i++) rM[i]=-1;
    for (Int i = 0; i < m+nf+1; i++) aParent[i]=-1;
    for (Int i = 0; i < m+nf+1; i++) aChild[i]=-1;
    for (Int i = 0; i < m+nf+2; i++) aChildp[i]=-1;

    aChildp[0]=0;
    Int offset=0; //number of rows visited in each iteration orig front+ rows
    Int lastChildFlag=0;
    Int childpointer=0;

    for (Int f = 0; f < nf; f++) {
        PR (("Front %ld\n", f)) ;
        PR (("pivot columns [ %ld to %ld ] n: %ld \n",
            Super [f], Super [f+1]-1, n)) ;
        ASSERT(Super[f+1] <= n);
        Int numRow =Sleft[Super[f+1]]-Sleft[Super[f]] ;
#ifndef NDEBUG
        PR (("~numRow=%ld",numRow));
        PR (("\n#offset=%ld\n",offset));
#endif

        Int numoforiginalChild=0;
        if (lastChildFlag){  // the current node is the parent
            PR (("Childs of %ld: ",f)) ;
            numoforiginalChild=Child[Childp[f+1]-1]-Child[Childp[f]]+1;
            for (Int i = Child[Childp[f]]; i < Child[Childp[f+1]]; i++){
                PR (("%ld,", i));
                ASSERT(snM[i] < m+nf+1);
                aParent[ snM[i]]=offset+numRow;
                ASSERT(childpointer < m+nf+1);
                aChild[childpointer++]=snM[i];
            }
        }

        for(Int i=offset ; i < offset+numRow ; i++)
            aChildp[i+1]=aChildp[i];

        for (Int i = Sleft[Super[f]]; i < Sleft[Super[f]+1]; i++){ // number of rows
            ASSERT(i < m);
            rM[i]=i+f;
            ASSERT(i+f < m+nf+1);
            aParent[i+f]=offset+numRow;
            ASSERT(childpointer < m+nf+1);
            aChild[childpointer++]=i+f;
        }

        offset+=numRow;
        snM[f]=offset++;
        ASSERT(offset < m+nf+1);
        aChildp[offset] = aChildp[offset-1]+numRow+numoforiginalChild;

        if( Parent[f] == f+1)  //last child
            lastChildFlag=1;  
        else
            lastChildFlag=0;  

    }

     spqr_freesym (&QRsym, cc) ;


    LUsym->aParent=     aParent;
    LUsym->aChildp=     aChildp;
    LUsym->aChild=      aChild;
    LUsym->row2atree=   rM;
    LUsym->super2atree= snM;

#ifndef NDEBUG
    PR (("\nsuper node->aP ")); 
    for (Int f = 0; f < nf; f++) PR (("%ld ",snM[f]));           PR (("\n"));
    PR (("row->aP "));  
    for (Int i = 0; i < m; i++)  PR (("%ld ",rM[i]));            PR (("\n"));
    PR (("aP: ")); 
    for (Int i = 0; i < m+nf+1; i++) PR (("%ld ",aParent[i]));   PR (("\n"));
    PR (("aChildp: "));
    for (Int i = 0; i < m+nf+1; i++) PR (("%ld ",aChildp[i]));   PR (("\n"));
    PR (("aChild: ")); 
    for (Int i = 0; i < m+nf+1; i++) PR (("%ld ",aChild[i]));    PR (("\n"));

    for(Int i=0; i< m+nf; i++){
        PR (("anode:%ld",i));
        for(Int c=aChildp[i]; c< aChildp[i+1]; c++)
            PR ((" %ld,",aChild[c]));                            PR (("\n"));
    }
#endif

    return (LUsym) ;
}
