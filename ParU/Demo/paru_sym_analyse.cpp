// =============================================================================
// ===paru_sym_analyse============================================================
// =============================================================================
#include "../../SPQR/Include/spqr.hpp"
//#include "../../CHOLMOD/Include/cholmod_internal.h"
//#include "../Include/Parallel_LU.hpp"
#include "./paru_mem.h"
paru_symbolic *paru_sym_analyse
( cholmod_sparse *A, cholmod_common *cc) ;

#define Int SuiteSparse_long

// =============================================================================

int main (int argc, char **argv)
{
    // DEBUG(3);
    cholmod_common Common, *cc ;
    cholmod_sparse *A ;
    int mtype ;
    paru_symbolic *LUsym;

    // start CHOLMOD
    cc = &Common ;
    cholmod_l_start (cc) ;

    // A = mread (stdin) ; read in the sparse matrix A
    A = (cholmod_sparse *) cholmod_l_read_matrix (stdin, 1, &mtype, cc) ;
    if (mtype != CHOLMOD_SPARSE)
    {
        printf ("input matrix must be sparse\n") ;
        exit (1) ;
    }
    // paru_sym_analyse(A,cc,LUsym);
    LUsym = paru_sym_analyse (A, cc) ;

    cholmod_l_free_sparse (&A, cc) ;
    paru_freesym(&LUsym,cc);
    ASSERT (LUsym == NULL) ;

    printf ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse);

    cholmod_l_finish (cc) ;
    printf ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse);

}

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

    SuiteSparse_long m, n, anz ; 
    m=LUsym->m=QRsym->m;
    n=LUsym->n=QRsym->n;
    anz=LUsym->anz=QRsym->anz;

    Long nf, *Parent, *Child, *Childp, *Rp, *ColCount, *Super;
    nf = QRsym->nf;
    Parent = QRsym->Parent;
    Child= QRsym->Child;
    Childp= QRsym->Childp;
    Super=QRsym->Super;

    //Staircase structure
    Long *Sp, *Sj, *Sleft, *Qfill, *PLinv;
    Sp=QRsym->Sp;
    Sj=QRsym->Sj;
    Qfill=QRsym->Qfill;
    PLinv=QRsym->PLinv;
    Sleft = QRsym->Sleft;

#ifndef NDEBUG
    /* print fronts*/
    for (SuiteSparse_long f = 0; f < nf; f++) {
        SuiteSparse_long fm, fn, fp;
        fm = QRsym->Fm[f];
        fn = QRsym->Rp[f+1]-QRsym->Rp[f];
        fp = QRsym->Super[f+1]-QRsym->Super[f];

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
    SuiteSparse_long *aParent; //augmented tree size m+nf+1
    SuiteSparse_long *aChild; // size m+nf+1
    SuiteSparse_long *aChildp; // size m+nf+2
    aParent=(SuiteSparse_long *)    paralloc(m+nf+1 ,sizeof(SuiteSparse_long),cc);
    aChild=(SuiteSparse_long *)     paralloc(m+nf+1 ,sizeof(SuiteSparse_long),cc);
    aChildp=(SuiteSparse_long *)    paralloc(m+nf+2 ,sizeof(SuiteSparse_long),cc);

    SuiteSparse_long *rM, *snM; // row map and supernode map
    rM=(SuiteSparse_long *)         paralloc(m      ,sizeof(SuiteSparse_long),cc);
    snM=(SuiteSparse_long *)        paralloc(nf     ,sizeof(SuiteSparse_long),cc);


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

    //brain transplant
    LUsym->nf = QRsym->nf;
    LUsym->Parent = QRsym->Parent; QRsym->Parent=NULL;
    LUsym->Child= QRsym->Child; QRsym->Child=NULL;
    LUsym->Childp= QRsym->Childp; QRsym->Childp=NULL;
    LUsym->Super=QRsym->Super;  QRsym->Super=NULL;

    LUsym->Sp=QRsym->Sp; QRsym->Sp=NULL;
    LUsym->Sj=QRsym->Sj; QRsym->Sj=NULL;
    LUsym->Sleft = QRsym->Sleft; QRsym->Sleft=NULL;

    spqr_freesym (&QRsym, cc) ;


    LUsym->aParent=aParent;
    LUsym->aChildp=aChildp;
    LUsym->aChild=aChild;
    LUsym->row2atree=rM;
    LUsym->super2atree=snM;

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
