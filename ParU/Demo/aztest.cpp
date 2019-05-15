// ===Aznaveh Working on this as a test=========================================
#include "../../SPQR/Include/spqr.hpp"
#include "./paru_mem.h"
int main (int argc, char **argv)
{
    printf("\n=========== Start ======================\n");
    cholmod_common Common, *cc ;
    cholmod_sparse *A ;
    int mtype ;
    spqr_symbolic *QRsym;

    // start CHOLMOD
    cc = &Common ;
    cholmod_l_start (cc) ;

    // A = mread (stdin) ; read in the sparse matrix A
    A = (cholmod_sparse *) cholmod_l_read_matrix (stdin, 1, &mtype, cc) ;
    cc->SPQR_grain = 1;
    cc->useGPU = -1;
    QRsym = spqr_analyze (A, SPQR_ORDERING_CHOLMOD, FALSE,FALSE , FALSE, cc) ;


    if (mtype != CHOLMOD_SPARSE)
    {
        printf ("%%input matrix must be sparse\n") ;
        exit (1) ;
    }

    Long nf, *Parent, *Child, *Childp, *Rp, *ColCount, *Super;
    nf = QRsym->nf;
    Parent = QRsym->Parent;
    Child= QRsym->Child;
    Childp= QRsym->Childp;
    Super=QRsym->Super;

    //Staircase structure
    Long *Sp, *Sj, *Sleft, *Qfill, *Pinv;
    Sp=QRsym->Sp;
    Sj=QRsym->Sj;
    Qfill=QRsym->Qfill;
    Pinv=QRsym->Pinv;
    Sleft = QRsym->Sleft;

    printf ("%ldx%ld nnz:%ld\n",QRsym->m, QRsym->n, QRsym->anz) ; // mxn nnz

    /* print first nonzeros in each column */
    printf("Sleft: ");      //n+2 elements
    for (int j = 0; j < QRsym->n+2; j++) 
        printf("%d ",Sleft[j]); printf("\n");



    /* print Super */
    printf("Super: ");
    for (int j = 0; j < QRsym->nf+1; j++) 
        printf("%d ",Super[j]); printf("\n");


    /* print etree*/
    printf("etree:\t");
    for (int i = 0; i < nf+1 ; i++)
        printf("%d,",Parent[i]); printf("\n");

    /* print Child matrix */
    printf("\n***Child: ");
    for (int f = 0; f < nf+1; f++) {
        printf("%d ",Child[f]);    
    }
    printf("\n$$$Childp: ");
    for (int f = 0; f < nf+2; f++) {
        printf("%d ",Childp[f]);    
    }



    /* print matrix */

    //     for (int j = 0; j < A->ncol; j++) {
    //         for (int p = ((int*)A->p)[j]; p < ((int*)A->p)[j+1]; p++) {
    //             printf("(%d,%d) %lf\n",p+1,j+1, ((double*)A->x)[p]);

    //         }
    //     }



    /* print fronts*/
    //   printf ("nf=%ld \n",QRsym->nf);
    //   for (int f = 0; f < nf; f++) {
    //       Long fm, fn, fp;
    //       fm = QRsym->Fm[f];
    //       fn = QRsym->Rp[f+1]-QRsym->Rp[f];
    //       fp = QRsym->Super[f+1]-QRsym->Super[f];

    //       printf ("Front=%ld #col=%ld #row=%ld #pivotCol=%ld Par=%ld", 
    //               f, fn, fm, fp,Parent[f]);
    //       printf(" #pivot col= %ld",Super[f+1]-Super[f]);
    //       printf ("\nlist of children:\t");
    //       for (int i = Childp[f]; i <= Childp[f+1]-1; i++) 
    //           printf ("%ld ",Child[i]);
    //       printf ("\n\n");
    //   }

    /*Constructing augmented tree */
    int m=QRsym->m;
    int *aParent; //augmented tree size m+nf+1
    int *aChild; // size m+nf+1
    int *aChildp; // size m+nf+2
    aParent=(int *)paralloc(m+nf+1,sizeof(int));
    aChild=(int *)paralloc(m+nf+1,sizeof(int));
    aChildp=(int *)paralloc(m+nf+2,sizeof(int));
    int *rM, *snM; // row map and supernode map
    rM=(int *)paralloc(m,sizeof(int));
    snM=(int *)paralloc(nf,sizeof(int));

    for (int f = 0; f < nf; f++) snM[f]=-1;
    for (int i = 0; i < m; i++) rM[i]=-1;
    for (int i = 0; i < m+nf+1; i++) aParent[i]=-1;
    for (int i = 0; i < m+nf+1; i++) aChild[i]=-1;
    for (int i = 0; i < m+nf+2; i++) aChildp[i]=-1;


    aChildp[0]=0;
    int offset=0; //number of rows visited in each iteration orig front+ rows
    int lastChildFlag=0;
    int childpointer=0;
    for (int f = 0; f < nf; f++) {

        int numRow =Sleft[Super[f]+1]-Sleft[Super[f]] ;
 //       printf("~%d",numRow);
 //       printf("\n#%d\n",offset);

        int numoforiginalChild=0;
        if (lastChildFlag){  // the current node is the parent
//            printf("Childs of %d: ",f);
            numoforiginalChild=Child[Childp[f+1]-1]-Child[Childp[f]]+1;
            for (int i = Child[Childp[f]]; i < Child[Childp[f+1]]; i++){
 //               printf("%d,", i);
                aParent[ snM[i]]=offset+numRow;
                aChild[childpointer++]=snM[i];
            }
        }


        for(int i=offset ; i < offset+numRow ; i++)
            aChildp[i+1]=aChildp[i];


        for (int i = Sleft[Super[f]]; i < Sleft[Super[f]+1]; i++){ // number of rows
            rM[i]=i+f;
            aParent[i+f]=offset+numRow;
            aChild[childpointer++]=i+f;
        }

        offset+=numRow;
        snM[f]=offset++;
        aChildp[offset] = aChildp[offset-1]+numRow+numoforiginalChild;

        if( Parent[f] == f+1)  //last child
            lastChildFlag=1;  
        else
            lastChildFlag=0;  

    }
    printf("\nsuper node->aP ");    for (int f = 0; f < nf; f++) printf("%d ",snM[f]); printf("\n");
    printf("row->aP ");    for (int i = 0; i < m; i++) printf("%d ",rM[i]); printf("\n");
    printf("aP: ");    for (int i = 0; i < m+nf+1; i++) printf("%d ",aParent[i]); printf("\n");
    printf("aChildp: ");    for (int i = 0; i < m+nf+1; i++) printf("%d ",aChildp[i]); printf("\n");
    printf("aChild: ");    for (int i = 0; i < m+nf+1; i++) printf("%d ",aChild[i]); printf("\n");

    for(int i=0; i< m+nf; i++){
        printf("anode:%d",i);
        for(int c=aChildp[i]; c< aChildp[i+1]; c++)
            printf(" %d,",aChild[c]);
        printf("\n");
    }




    parfree(aParent);
    parfree(aChild);
    parfree(aChildp);
    parfree(rM);
    parfree(snM);


    cholmod_l_free_sparse (&A, cc) ;
    cholmod_l_finish (cc) ;
    printf("\n~~~~~~~~The End~~~~~~~~~~~~~~\n");

    return (0) ;
}
