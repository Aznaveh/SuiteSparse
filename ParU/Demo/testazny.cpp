#include "Parallel_LU.hpp"
#define PRINTCBsTUPLES 0

// =============================================================================
void paru_print_element (paru_matrix *paruMatInfo, Int e){
    DEBUGLEVEL(0);
    // print out contribution blocks
    Element **elementList; 
    elementList = paruMatInfo->elementList;
    Element *curEl = elementList[e];

    Int morign = paruMatInfo->m;
    Int nf = paruMatInfo->LUsym->nf;

    if ( e > morign + nf +1){
        printf("Element %ld is out of range; just %ld elements \n", 
                e,  morign + nf +1);
        return;
    }

    if (curEl == NULL){
        printf("Element %ld is empty\n",e );
        return;
    }

    Int m,n;
    m = curEl->nrows;
    n = curEl->ncols;
   
    Int *el_colrowIndex = (Int*)(curEl+1);     // pointers to element index 
    double *el_colrowNum = (double*)(el_colrowIndex + m + n); //and values

    PRLEVEL (1, ("el_colrowIndex =%p, el_colrowNum = %p \n", 
                el_colrowIndex, el_colrowNum));

    printf("\n"); 
    printf("Element %ld is %ld x %ld:\n", e, m, n);


    printf("\t"); 
    for (int j = 0; j < n; j++) 
        printf("%ld\t", el_colrowIndex [j] );
    printf("\n"); 
    for (int i = 0; i < m; i++) {
        printf("%ld\t",el_colrowIndex [n+i] );
        for (int j = 0; j < n; j++) {
            double value =  el_colrowNum [i*m + j];
            printf("%2.4lf\t",value );
        }
        printf("\n"); 
    }

}

void paru_print_tupleList (tupleList *listSet, Int index){
    DEBUGLEVEL(0);
    PRLEVEL (1, ("listSet =%p\n", listSet));

    if (listSet == NULL) {
       printf("Empty tuple\n"); 
       return;
    }

    tupleList cur= listSet [index];
    Int numTuple = cur.numTuple;
    Tuple *l = cur.list;

    printf(" There are %ld tuples in this list:\n", numTuple);
    for (Int i = 0; i < numTuple; i++) {
       Tuple curTpl = l [i];
        printf(" (%ld,%ld)", curTpl.e, curTpl.f);
    }
    printf("\n"); 
}
int main (int argc, char **argv)
{
    DEBUGLEVEL(0); 
    printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    cholmod_common Common, *cc;
    cholmod_sparse *A;
    int mtype;
    //paru_symbolic *LUsym;

    // start CHOLMOD
    cc = &Common;
    cholmod_l_start (cc);

    // A = mread (stdin) ; read in the sparse matrix A
    A = (cholmod_sparse *) cholmod_l_read_matrix (stdin, 1, &mtype, cc);
    if (mtype != CHOLMOD_SPARSE)
    {
        printf ("input matrix must be sparse\n");
        exit (1);
    }

//    LUsym = paru_sym_analyse (A, cc);
    paru_symbolic LUsym (A,cc);

    if (LUsym == NULL) {
        exit(0);
    }


    paru_matrix *paruMatInfo = paru_init_rowFronts (A, LUsym, cc);
    if (paruMatInfo == NULL) {
        exit(0);
    }

    Int m,n,nf;
    m = paruMatInfo-> m;
    n = paruMatInfo-> n;
    nf = paruMatInfo->LUsym->nf;
 
   for (Int i = 0; i < nf; i++) {
    paru_assemble (paruMatInfo, i, cc);
   }
    
    if (PRINTCBsTUPLES){
        for (int i = 0; i < m+nf+1; ++i) 
            paru_print_element (paruMatInfo, i);

        tupleList *RowList = paruMatInfo -> RowList;
        tupleList *ColList = paruMatInfo -> ColList;  

        printf ("RowList:\n");
        for (int i = 0; i < m; ++i){ 
            printf("row %ld :",i);
            paru_print_tupleList (RowList , i);
        }

        printf ("ColList:\n");
        for (int i = 0; i < n; ++i) {
            printf("col %ld :",i);
            paru_print_tupleList (ColList , i);
        }

    }


    cholmod_l_free_sparse (&A, cc);
    paru_freemat (&paruMatInfo, cc);

    paru_freesym (&LUsym,cc);


    //  spqr_symbolic *QRsym = 
    //      spqr_analyze (A, SPQR_ORDERING_CHOLMOD, FALSE,FALSE , FALSE, cc);
    //           spqr_freesym (&QRsym, cc);
    //  PRLEVEL (1, ("malloc_count %ld inuse %ld\n", 
    //cc->malloc_count, cc->memory_inuse));

    cholmod_l_finish (cc);
    printf("***************************************************************\n");
    //   printf 
    //   ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse);
}
