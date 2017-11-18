#include "Parallel_LU.hpp"
#define PRINTCBsTUPLES 0

int main (int argc, char **argv)
{
    DEBUGLEVEL(0); 
    printf("%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
    cholmod_common Common, *cc;
    cholmod_sparse *A;
    int mtype;
    paru_symbolic *LUsym;

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

    LUsym = paru_sym_analyse (A, cc);
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

    if (PRINTCBsTUPLES){
        for (int i = 0; i < m+nf+1; ++i) 
            paru_print_element (paruMatInfo, i);

        tupleList *RowList = paruMatInfo -> RowList;
        tupleList *ColList = paruMatInfo -> ColList;  

        printf ("RowList:\n");
        for (Int i = 0; i < m; ++i){ 
            printf("row %ld :",i);
            paru_print_tupleList (RowList , i);
        }

        PRLEVEL (1, ("ColList =%p\n", ColList));
        printf ("ColList:\n");
        for (Int i = 0; i < n; ++i) {
            printf("col %ld :",i);
            paru_print_tupleList (ColList , i);
        }

    }

    for (Int i = 0; i < nf; i++) {
        paru_assemble (paruMatInfo, i, cc);
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
    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*\n");
    //   printf 
    //   ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse);
}
