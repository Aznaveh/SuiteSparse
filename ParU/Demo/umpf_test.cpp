/** =========================================================================  /
 * =======================  paru_test========================================  /
 * ==========================================================================  /
 * @brief    test to see how to call umfpack symbolic analysis
 * 
 * @author Aznaveh
 * */
#include "Parallel_LU.hpp"

int main (int argc, char **argv)
{
    DEBUGLEVEL(0); 

    cholmod_common Common, *cc;
    cholmod_sparse *A;
    int mtype;
    paru_symbolic *LUsym;


    //~~~~~~~~~Reading the input matrix and test if the format is OK~~~~~~~~~~~~
    // start CHOLMOD
    cc = &Common;
    cholmod_l_start (cc);

    // A = mread (stdin) ; read in the sparse matrix A
    A = (cholmod_sparse *) cholmod_l_read_matrix (stdin, 1, &mtype, cc);
    if (A == NULL){
        printf ("Paru: input matrix is invalid\n");
        exit (1);
    }

    if (mtype != CHOLMOD_SPARSE){
        printf ("Paru: input matrix must be sparse\n");
        exit (1);
    }

    if (A->xtype != CHOLMOD_REAL){
        printf ("Paru: input matrix must be real\n");
        exit (1);
    }

    //~~~~~~~~~~~~~~~~~~~Starting computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LUsym = paru_analyze (A, cc);
    if (LUsym == NULL) {
        exit(0);
    }



    int scale = 0;
    if (argc == 3){
        scale = atoi(argv[2]);
        if (scale) 
            PRLEVEL (1, ("The input matrix will be scaled\n"));
    }
        
//    paru_matrix *paruMatInfo = paru_init_rowFronts (A, scale, LUsym, cc);
//    if (paruMatInfo == NULL) {
//        exit(0);
//    }


//    Int m,n,nf;
//    m = paruMatInfo-> m;
//    n = paruMatInfo-> n;
//    nf = paruMatInfo->LUsym->nf;
//

//    for (Int i = 0; i < nf; i++) {
//        //paru_assemble (paruMatInfo, i, cc);
//        paru_front (paruMatInfo, i, cc);
//    }
    //~~~~~~~~~~~~~~~~~~~End computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    //~~~~~~~~~~~~~~~~~~~Free Everything~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cholmod_l_free_sparse (&A, cc);

//    paru_freemat (&paruMatInfo, cc);
    paru_freesym (&LUsym,cc);

    cholmod_l_finish (cc);

}
