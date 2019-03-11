#include "Parallel_LU.hpp"

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


    double F[9]= {1, 2, 3, 1, 8, 3, 4, 5, 6};
    Int fsRowList[3] = {5, 12, 7};
    Int rowCount = 3; Int fp = 3;


//    double F[4]= {1, 2, 3, 1};
//    Int fsRowList[2] = {5, 12};
//    Int rowCount = 2; Int fp = 2;
 
//    paru_panel_factorize ( F, fsRowList, rowCount, fp, 
//            fp, 0, rowCount, paruMatInfo);

    paru_factorize(F, fsRowList, rowCount, fp, paruMatInfo);

    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*\n");

    cholmod_l_free_sparse (&A, cc);
//    paru_freemat (&paruMatInfo, cc);

    paru_freesym (&LUsym,cc);
    cholmod_l_finish (cc);

    for (Int j=0; j < fp; j++){
        printf ("%ld\t", fsRowList[j]);
        for (Int i=0; i < rowCount; i++){
            printf ("%2.4lf\t", F[i*rowCount+j]);
        }
        printf("\n");
    }

}
