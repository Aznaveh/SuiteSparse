#include "Parallel_LU.hpp"

// =============================================================================
int main (int argc, char **argv)
{
    DEBUGLEVEL(0); 
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
    paru_matrix *paruMatInfo = paru_init_rowFronts (A, LUsym, cc);


    cholmod_l_free_sparse (&A, cc);
    paru_freemat (&paruMatInfo, cc);

    //paru_freesym (&LUsym,cc);

 //  spqr_symbolic *QRsym = 
 //      spqr_analyze (A, SPQR_ORDERING_CHOLMOD, FALSE,FALSE , FALSE, cc);
 //           spqr_freesym (&QRsym, cc);
  //  PRLEVEL (1, ("malloc_count %ld inuse %ld\n", 
   //cc->malloc_count, cc->memory_inuse));

    cholmod_l_finish (cc);
    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    //   printf 
    //   ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse);
}
