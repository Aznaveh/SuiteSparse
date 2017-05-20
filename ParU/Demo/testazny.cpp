#include "Parallel_LU.hpp"

// =============================================================================
int main (int argc, char **argv)
{
    DEBUGLEVEL(0); 
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


    PRLEVEL (1, ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse));


    cholmod_l_finish (cc) ;
 //   printf ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse);
//Valgrind :)
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
}

