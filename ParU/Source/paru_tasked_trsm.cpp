////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_tasked_trsm  //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief      a wrapper around  BLAS_TRSM for tasking
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define L 256
void paru_tasked_trsm(Int f, char *side, char *uplo, char *transa, char *diag,
                      int *m, int *n, double *alpha, double *a, int *lda,
                      double *b, int *ldb)
{
    DEBUGLEVEL(0);
    if (*n < L)
    {
        PRLEVEL(1, ("%% No tasking for TRSM (%dx%d) in %ld\n", *m, *n, f));
        BLAS_DTRSM (side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
    else
    {
        PRLEVEL(1, ("%% YES tasking for BLAS (%dx%d) in %ld", *m, *n, f));
        Int num_blocks = (*n %L == 0) ? *n / L : *n / L + 1;
        PRLEVEL(1, ("%%  num_blocks = %ld", num_blocks));
        #pragma omp parallel
        {
            #pragma omp single
            {
                for (Int J = 0; J < num_blocks; J++)
                {
                    BLAS_INT n_b = (J + 1) * L > *n ? (*n - J * L) : L;
                    PRLEVEL(1, ("%%  n_b= %d\n", n_b));
                    #pragma omp task
                    BLAS_DTRSM(side, uplo, transa, diag, m, &n_b, alpha, a, lda, 
                            (b + J * L* *ldb), ldb);
                }
            }  // end of single region
        }      // end of parallel region
    }
}
