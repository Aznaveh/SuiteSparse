////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_tasked_dgemm //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief      a wrapper around  BLAS_DGEMM
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define L 128
void paru_tasked_dgemm(Int f, char *transa, char *transb, BLAS_INT *M,
                       BLAS_INT *N, BLAS_INT *K, double *alpha, double *A,
                       BLAS_INT *lda, double *B, BLAS_INT *ldb, double *beta,
                       double *C, BLAS_INT *ldc)
{
    DEBUGLEVEL(1);
    if (*M < L && *N < L)
    {
        PRLEVEL(1, ("%% No tasking for BLAS (%dx%d) in %ld\n", *M, *N, f));
        BLAS_DGEMM(transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C,
                   ldc);
    }
    else
    {
        PRLEVEL(1, ("%% YES tasking for BLAS (%dx%d) in %ld", *M, *N, f));
        Int num_col_blocks = (*M % L == 0) ? *M / L : *M / L + 1;
        Int num_row_blocks = (*N % L == 0) ? *N / L : *N / L + 1;
        PRLEVEL(1, ("%% col-blocks=%ld,row-blocks=%ld) \n", 
                    num_col_blocks, num_row_blocks));
        #pragma omp parallel
        {
            #pragma omp single
            {
                for (Int I = 0; I < num_col_blocks; I++)
                {
                    BLAS_INT m = (I + 1) * L > *M ? (*M - I * L) : L;

                    for (Int J = 0; J < num_row_blocks; J++)
                    {
                        BLAS_INT n = (J + 1) * L > *N ? (*N - J * L) : L;
                        PRLEVEL(1, ("%% I=%ld J=%ld m=%ld n=%ld in %ld\n", I, J,
                                    m, n, f));
                        #pragma omp task
                        BLAS_DGEMM(transa, transb, &m, &n, K, alpha, A + (I * L),
                                   lda, B + (J * *K * L), ldb, beta,
                                   C + (J * *M * L + I * L), ldc);
                    }
                }

            }  // end of single region
        }      // end of parallel region
    }
}
