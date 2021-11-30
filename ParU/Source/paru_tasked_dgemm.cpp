////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_tasked_dgemm //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief      a wrapper around  BLAS_DGEMM for tasked base dgemmed 
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define L 4096
void paru_tasked_dgemm(Int f, char *transa, char *transb, BLAS_INT *M,
                       BLAS_INT *N, BLAS_INT *K, double *alpha, double *A,
                       BLAS_INT *lda, double *B, BLAS_INT *ldb, double *beta,
                       double *C, BLAS_INT *ldc)
{
    DEBUGLEVEL(0);
    if (*M < L && *N < L)
    { //TODO use a nested loop for very small dgemms?
        PRLEVEL(1, ("%% No tasking for DGEMM (%dx%d) in %ld\n", *M, *N, f));
        BLAS_DGEMM(transa, transb, M, N, K, alpha, A, lda, B, ldb, beta, C,
                   ldc);
    }
    else
    {
        PRLEVEL(1, ("%% YES tasking for DGEMM (%dx%d) in %ld", *M, *N, f));

        Int num_col_blocks =  *N / L + 1 ;
        Int num_row_blocks =  *M / L + 1 ;

        Int len_col = *N / num_col_blocks;
        Int len_row = *M / num_row_blocks;
 
        PRLEVEL(1, ("%% col-blocks=%ld,row-blocks=%ld) \n", num_col_blocks,
                    num_row_blocks));
 
        #pragma omp parallel proc_bind(close)
        #pragma omp single
        {
            for (Int I = 0; I < num_row_blocks; I++)
            {
                //BLAS_INT m = (I + 1) * L > *M ? (*M - I * L) : L;
                BLAS_INT m = (I + 1) == num_row_blocks ? 
                    (*M - I * len_row) : len_row;

                for (Int J = 0; J < num_col_blocks; J++)
                {
                    //BLAS_INT n = (J + 1) * L > *N ? (*N - J * L) : L;
                    BLAS_INT n = (J + 1) == num_col_blocks ? 
                        (*N - J * len_col) : len_col;
                    PRLEVEL(1, ("%% I=%ld J=%ld m=%d n=%d in %ld\n", I, J,
                                m, n, f));
                    #pragma omp task
                    BLAS_DGEMM(transa, transb, &m, &n, K, alpha,
                            A + (I * len_row), lda, B + (J * len_col * *ldb), ldb,
                            beta, C + (J * *ldc * len_col + I * len_row), ldc);
                }
            }
        }  
    }
}
