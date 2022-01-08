////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_tasked_trsm  //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief      a wrapper around  BLAS_TRSM for tasking
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define L 4096
void paru_tasked_trsm(Int f, int *m, int *n, double *alpha, double *a, int *lda,
                      double *b, int *ldb)
{
    DEBUGLEVEL(0);
    if (*n < L)
    //if(1)
    {
        PRLEVEL(1, ("%% No tasking for TRSM (%dx%d) in %ld\n", *m, *n, f));
        //BLAS_DTRSM ("L", "L", "N", "U", m, n, alpha, a, lda, b, ldb);
        cblas_dtrsm (CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, 
                CblasUnit,
                *m, *n, *alpha, a, *lda, b, *ldb);
    }
    else
    {
        PRLEVEL(1, ("%% YES tasking for TRSM (%dx%d) in %ld", *m, *n, f));
        Int num_blocks = *n / L + 1;
        Int len_bloc = *n / num_blocks;
        PRLEVEL(1, ("%%  num_blocks = %ld", num_blocks));
        #pragma omp parallel proc_bind(close)
        #pragma omp single nowait
        {
            for (Int J = 0; J < num_blocks; J++)
            {
                BLAS_INT n_b = (J + 1) == num_blocks  ? 
                    (*n - J * len_bloc) : len_bloc;
                PRLEVEL(1, ("%%  n_b= %d\n", n_b));
               #pragma omp task 
                //BLAS_DTRSM("L", "L", "N", "U", m, &n_b, alpha, a, lda, 
                //        (b + J * len_bloc* *ldb), ldb);
                cblas_dtrsm (CblasColMajor, CblasLeft, CblasLower, 
                        CblasNoTrans, CblasUnit,
                        *m, n_b, *alpha, a, *lda, 
                        (b + J * len_bloc* *ldb), *ldb);
 
            }
        }  
    }
}
