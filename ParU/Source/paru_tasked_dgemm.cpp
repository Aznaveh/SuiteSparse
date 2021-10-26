////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_tasked_dgemm //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief      a wrapper around  BLAS_DGEMM
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
void paru_tasked_dgemm(char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
                       BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda,
                       double *B, BLAS_INT *ldb, double *beta, double *C,
                       BLAS_INT *ldc)
{
    BLAS_DGEMM(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

