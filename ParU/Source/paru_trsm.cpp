/** =========================================================================  /
 * =======================  paru_trsm =======================================  /
 * ========================================================================== */
/*! @brief trsm wrapper
 *     l11*u12=a12 and u12 is unkonwn  so we need this:
 *              op( A ) * X = alpha*B
 *             part(pF) * X = 1 * upart
 *             part(pF) * upart = 1 * upart
 *      SIDE = 'L' or 'l'
 *      UPLO = 'L' or 'l'; A is a lower triangular matrix.
 *      TRANSA = 'N' or 'n'   op( A ) = A.
 *      DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *      M     M specifies the number of rows of B.    (fp)
 *      N     N specifies the number of columns of B. (colCount)
 *      ALPHA, (alpha = 1.0)
 *      A (pF)
 *      LDA  leading dimension of A. (rowCount)
 *      B (upart)
 *      LDB  leading dimension of B.  (fp)
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

/// Already defined in /SuiteSparse/CHOLMOD/Include/cholmod_blas.h
//          use BLAS_DTRSM insted
// extern "C"  int dtrsm_(char *side, char *uplo, char *transa, char *diag,
//                           int *m, int *n, double *alpha, double *a, int *lda,
//                                                         double *b, int *ldb);

Int paru_trsm(double *pF, double *uPart, Int fp, Int rowCount, Int colCount)
{
    DEBUGLEVEL(0);
    BLAS_INT mB = (BLAS_INT)fp;
    BLAS_INT nB = (BLAS_INT)colCount;
    double alpha = 1.0;
    BLAS_INT lda = (BLAS_INT)rowCount;
    BLAS_INT ldb = (BLAS_INT)fp;

#ifndef NDEBUG  // Printing the  U part
    Int p = 1;
    PRLEVEL(p, ("mB=%d nB = %d alpha = %f \n", mB, nB, alpha));
    PRLEVEL(p, ("lda =%d ldb =%d\n", lda, ldb));
    PRLEVEL(p, ("(I)U Before Trsm: %ld x %ld\n", fp, colCount));
    for (Int i = 0; i < fp; i++)
    {
        for (Int j = 0; j < colCount; j++)
            PRLEVEL(p, (" %2.5lf\t", uPart[j * fp + i]));
        PRLEVEL(p, ("\n"));
    }
#endif

    BLAS_DTRSM("L", "L", "N", "U", &mB, &nB, &alpha, pF, &lda, uPart, &ldb);

#ifndef NDEBUG  // Printing the  U part
    PRLEVEL(p, ("(I)U After Trsm: %ld x %ld\n", fp, colCount));
    for (Int i = 0; i < fp; i++)
    {
        for (Int j = 0; j < colCount; j++)
            PRLEVEL(p, (" %2.5lf\t", uPart[j * fp + i]));
        PRLEVEL(p, ("\n"));
    }
#endif

    return 0;
}
