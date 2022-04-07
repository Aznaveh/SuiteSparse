////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_trsm //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*! @brief trsm wrapper
 *
 *       l11*u12=a12 and u12 is unkonwn  so we need this:
 *                op( A ) * X = alpha*B
 *               part(pF) * X = 1 * upart
 *               part(pF) * upart = 1 * upart
 *        SIDE = 'L' or 'l'
 *        UPLO = 'L' or 'l'; A is a lower triangular matrix.
 *        TRANSA = 'N' or 'n'   op( A ) = A.
 *        DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *        M     M specifies the number of rows of B.    (fp)
 *        N     N specifies the number of columns of B. (colCount)
 *        ALPHA, (alpha = 1.0)
 *        A (pF)
 *        LDA  leading dimension of A. (rowCount)
 *        B (upart)
 *        LDB  leading dimension of B.  (fp)
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

Int paru_trsm(Int f, double *pF, double *uPart, Int fp, Int rowCount,
              Int colCount, paru_work *Work, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);
    BLAS_INT mB = (BLAS_INT)fp;
    BLAS_INT nB = (BLAS_INT)colCount;
    double alpha = 1.0;
    BLAS_INT lda = (BLAS_INT)rowCount;
    BLAS_INT ldb = (BLAS_INT)fp;

#ifndef NDEBUG  // Printing the  U part
    PRLEVEL(1, ("TRSM (%dx%d) (%dx%d) \n", mB, mB, mB, nB));
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

    paru_tasked_trsm(f, mB, nB, alpha, pF, lda, uPart, ldb, Work, Num);

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
