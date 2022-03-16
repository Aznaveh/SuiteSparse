////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_dgemm /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief      A wraper around dgemm for outer product.
 *
 *          This does the outer product
 *              C := alpha*op( A )*op( B ) + beta*C,
 *              el:= -1*part(pF) * uPart + zero
 *              TRANSA = 'N'   op(A) = A
 *              TRANSB = 'N'   op(B) = B
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

Int paru_dgemm(Int f, double *pF, double *uPart, double *el, Int fp,
               Int rowCount, Int colCount, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% rowCount =%ld  ", rowCount));
    PRLEVEL(1, ("%% colCount =%ld  ", colCount));
    PRLEVEL(1, ("%% fp =%ld\n", fp));

    BLAS_INT mA = (BLAS_INT)(rowCount - fp);
    BLAS_INT nB = (BLAS_INT)colCount;
    BLAS_INT nA = (BLAS_INT)fp;

    PRLEVEL(1, ("%% mA =%d  ", mA));
    PRLEVEL(1, ("%% nB =%d  ", nB));
    PRLEVEL(1, ("%% nA =%d\n", nA));

#ifndef NDEBUG
    double *Ap = pF + fp;
    PRLEVEL(1, ("%% A =\n"));
    for (Int i = 0; i < mA; i++)
    {
        PRLEVEL(1, ("%% "));
        for (Int j = 0; j < nA; j++)
            PRLEVEL(1, ("%2.4lf\t", Ap[j * rowCount + i]));
        PRLEVEL(1, ("\n"));
    }

    Int mB = nA;
    double *Bp = uPart;
    PRLEVEL(1, ("%% B =\n"));
    for (Int i = 0; i < mB; i++)
    {
        PRLEVEL(1, ("%% "));
        for (Int j = 0; j < nB; j++) PRLEVEL(1, ("%2.4lf\t", Bp[j * mB + i]));
        PRLEVEL(1, ("\n"));
    }

    double *Cp = el;
    PRLEVEL(1, ("%%Before DGEMM C =\n"));
    for (Int i = 0; i < mA; i++)
    {
        PRLEVEL(1, ("%% "));
        for (Int j = 0; j < nB; j++)
        {
            PRLEVEL(1, ("%2.4lf\t", Cp[j * mA + i]));
        }
        PRLEVEL(1, ("\n"));
    }
#endif

    // alpha = -1;
    BLAS_INT lda = (BLAS_INT)rowCount;
    BLAS_INT ldb = (BLAS_INT)fp;
    BLAS_INT ldc = (BLAS_INT)(rowCount - fp);

    PRLEVEL(1, ("%% lda =%d  ", lda));
    PRLEVEL(1, ("%% ldb =%d  ", ldb));
    PRLEVEL(1, ("%% ldc =%d\n", ldc));

    // double beta = 0;  // U part is not initialized

    paru_tasked_dgemm(f, mA, nB, nA, pF + fp, lda, uPart, ldb, 0, el, ldc, Num);

#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("%%After DGEMM C =\n"));
    for (Int i = 0; i < mA; i++)
    {
        PRLEVEL(PR, ("%% "));
        for (Int j = 0; j < nB; j++) PRLEVEL(PR, ("%2.4lf\t", Cp[j * mA + i]));
        PRLEVEL(PR, ("\n"));
    }
#endif

    return 0;
}
