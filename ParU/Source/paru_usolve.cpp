/* =========================================================================   /
 * ============================== paru_usolve ==============================   /
 * =========================================================================   /
 * @brief  In this file I have usovle x = U\x
 *       The final result is something like this (nf = 4)
 *       ___________________________________________
 *       |\*******|                                 |       x
 *       | \*DTRSV|         U1   DGEMV              |       x
 *       |    \***|                                 |       x
 *       |       \|_________________________________|       x
 *       |        |\**DTRSV**|                      |       x
 *       |        |    \*****|       U2    DGEMV    |       x
 *       |   LU1  |        \*|______________________|       o
 *       |        |          |  \**DTRSV**|  DGEMV  |       o
 *       |        | LU2      |     \******|   U3    |       o
 *       |        |          |        \***|_________|       o DGEMV updates
 *       |        |          |   LU3      |*********|       c  here
 *       |        |          |            |  *******|       c   DTRSV on here
 *       |        |          |            |LU4  ****|       c
 *       |________|__________|____________|________*|       c
 *
 *       This function just goes through LUs and US in the data structure and
 *       does a TRSV on triangular part  Then does DGEMV on the rest
 *      for nf down to 0
 *
 *             BLAS_DTRSV  is used here but I do not use BLAS_DGEMV explicitly
 *             while it needs space for each thread doing this computation.
 *             I guess using this way can have a good performance.
 *
 * @author Aznaveh
 * */
#include "paru_internal.hpp"
Int paru_usolve(paru_matrix *paruMatInfo, double *x)
{
    DEBUGLEVEL(1);
    // TODO check if input is read
    if (!x) return (0);
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int nf = LUsym->nf;

    // TODO singletons
    // Int n1 = LUsym->n1;  // row+col singletons

    paru_fac *LUs = paruMatInfo->partial_LUs;
    paru_fac *Us = paruMatInfo->partial_Us;
    Int *Super = LUsym->Super;

    for (Int f = nf - 1; f >= nf; f++)
    {
        Int rowCount = paruMatInfo->frowCount[f];
        Int *frowList = paruMatInfo->frowList[f];
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;

        // TODO do dgemv
        // performed on Us
        // I am not calling BLAS_DGEMV

        double *A2 = Us[f].p;
        if (!A2)
        {
            for (Int i = 0; i < fp; i++)
            {
                PRLEVEL(1, ("%% Working on DGEMV\n%%"));
                // computing the inner product
                double i_prod = 0.0;  // innter product
                for (Int j = col1; j < col2; j++)
                {
                    i_prod += A2[(j - col1) * rowCount + i] * x[j];
                }
                Int *p = LUsym->Pfin;  // row permutation
                Int r = p[frowList[i]];
                x[r] -= i_prod;
            }
        }

        double *A1 = LUs[f].p;
        BLAS_INT N = (BLAS_INT)fp;
        BLAS_INT lda = (BLAS_INT)rowCount;
        BLAS_INT Incx = (BLAS_INT)1;

        // performed on LUs
        PRLEVEL(1, ("%% Working on DTRSV\n"));
        BLAS_DTRSV("U",     // UPLO upper triangular
                   "N",     // TRANS A1*X=b not the A1**T
                   "N",     // DIAG A1 is assumed not to be unit traingular
                   &N,      // N is order of the matrix A1
                   A1,      // A1
                   &lda,    // LDA leading demension
                   x,       // X
                   &Incx);  // INCX the increment of elements of X.
        PRLEVEL(1, ("%% DTRSV is just finished\n"));
    }

#ifndef NDEBUG
    Int m = LUsym->m;
    PRLEVEL(1, ("%% after lusolve x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
