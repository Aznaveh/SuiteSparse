/* =========================================================================   /
 * ============================== paru_lsolve ==============================   /
 * =========================================================================   /
 * @brief  In this file I have lsovle x = L\x
 *       The final result is something like this (nf = 4)
 *       ___________________________________________
 *       |\       |                                 |       c
 *       |*\      |         U1                      |       c
 *       |****\   |                                 |       c  DTRSV on here
 *       |*DTRSV*\|_________________________________|       c
 *       |******* |\         |                      |       c
 *       |        |****\     |       U2             |       x
 *       |   LU1  |*DTRSV**\ |______________________|       x DGEMV updates
 *       |        |          |**\         |         |       x   here
 *       |        | LU2      |*****\      |   U3    |       x
 *       | DGEMV  |          |*DTRSV**\   |_________|       x
 *       |        |          |   LU3      |* LU4    |       x
 *       |        | DGEMV    |  DGEMV     |****     |       x
 *       |        |          |            |DTRSV*   |       x
 *       |________|__________|____________|_________|       x
 *
 *       This function just goes through LUs in the data structure and does a
 *       TRSV on triangular part
 *       Then does DGEMV on the rest for 0 to nf
 *
 *             BLAS_DTRSV  is used here but I do not use BLAS_DGEMV explicitly
 *             while it needs space for each thread doing this computation.
 *             I guess using this way can have a good performance.
 *
 * @author Aznaveh
 * */
#include "paru_internal.hpp"
Int paru_lsolve(paru_matrix *paruMatInfo, double *x)
{
    DEBUGLEVEL(0);
    // TODO check if input is read
    if (!x) return (0);
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int nf = LUsym->nf;

#ifndef NDEBUG
    Int PR = 1;
#endif
    Int n1 = LUsym->n1;   // row+col singletons
    Int *Ps = LUsym->Ps;  // row permutation S->LU

    // singletons
    Int rs1 = LUsym->rs1;
    if (rs1 > 0)
    {
        Int cs1 = LUsym->cs1;

        for (Int j = cs1; j < n1; j++)
        {
            PRLEVEL(PR, ("j = %ld\n", j));
            Int *Slp = LUsym->lstons.Slp;
            Int *Sli = LUsym->lstons.Sli;
            double *Slx = LUsym->lstons.Slx;
            ASSERT(Sli != NULL && Slx != NULL && Slp != NULL);
            Int diag = Slp[j - cs1];
            PRLEVEL(PR, (" x[%ld]=%.2lf Slx[%ld]=%.2lf\n", j, x[j], diag,
                         Slx[diag]));
            x[j] /= Slx[diag];
            PRLEVEL(PR, (" After x[%ld]=%.2lf \n", j, x[j]));

            for (Int p = Slp[j - cs1] + 1; p < Slp[j - cs1 + 1]; p++)
            {
                Int r = Sli[p] < n1 ? Sli[p] : Ps[Sli[p] - n1] + n1;
                PRLEVEL(PR, (" r=%ld\n", r));
                x[r] -= Slx[p] * x[j];
                PRLEVEL(PR, ("A x[%ld]=%.2lf\n", Sli[p], x[Sli[p]]));
            }
        }
    }
    PRLEVEL(PR, ("%%lsove singletons finished.\n"));

    paru_fac *LUs = paruMatInfo->partial_LUs;
    Int *Super = LUsym->Super;

    for (Int f = 0; f < nf; f++)
    {
        Int rowCount = paruMatInfo->frowCount[f];
        Int *frowList = paruMatInfo->frowList[f];
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;
        double *A = LUs[f].p;
        double *X = x + col1 + n1;

        BLAS_INT N = (BLAS_INT)fp;
        BLAS_INT lda = (BLAS_INT)rowCount;
        BLAS_INT Incx = (BLAS_INT)1;

        PRLEVEL(1, ("%% Working on DTRSV\n"));
        BLAS_DTRSV("L",     // UPLO lower triangular
                   "N",     // TRANS A*X=b not the A**T
                   "U",     // DIAG A is assumed to be unit traingular
                   &N,      // N is order of the matrix A
                   A,       // A
                   &lda,    // LDA leading demension
                   X,       // X
                   &Incx);  // INCX the increment of elements of X.
        PRLEVEL(1, ("%% DTRSV is just finished\n"));
#ifndef NDEBUG
        PR = 2;
        PRLEVEL(PR, ("%% LUs:\n%%"));
        for (Int r = 0; r < rowCount; r++)
        {
            PRLEVEL(PR, ("%% %ld\t", frowList[r]));
            for (Int c = col1; c < col2; c++)
                PRLEVEL(PR, (" %2.5lf\t", A[(c - col1) * rowCount + r]));
            PRLEVEL(PR, ("\n"));
        }

        PRLEVEL(2, ("%% lda = %d\n%%", lda));
        PRLEVEL(2, ("%% during lsolve x [%ld-%ld)is:\n%%", col1, col2));
        // for (Int k = col1; k < col2; k++)
        Int m = LUsym->m;
        for (Int k = 0; k < m; k++)
        {
            PRLEVEL(1, (" %.2lf, ", x[k]));
        }
        PRLEVEL(1, (" \n"));
#endif

        // TODO do dgemv
        // I am not calling BLAS_DGEMV

        for (Int i = fp; i < rowCount; i++)
        {
            PRLEVEL(1, ("%% Working on DGEMV\n%%"));
            // computing the inner product
            double i_prod = 0.0;  // innter product
            for (Int j = col1; j < col2; j++)
            {
                i_prod += A[(j - col1) * rowCount + i] * x[j + n1];
            }
            Int r = Ps[frowList[i]] + n1;
            x[r] -= i_prod;
        }
    }

#ifndef NDEBUG
    Int m = LUsym->m;
    PRLEVEL(1, ("%% after lsolve x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
