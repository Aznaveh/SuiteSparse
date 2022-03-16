////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// paru_usolve //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  In this file I have usovle x = U\x
 *
 *
 ********       The final result is something like this (nf = 4)
 *         ___________________________________________
 *        |\*******|                                 |       x
 *        | \*DTRSV|         U1   DGEMV              |       x
 *        |    \***|                                 |       x
 *        |       \|_________________________________|       x
 *        |        |\**DTRSV**|                      |       x
 *        |        |    \*****|       U2    DGEMV    |       x
 *        |   LU1  |        \*|______________________|       o
 *        |        |          |  \**DTRSV**|  DGEMV  |       o
 *        |        | LU2      |     \******|   U3    |       o
 *        |        |          |        \***|_________|       o DGEMV updates
 *        |        |          |   LU3      |*********|       c  here
 *        |        |          |            |  *******|       c   DTRSV on here
 *        |        |          |            |LU4  ****|       c
 *        |________|__________|____________|________*|       c
 *
 *        This function just goes through LUs and US in the data structure and
 *        does a TRSV on triangular part  Then does DGEMV on the rest
 *       for nf down to 0
 *
 *              BLAS_DTRSV  is used here but I do not use BLAS_DGEMV explicitly
 *              while it needs space for each thread doing this computation.
 *              I guess using this way can have a good performance.
 *
 * @author Aznaveh
 * */
#include "paru_internal.hpp"
Int paru_usolve(double *x, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);
    // check if input is read
    if (!x) return (0);
    PARU_DEFINE_PRLEVEL;
#ifndef NTIME
    double start_time = PARU_OPENMP_GET_WTIME;
#endif
    ParU_Symbolic *Sym = Num->Sym;
    Int nf = Sym->nf;

    Int n1 = Sym->n1;   // row+col singletons
    Int *Ps = Sym->Ps;  // row permutation

    ParU_Factors *LUs = Num->partial_LUs;
    ParU_Factors *Us = Num->partial_Us;
    Int *Super = Sym->Super;

    const Int max_threads = Num->paru_max_threads;
    BLAS_set_num_threads(max_threads);
    std::vector<double> work(Num->max_col_count);

    for (Int f = nf - 1; f >= 0; --f)
    {
        Int *frowList = Num->frowList[f];
        Int *fcolList = Num->fcolList[f];
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;
        Int colCount = Num->fcolCount[f];

        // do dgemv
        // performed on Us
        // I am not calling BLAS_DGEMV while the column permutation is different

        double *A2 = Us[f].p;
        if (A2 != NULL)
        {
            PRLEVEL(2, ("%% usolve: Working on DGEMV\n%%"));

            double *xg = &work[0] + fp;         // size Xg is colCount
            for (Int j = 0; j < colCount; j++)  // gathering x in Xg
            {
                xg[j] = x[fcolList[j] + n1];
            }

            BLAS_INT mm = (BLAS_INT)fp;
            BLAS_INT nn = (BLAS_INT)colCount;
            BLAS_INT lda = (BLAS_INT)fp;
            cblas_dgemv(CblasColMajor, CblasNoTrans, mm, nn, 1, A2, lda, xg, 1,
                        0, &work[0], 1);

            for (Int i = 0; i < fp; i++)  // scattering the back in to x
            {
                Int r = Ps[frowList[i]] + n1;
                x[r] -= work[i];
            }

            // pragma omp parallel for
            // for (Int i = 0; i < fp; i++)
            //{
            //    PRLEVEL(2, ("%% Usolve: Working on DGEMV\n"));
            //    // computing the inner product
            //    double i_prod = 0.0;  // innter product
            //    for (Int j = 0; j < colCount; j++)
            //    {
            //        i_prod += A2[fp * j + i] * x[fcolList[j] + n1];
            //    }
            //    Int r = Ps[frowList[i]] + n1;
            //    PRLEVEL(2, ("i_prod[%ld]=%lf  r=%ld\n", i, i_prod,  r));
            //    x[r] -= i_prod;
            //}
        }

        Int rowCount = Num->frowCount[f];

        double *A1 = LUs[f].p;
        BLAS_INT N = (BLAS_INT)fp;
        BLAS_INT lda = (BLAS_INT)rowCount;
        // performed on LUs
        PRLEVEL(2, ("%% Usolve: Working on DTRSV\n"));
        cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, N,
                    A1, lda, x + col1 + n1, 1);

        PRLEVEL(2, ("%% DTRSV is just finished\n"));
    }

#ifndef NDEBUG
    Int m = Sym->m;
    PRLEVEL(1, ("%% before singleton x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
    PR = 1;
#endif
    Int cs1 = Sym->cs1;
    if (cs1 > 0)
    {
        for (Int i = cs1 - 1; i >= 0; i--)
        {
            PRLEVEL(PR, ("i = %ld\n", i));
            Int *Sup = Sym->ustons.Sup;
            Int *Suj = Sym->ustons.Suj;
            double *Sux = Sym->ustons.Sux;
            ASSERT(Suj != NULL && Sux != NULL && Sup != NULL);
            PRLEVEL(PR, (" Before computation x[%ld]=%.2lf \n", i, x[i]))
            for (Int p = Sup[i] + 1; p < Sup[i + 1]; p++)
            {
                Int r = Suj[p];
                PRLEVEL(PR, (" r=%ld\n", r));
                x[i] -= Sux[p] * x[r];
                PRLEVEL(PR, ("Suj[%ld]=%ld\n", p, Suj[p]));
                PRLEVEL(PR, (" x[%ld]=%.2lf x[%ld]=%.2lf\n", r, x[r], i, x[i]));
            }
            Int diag = Sup[i];
            x[i] /= Sux[diag];
            PRLEVEL(PR, (" After computation x[%ld]=%.2lf \n", i, x[i]))
            PRLEVEL(PR, ("\n"));
        }
    }

#ifndef NTIME
    double time = PARU_OPENMP_GET_WTIME;
    time -= start_time;
    PRLEVEL(-1, ("%% usolve took %1.1lf\n", time));
#endif
#ifndef NDEBUG
    PRLEVEL(-1, ("%%after usolve x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
///////////////////////////////// paru_usolve ///multiple
///mRHS///////////////////
Int paru_usolve(double *X, Int n, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);
    // check if input is read
    if (!X) return (0);
    PARU_DEFINE_PRLEVEL;
    ParU_Symbolic *Sym = Num->Sym;
    Int m = Sym->m;
    Int nf = Sym->nf;
#ifndef NDEBUG
    PRLEVEL(1, ("%% mRHS inside USolve X is:\n"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < n; l++)
        {
            PRLEVEL(1, (" %.2lf, ", X[l * m + k]));
            // PRLEVEL(1, (" %.2lf, ", X[k*n+l])); X row major
        }
        PRLEVEL(1, (" \n"));
    }
    PRLEVEL(1, (" \n"));
#endif
#ifndef NTIME
    double start_time = PARU_OPENMP_GET_WTIME;
#endif
    Int n1 = Sym->n1;   // row+col singletons
    Int *Ps = Sym->Ps;  // row permutation

    ParU_Factors *LUs = Num->partial_LUs;
    ParU_Factors *Us = Num->partial_Us;
    Int *Super = Sym->Super;

    const Int max_threads = Num->paru_max_threads;
    BLAS_set_num_threads(max_threads);
    std::vector<double> work(Num->max_col_count * n);

    for (Int f = nf - 1; f >= 0; --f)
    {
        Int *frowList = Num->frowList[f];
        Int *fcolList = Num->fcolList[f];
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;
        Int colCount = Num->fcolCount[f];

        // do dgemm
        // performed on Us

        double *A2 = Us[f].p;
        if (A2 != NULL)
        {
            PRLEVEL(2, ("%% mRHS usolve: Working on DGEMM f=%ld\n%%", f));
            double *Xg = &work[0] + fp * n;     // size Xg is colCount x n
            for (Int j = 0; j < colCount; j++)  // gathering X in Xg
            {
                for (Int l = 0; l < n; l++)
                {
                    Xg[l * colCount + j] = X[l * m + fcolList[j] + n1];
                }
            }

            BLAS_INT mm = (BLAS_INT)fp;
            BLAS_INT nn = (BLAS_INT)n;
            BLAS_INT kk = (BLAS_INT)colCount;
            BLAS_INT lda = (BLAS_INT)fp;
            BLAS_INT ldb = (BLAS_INT)colCount;
            BLAS_INT ldc = (BLAS_INT)fp;
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, mm, nn, kk,
                        1, A2, lda, Xg, ldb, 0, &work[0], ldc);
            for (Int i = 0; i < fp; i++)  // scattering the back in to X
            {
                Int r = Ps[frowList[i]] + n1;
                for (Int l = 0; l < n; l++)
                {
                    X[l * m + r] -= work[l * fp + i];
                }
            }

            // algternative way for dgemm
            // pragma omp parallel for schedule(static)
            // for (Int i = 0; i < fp; i++)
            //{
            //    // computing the inner product
            //    double i_prod[n] = {0.0};  // inner product
            //    for (Int j = 0; j < colCount; j++)
            //    {
            //        for (Int l = 0; l < n; l++)
            //            i_prod[l] += A2[fp * j + i] * X[l*m + fcolList[j] +
            //            n1];
            //    }
            //    Int r = Ps[frowList[i]] + n1;
            //    for (Int l = 0; l < n; l++)
            //    {
            //        PRLEVEL(2, ("i_prod[%ld]=%lf  r=%ld\n", i, i_prod[i], r));
            //        X[l*m+r] -= i_prod[l];
            //    }
            //}
        }

        Int rowCount = Num->frowCount[f];

        double *A1 = LUs[f].p;
        BLAS_INT mm = (BLAS_INT)fp;
        BLAS_INT lda = (BLAS_INT)rowCount;
        BLAS_INT nn = (BLAS_INT)n;
        BLAS_INT ldb = (BLAS_INT)m;

        PRLEVEL(2, ("%% mRHS Usolve: Working on DTRSM\n"));
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
                    CblasNonUnit, mm, nn, 1, A1, lda, X + col1 + n1, ldb);
        PRLEVEL(2, ("%% mRHS DTRSM is just finished\n"));
    }

    PRLEVEL(1, ("%% mRHS Usolve working on singletons \n"));
    Int cs1 = Sym->cs1;
    if (cs1 > 0)
    {
        for (Int i = cs1 - 1; i >= 0; i--)
        {
            PRLEVEL(PR, ("i = %ld\n", i));
            Int *Sup = Sym->ustons.Sup;
            Int *Suj = Sym->ustons.Suj;
            double *Sux = Sym->ustons.Sux;
            ASSERT(Suj != NULL && Sux != NULL && Sup != NULL);
            for (Int p = Sup[i] + 1; p < Sup[i + 1]; p++)
            {
                Int r = Suj[p];
                PRLEVEL(PR, (" r=%ld\n", r));
#pragma omp simd
                for (Int l = 0; l < n; l++)
                {
                    X[l * m + i] -= Sux[p] * X[l * m + r];
                }
            }
            Int diag = Sup[i];
#pragma omp simd
            for (Int l = 0; l < n; l++)
            {
                X[l * m + i] /= Sux[diag];
            }
        }
    }
#ifndef NTIME
    double time = PARU_OPENMP_GET_WTIME;
    time -= start_time;
    PRLEVEL(-1, ("%% mRHS usolve took %1.1lfs\n", time));
#endif
#ifndef NDEBUG
    PRLEVEL(1, ("%%after usolve X is:\n"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < n; l++)
        {
            PRLEVEL(1, (" %.2lf, ", X[l * m + k]));
            // PRLEVEL(1, (" %.2lf, ", X[k*n+l])); X row major
        }
        PRLEVEL(1, (" \n"));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
