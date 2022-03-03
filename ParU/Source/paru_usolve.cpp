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
Int paru_usolve(paru_matrix *paruMatInfo, double *x)
{
    DEBUGLEVEL(0);
    // check if input is read
    if (!x) return (0);
#ifndef NDEBUG
    Int PR = 1;
    double start_time = omp_get_wtime();
#endif
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int nf = LUsym->nf;

    Int n1 = LUsym->n1;   // row+col singletons
    Int *Ps = LUsym->Ps;  // row permutation

    paru_fac *LUs = paruMatInfo->partial_LUs;
    paru_fac *Us = paruMatInfo->partial_Us;
    Int *Super = LUsym->Super;

    const Int max_threads = paruMatInfo->paru_max_threads;
    BLAS_set_num_threads(max_threads);

    for (Int f = nf - 1; f >= 0; --f)
    {
        Int *frowList = paruMatInfo->frowList[f];
        Int *fcolList = paruMatInfo->fcolList[f];
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;
        Int colCount = paruMatInfo->fcolCount[f];

        // do dgemv
        // performed on Us
        //I am not calling BLAS_DGEMV while the column permutation is different

        double *A2 = Us[f].p;
        if (A2 != NULL)
        {
            PRLEVEL(1, ("%% usolve: Working on DGEMV\n%%"));
            #pragma omp parallel for
            for (Int i = 0; i < fp; i++)
            {
                PRLEVEL(1, ("%% Usolve: Working on DGEMV\n"));
                // computing the inner product
                double i_prod = 0.0;  // innter product
                for (Int j = 0; j < colCount; j++)
                {
                    i_prod += A2[fp * j + i] * x[fcolList[j] + n1];
                }
                Int r = Ps[frowList[i]] + n1;
                PRLEVEL(2, ("i_prod[%ld]=%lf  r=%ld\n", i, i_prod,  r));
                x[r] -= i_prod;
            }
        }

        Int rowCount = paruMatInfo->frowCount[f];

        double *A1 = LUs[f].p;
        BLAS_INT N = (BLAS_INT)fp;
        BLAS_INT lda = (BLAS_INT)rowCount;
        // performed on LUs
        PRLEVEL(1, ("%% Usolve: Working on DTRSV\n"));
        cblas_dtrsv (CblasColMajor, CblasUpper, CblasNoTrans, 
                CblasNonUnit, N, A1, lda, x+col1+n1, 1);

        PRLEVEL(1, ("%% DTRSV is just finished\n"));
    }

#ifndef NDEBUG
    Int m = LUsym->m;
    PRLEVEL(1, ("%% before singleton x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
    PR = 1;
#endif
    Int cs1 = LUsym->cs1;
    if (cs1 > 0)
    {
        for (Int i = cs1 - 1; i >= 0; i--)
        {
            PRLEVEL(PR, ("i = %ld\n", i));
            Int *Sup = LUsym->ustons.Sup;
            Int *Suj = LUsym->ustons.Suj;
            double *Sux = LUsym->ustons.Sux;
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
#ifndef NDEBUG
    double time = omp_get_wtime() - start_time;  
    PRLEVEL(-1, ("%% usolve took %1.1lf s; after usolve x is:\n%%", time));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
///////////////////////////////// paru_usolve ///multiple mRHS///////////////////
Int paru_usolve(paru_matrix *paruMatInfo, double *X, Int n)
{
    DEBUGLEVEL(1);
    // check if input is read
    if (!X) return (0);
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int m = paruMatInfo->m;
    Int nf = LUsym->nf;
#ifndef NDEBUG
    Int PR = 1;
    double start_time = omp_get_wtime();
    PRLEVEL(1, ("%% mRHS inside USolve X is:\n"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < n; l++)
        {
            PRLEVEL(1, (" %.2lf, ", X[l*m+k]));
            // PRLEVEL(1, (" %.2lf, ", X[k*n+l])); X row major
        }
        PRLEVEL(1, (" \n"));
    }
    PRLEVEL(1, (" \n"));
#endif
    Int n1 = LUsym->n1;   // row+col singletons
    Int *Ps = LUsym->Ps;  // row permutation

    paru_fac *LUs = paruMatInfo->partial_LUs;
    paru_fac *Us = paruMatInfo->partial_Us;
    Int *Super = LUsym->Super;

    const Int max_threads = paruMatInfo->paru_max_threads;
    BLAS_set_num_threads(max_threads);

    for (Int f = nf - 1; f >= 0; --f)
    {
        Int *frowList = paruMatInfo->frowList[f];
        Int *fcolList = paruMatInfo->fcolList[f];
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;
        Int colCount = paruMatInfo->fcolCount[f];

        // do dgemm
        // performed on Us
        //I cannot call BLAS_DGEMM while the column permutation is different

        double *A2 = Us[f].p;
        if (A2 != NULL)
        {
            PRLEVEL(1, ("%% mRHS usolve: Working on DGEMM f=%ld\n%%", f));
            #pragma omp parallel for
            for (Int i = 0; i < fp; i++)
            {
                // computing the inner product
                double i_prod[n] = {0.0};  // inner product
                for (Int j = 0; j < colCount; j++)
                {
                    for (Int l = 0; l < n; l++)
                        i_prod[l] += A2[fp * j + i] * X[l*m + fcolList[j] + n1];
                }
                Int r = Ps[frowList[i]] + n1;
                for (Int l = 0; l < n; l++)
                {
                    PRLEVEL(2, ("i_prod[%ld]=%lf  r=%ld\n", i, i_prod[i],  r));
                    X[l*m+r] -= i_prod[l];
                }
            }
        }

        Int rowCount = paruMatInfo->frowCount[f];

        double *A1 = LUs[f].p;
        BLAS_INT mm = (BLAS_INT)fp;
        BLAS_INT lda = (BLAS_INT)rowCount;
        BLAS_INT nn = (BLAS_INT)n;

        PRLEVEL(1, ("%% mRHS Usolve: Working on DTRSM\n"));
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, 
                CblasNoTrans,  CblasNonUnit, 
                mm, nn, 1, A1, lda, X+col1+n1, mm);
        PRLEVEL(1, ("%% mRHS DTRSM is just finished\n"));
    }

    PRLEVEL(1, ("%% mRHS Usolve working on singletons \n"));
    Int cs1 = LUsym->cs1;
    if (cs1 > 0)
    {
        for (Int i = cs1 - 1; i >= 0; i--)
        {
            PRLEVEL(PR, ("i = %ld\n", i));
            Int *Sup = LUsym->ustons.Sup;
            Int *Suj = LUsym->ustons.Suj;
            double *Sux = LUsym->ustons.Sux;
            ASSERT(Suj != NULL && Sux != NULL && Sup != NULL);
            for (Int p = Sup[i] + 1; p < Sup[i + 1]; p++)
            {
                Int r = Suj[p];
                PRLEVEL(PR, (" r=%ld\n", r));
                #pragma omp simd
                for (Int l = 0; l < n; l++)
                {
                    X[l*m+i] -= Sux[p] * X[l*m+r];
                }
 
            }
            Int diag = Sup[i];
            #pragma omp simd
            for (Int l = 0; l < n; l++)
            {
                X[l*m+i] /= Sux[diag];
            }
        }
    }
#ifndef NDEBUG
    double time = omp_get_wtime() - start_time;  
    PRLEVEL(1, 
            ("%% mRHS usolve took %1.1lfs; after usolve X is:\n", time));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < n; l++)
        {
            PRLEVEL(1, (" %.2lf, ", X[l*m+k]));
            // PRLEVEL(1, (" %.2lf, ", X[k*n+l])); X row major
        }
        PRLEVEL(1, (" \n"));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
