////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// paru_lsolve //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  In this file I have lsovle x = L\x
 * 
 *
 ********    The final result is something like this (nf = 4)                     
 *     col1    col2
 *      ___________________________________________                          
 *     |\       |                                 |       c                 
 *     |*\      |         U1                      |       c                 
 *     |****\   |                                 |       c  DTRSV on here  
 *     |*DTRSV*\|_________________________________|       c                 
 *     |******* |\         |                      |       c                 
 *     |        |****\     |       U2             |       x                 
 *     |   LU1  |*DTRSV**\ |______________________|       x DGEMV updates   
 *     |        |          |**\         |         |       x   here          
 *     |        | LU2      |*****\      |   U3    |       x                 
 *     | DGEMV  |          |*DTRSV**\   |_________|       x                 
 *     |        |          |   LU3      |* LU4    |       x                 
 *     |        | DGEMV    |  DGEMV     |****     |       x                 
 *     |        |          |            |DTRSV*   |       x                 
 *     |________|__________|____________|_________|       x                 
 *                                                                          
 *     This function just goes through LUs in the data structure and does a 
 *     TRSV on triangular part                                              
 *     Then does DGEMV on the rest for 0 to nf                              
 *  
 *           BLAS_DTRSV  is used here but I do not use BLAS_DGEMV explicitly
 *           while it needs space for each thread doing this computation.
 *           I guess using this way can have a good performance.
 *  
 * @author Aznaveh
 * */
#include "paru_internal.hpp"
Int paru_lsolve(paru_matrix *paruMatInfo, double *x)
{
    DEBUGLEVEL(1);
    if (!x) return (0);
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int nf = LUsym->nf;

#ifndef NDEBUG
    Int m = LUsym->m;
    Int PR = 1;
    double start_time = omp_get_wtime();
    PRLEVEL(1, ("%%inside lsolve x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
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
#ifndef NDEBUG
    PRLEVEL(PR, ("%%lsove singletons finished and  x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
 
    const Int max_threads = paruMatInfo->paru_max_threads;
    BLAS_set_num_threads(max_threads);
    double work[LUsym->m]; //gather scatter space for dgemv 

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

        BLAS_INT lda = (BLAS_INT)rowCount;
        {

            BLAS_INT N = (BLAS_INT)fp;
            PRLEVEL(2, ("%% Working on DTRSV\n"));
            cblas_dtrsv (CblasColMajor, CblasLower, CblasNoTrans, 
                    CblasUnit, N, A, lda, X, 1);
            PRLEVEL(2, ("%% DTRSV is just finished\n"));
        }
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

        PRLEVEL(PR, ("%% lda = %d\n%%", lda));
        PRLEVEL(PR, ("%% during lsolve x [%ld-%ld)is:\n%%", col1, col2));
        // for (Int k = col1; k < col2; k++)
        Int m = LUsym->m;
        for (Int k = 0; k < m; k++)
        {
            PRLEVEL(PR, (" %.2lf, ", x[k]));
        }
        PRLEVEL(PR, (" \n"));
#endif

        if (rowCount > fp)
        {
            PRLEVEL(2, ("%% lsolve: Working on DGEMV\n%%"));
            PRLEVEL(2, ("fp=%ld  rowCount=%ld\n", fp, rowCount));
            BLAS_INT m = (BLAS_INT)(rowCount-fp);
            BLAS_INT n = (BLAS_INT)fp;
            cblas_dgemv (CblasColMajor, CblasNoTrans, m, n, 1, A+fp, lda, 
                    x+n1+col1, 1, 0, work, 1);
        }

        //don't use parallel loop if using dgemv
        //pragma omp parallel for
        for (Int i = fp; i < rowCount; i++)
        {
            //alternative to dgemv; do not need work if using this
            // computing the inner product
            //double i_prod = 0.0;  // inner product
            //for (Int j = col1; j < col2; j++)
            //{
            //    i_prod += A[(j - col1) * rowCount + i] * x[j + n1];
            //}

            double i_prod = work[i-fp];
            Int r = Ps[frowList[i]] + n1;
            PRLEVEL(2, ("i_prod[%ld]=%lf  work=%lf r=%ld\n", 
                        i, i_prod,  work[i-fp], r));
            x[r] -= i_prod;
        }
    }

#ifndef NDEBUG
    double time = omp_get_wtime() - start_time;  
    PRLEVEL(1, ("%% lsolve took %1.1lf s; after lsolve x is:\n%%", time));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
///////////////////////////////// paru_lsolve ///multiple mRHS///////////////////
Int paru_lsolve(paru_matrix *paruMatInfo, double *X, Int n)
{
    DEBUGLEVEL(1);
    if (!X) return (0);
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int m = LUsym->m;
    Int nf = LUsym->nf;

#ifndef NDEBUG
    Int PR = 1;
    double start_time = omp_get_wtime();
    PRLEVEL(1, ("%% mRHS inside LSolve X is:\n"));
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
            PRLEVEL(PR, (" X[%ld]=%.2lf Slx[%ld]=%.2lf\n", j, X[j*n], diag,
                        Slx[diag]));
            #pragma omp simd
            for (Int l = 0; l < n; l++)
            {
                //  X[j*n+l] /= Slx[diag]; row major
                X[l*m+j] /= Slx[diag];
            }
            PRLEVEL(PR, (" After X[%ld]=%.2lf \n", j, X[j*n]));

            for (Int p = Slp[j - cs1] + 1; p < Slp[j - cs1 + 1]; p++)
            {
                Int r = Sli[p] < n1 ? Sli[p] : Ps[Sli[p] - n1] + n1;
                PRLEVEL(PR, (" r=%ld\n", r));
                #pragma omp simd
                for (Int l = 0; l < n; l++)
                {
                    //X[r*n+l] -= Slx[p] * X[j*n+l];   row major
                    X[l*m+r] -= Slx[p] * X[l*m+j];
                }
                PRLEVEL(PR, ("A X[%ld]=%.2lf\n", Sli[p], X[Sli[p]]*n));
            }
        }
    }
#ifndef NDEBUG
    PRLEVEL(1, ("%% mRHS lsovle singletons finished and X is:\n"));
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
    const Int max_threads = paruMatInfo->paru_max_threads;
    BLAS_set_num_threads(max_threads);
    double work[m*n]; //gather scatter space for dgemm

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
        //double *Xp = X + (col1 + n1)*n; row major
        BLAS_INT lda = (BLAS_INT)rowCount;
        BLAS_INT ldb = (BLAS_INT)m;
        {
            BLAS_INT mm = (BLAS_INT)fp;
            BLAS_INT nn = (BLAS_INT)n;
            PRLEVEL(2, ("%% mRHS Working on DTRSM f=%ld\n", f));
            cblas_dtrsm (CblasColMajor, CblasLeft, CblasLower, 
                    CblasNoTrans,  CblasUnit, 
                    mm, nn, 1, A, lda,  X+n1+col1, ldb);
            PRLEVEL(2, ("%% mRHS DTRSM is just finished f=%ld\n",f));
        }
        #ifndef NDEBUG
        PR = 2;
        PRLEVEL(PR, ("%% LUs:\n%%"));
        for (Int r = 0; r < rowCount; r++)
        {
            PRLEVEL(PR, ("%% %ld\t", frowList[r]));
            for (Int c = col1; c < col2; c++)
                PRLEVEL(PR, (" %2.5lf\t", A[(c - col1)*rowCount + r]));
            PRLEVEL(PR, ("\n"));
        }

        PRLEVEL(PR, ("%% lda = %d\n%%", lda));
        PR = 1;
        PRLEVEL(PR, ("%% during lsolve X f=%ld[%ld-%ld)is:\n%%", 
                    f, col1, col2));
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

        if (rowCount > fp)
        {
            PRLEVEL(2, ("%% mRHS lsolve: Working on DGEMM\n%%"));
            PRLEVEL(2, ("fp=%ld  rowCount=%ld\n", fp, rowCount));
            BLAS_INT mm = (BLAS_INT)(rowCount-fp);
            BLAS_INT kk = (BLAS_INT)fp;
            BLAS_INT nn = (BLAS_INT)n;
            cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, mm, nn, kk,
                    1, A+fp, lda, X+n1+col1, ldb, 0, work, ldb);
 
        }

        //don't use parallel loop if using dgemm
        //pragma omp parallel for
        for (Int i = fp; i < rowCount; i++)
        {
            //alternative to dgemm; do not need work if using this
            // computing the inner product
            //double i_prod[n] = {0.0};  // inner product
            //for (Int j = col1; j < col2; j++)
            //{
            //    for (Int l = 0; l < n; l++)
            //        i_prod[l] += A[(j - col1) * rowCount + i] * X[l*m +j+n1];
            //}
        

            double* i_prod = work+i-fp;
            Int r = Ps[frowList[i]] + n1;
            for (Int l = 0; l < n; l++)
            {
                PRLEVEL(2, ("i_prod[%ld]=%lf  work=%lf r=%ld\n", 
                            i, i_prod[i],  work[i-fp], r));
                X[l*m+r] -= i_prod[l*m];
            }

        }
    }

    #ifndef NDEBUG
    double time = omp_get_wtime() - start_time;  
    PRLEVEL(1, 
            ("%% mRHS lsolve took %1.1lfs; after lsolve X is:\n", time));
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
