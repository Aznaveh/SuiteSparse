////////////////////////////////////////////////////////////////////////////////
//////////////////////////  ParU_Residual //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief    get a factorized matrix A and a vector b
 *
 *          compute ||Ax -b||
 *
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp" 
ParU_Ret ParU_Residual(double *b, double &resid, double &norm,
                       cholmod_sparse *A, ParU_Numeric *Num,
                       ParU_Control *Control)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside residual\n"));
    ParU_Symbolic *Sym = Num->Sym;
    Int m = Sym->m;
#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("%% before everything b is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(PR, (" %.2lf, ", b[k]));
    }
    PRLEVEL(PR, (" \n"));
#endif
    double *x = (double *)paru_alloc(m, sizeof(double));
    if (x == NULL)
    {
        printf("Paru: memory problem inside residual\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(x, b, m * sizeof(double), Control);

#ifndef NDEBUG
    PRLEVEL(1, ("%% after copying x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif

    ParU_Ret info;
    info = ParU_Solve(Sym, Num, x, Control);
    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% A problem happend during factorization\n"));
        return info;
    }

#ifndef NDEBUG
    PR = -1;
    PRLEVEL(PR, ("x = [ "));
    for (Int i = 0; i < MIN(m, 10); ++i) PRLEVEL(PR, ("%lf ", x[i]));
    PRLEVEL(PR, (" ...]\n"));
#endif

    PRLEVEL(1, ("%% gaxpy\n"));
    paru_gaxpy(A, x, b, -1);
    resid = paru_vec_1norm(b, m); // |Ax-b|
    PRLEVEL(1, ("%% resid=%lf\n", resid));
    norm = resid / (paru_spm_1norm(A) * paru_vec_1norm(x, m));
    //    PRLEVEL(1, ("Residual is |%.2lf| and weigheted residual is |%.2f|.\n",
    //                resid == 0 ? 0 : log10(resid),
    //                resid == 0 ? 0 :log10(norm)));
    //
    printf("Residual is |%.2lf| and weigheted residual is |%.2f|.\n",
           resid == 0 ? 0 : log10(resid), resid == 0 ? 0 : log10(norm));

    paru_free(m, sizeof(double), x);
    return PARU_SUCCESS;
}
// resid = norm1(b-A*x) / norm1(A)
 ParU_Ret ParU_Residual (cholmod_sparse *A, double *x, double *b, Int m,
    double &resid, double &anorm, ParU_Control *Control)
{
    //FIXME check for the results
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside residual\n"));
#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("%% before everything b is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(PR, (" %.2lf, ", b[k]));
    }
    PRLEVEL(PR, (" \n"));
#endif
    double *ax_b = (double *)paru_alloc(m, sizeof(double));
    if (ax_b == NULL)
    {
        printf("Paru: memory problem inside residual\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(ax_b, b, m * sizeof(double), Control);

#ifndef NDEBUG
    PRLEVEL(1, ("%% after copying x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", ax_b[k]));
    }
    PRLEVEL(1, (" \n"));
#endif


    PRLEVEL(1, ("%% gaxpy\n"));
    paru_gaxpy(A, x, ax_b, -1);
    anorm = paru_spm_1norm(A) ;
    PRLEVEL(1, ("%% resid=%lf\n", resid));
    //resid =  paru_vec_1norm(ax_b, m)/ (anorm* paru_vec_1norm(b, m));
    resid =  paru_vec_1norm(ax_b, m)/ anorm;
    paru_free(m, sizeof(double), ax_b);
    return PARU_SUCCESS;
}


//////////////////////////  ParU_Residual //////////////mRHS////////////////////
/*! @brief  get a factorized matrix A and a  multiple right hand side matrix B
 *
 *          compute ||Ax -B||
 *
 * @author Aznaveh  for testing now
 * */
ParU_Ret ParU_Residual(cholmod_sparse *A, ParU_Numeric *Num, double *B,
                       double *Results, Int n, ParU_Control *Control)  // output
                                                //  0 residual
                                                //  1 weighted residual
                                                //  2 time
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% mRHS inside residual\n"));
    ParU_Symbolic *Sym = Num->Sym;
    Int m = Sym->m;
#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("%% mRHS before everything B is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < n; l++)
        {
            PRLEVEL(1, (" %.2lf, ", B[l * m + k]));
            // PRLEVEL(1, (" %.2lf, ", B[k*n+l])); B row-major
        }
        PRLEVEL(1, (" \n"));
    }

#endif
    double *X = (double *)paru_alloc(m * n, sizeof(double));
    if (X == NULL)
    {
        printf("Paru: memory problem inside residual\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(X, B, m * n * sizeof(double), Control);

#ifndef NDEBUG
    PRLEVEL(1, ("%% mRHS after copying X is:\n%%"));
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

    ParU_Ret info;
    info = ParU_Solve(Sym, Num, n, X, Control);
    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% A problem happend during factorization\n"));
        return info;
    }

#ifndef NDEBUG
    PRLEVEL(1, ("%% mRHS after solve X is:\n%%"));
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

    double res = 2000.0;  // They don't need to be initialized; just for some
    double weighted_res = 2000.0;  // compiler warinings
    for (Int l = 0; l < n; l++)
    {
        PRLEVEL(1, ("%% gaxpy\n"));
        paru_gaxpy(A, X + m * l, B + m * l, -1);
        res = paru_vec_1norm(B + m * l, m);
        PRLEVEL(1, ("%% res=%lf\n", res));
        weighted_res = res / (paru_spm_1norm(A) * paru_vec_1norm(X + m * l, m));
        printf("%ld |%.2lf| |%.2f|,", l, res == 0 ? 0 : log10(res),
               res == 0 ? 0 : log10(weighted_res));
        if ((l + 1) % 20 == 0) printf("\n");
    }
    printf("\n");
    paru_free(m * n, sizeof(double), X);
    Results[0] = res;
    Results[1] = weighted_res;
    return PARU_SUCCESS;
}
// /////////////////////////////////////////////////////////////////////////////
// resid = norm1(b-A*x) / norm1(A)
ParU_Ret ParU_Residual(
   // inputs:
   cholmod_sparse *A, double *X, double *B,
   Int m, Int nrhs,  
   // output:
   double &resid, double &anorm,
   // control:
   ParU_Control *Control)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% mRHS inside residual\n"));
#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("%% mRHS before everything B is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < n; l++)
        {
            PRLEVEL(1, (" %.2lf, ", B[l * m + k]));
            // PRLEVEL(1, (" %.2lf, ", B[k*n+l])); B row-major
        }
        PRLEVEL(1, (" \n"));
    }

#endif
    double *AX_B = (double *)paru_alloc(m * nrhs, sizeof(double));
    if (AX_B == NULL)
    {
        printf("Paru: memory problem inside residual\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(AX_B, B, m * nrhs * sizeof(double), Control);

#ifndef NDEBUG
    PRLEVEL(1, ("%% mRHS after copying X is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < nrhs; l++)
        {
            PRLEVEL(1, (" %.2lf, ", X[l * m + k]));
            // PRLEVEL(1, (" %.2lf, ", X[k*n+l])); X row major
        }
        PRLEVEL(1, (" \n"));
    }
    PRLEVEL(1, (" \n"));
#endif
    anorm = paru_spm_1norm(A) ;
    resid = 0;
    for (Int l = 0; l < nrhs; l++)
    {
        ///TODO I need to use spmv instead of gaxpy here?
        PRLEVEL(1, ("%% gaxpy\n"));
        paru_gaxpy(A, X + m * l, AX_B + m * l, -1);
        //    TODO: I also need to compute 1norm of a dense matrix?
        double res = paru_vec_1norm(AX_B + m * l, m);
        PRLEVEL(1, ("%% res=%lf\n", res));
        resid = MAX(resid, res / anorm);
    }
    paru_free(m*nrhs, sizeof(double), AX_B);
    return PARU_SUCCESS;
}
