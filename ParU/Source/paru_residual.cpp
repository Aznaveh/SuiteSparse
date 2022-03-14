////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_residual //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief    get a factorized matrix A and a vector b
 *
 *          compute ||Ax -b||
 *
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_Res paru_residual(double *b, double &resid, double &norm,        
        cholmod_sparse *A, paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside residual\n"));
    ParU_symbolic *Sym = paruMatInfo->Sym;
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
    paru_memcpy(x, b, m * sizeof(double));

#ifndef NDEBUG
    PRLEVEL(1, ("%% after copying x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif

    ParU_Res info;
    info = paru_solve(x, paruMatInfo);
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
    resid = paru_vec_1norm(b, m);
    PRLEVEL(1, ("%% resid=%lf\n", resid));
    norm = resid / (paru_spm_1norm(A) * paru_vec_1norm(x, m));
    //    PRLEVEL(1, ("Residual is |%.2lf| and weigheted residual is |%.2f|.\n",
    //                resid == 0 ? 0 : log10(resid),
    //                resid == 0 ? 0 :log10(norm)));
    //
    printf("Residual is |%.2lf| and weigheted residual is |%.2f|.\n",
           resid == 0 ? 0 : log10(resid), resid == 0 ? 0 : log10(norm));

    paru_free(m, sizeof(Int), x);
    return PARU_SUCCESS;
}

//////////////////////////  paru_residual //////////////mRHS////////////////////
/*! @brief  get a factorized matrix A and a  multiple right hand side matrix B 
 *
 *          compute ||Ax -B||
 *
 * @author Aznaveh  for testing now
 * */
ParU_Res paru_residual(cholmod_sparse *A, paru_matrix *paruMatInfo,
                              double *B,
                              double *Results, Int n)  // output
                                                //  0 residual
                                                //  1 weighted residual
                                                //  2 time
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% mRHS inside residual\n"));
    ParU_symbolic *Sym = paruMatInfo->Sym;
    Int m = Sym->m;
#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("%% mRHS before everything B is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < n; l++)
        {
            PRLEVEL(1, (" %.2lf, ", B[l*m+k]));
            // PRLEVEL(1, (" %.2lf, ", B[k*n+l])); B row-major
        }
        PRLEVEL(1, (" \n"));
    }

#endif
    double *X = (double *)paru_alloc(m*n, sizeof(double));
    if (X == NULL)
    {
        printf("Paru: memory problem inside residual\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(X, B, m*n * sizeof(double));

#ifndef NDEBUG
    PRLEVEL(1, ("%% mRHS after copying X is:\n%%"));
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

    ParU_Res info;
    info = paru_solve(X, n, paruMatInfo);
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
            PRLEVEL(1, (" %.2lf, ", X[l*m+k]));
            // PRLEVEL(1, (" %.2lf, ", X[k*n+l])); X row major
        }
        PRLEVEL(1, (" \n"));
    }
    PRLEVEL(1, (" \n"));

#endif

    double res = 2000.0;  //They don't need to be initialized; just for some 
    double weighted_res = 2000.0; //compiler warinings
    for (Int l = 0; l < n; l++)
    {
        PRLEVEL(1, ("%% gaxpy\n"));
        paru_gaxpy(A, X+m*l, B+m*l, -1);
        res = paru_vec_1norm(B+m*l, m);
        PRLEVEL(1, ("%% res=%lf\n", res));
        weighted_res = 
            res / (paru_spm_1norm(A) * paru_vec_1norm(X+m*l, m));
        printf("%ld |%.2lf| |%.2f|,", l,
                res == 0 ? 0 : log10(res), res == 0 ? 0 : log10(weighted_res));
        if ((l+1)%20 == 0) printf("\n");
    }
    printf("\n");
    paru_free(m*n, sizeof(Int), X);
    Results[0] = res;
    Results[1] = weighted_res;
    return PARU_SUCCESS;
}
