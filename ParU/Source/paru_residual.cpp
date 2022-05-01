////////////////////////////////////////////////////////////////////////////////
//////////////////////////  ParU_Residual //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*! @brief    get a factorized matrix A and a vector b
 *
 *          compute ||Ax -b||
 *
 *           resid = norm1(b-A*x) / norm1(A)
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp" 

 ParU_Ret ParU_Residual (cholmod_sparse *A, double *x, double *b, Int m,
    double &resid, double &anorm, ParU_Control *Control)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside residual\n"));
    if (A == NULL || x == NULL || b == NULL)
    {
        return PARU_INVALID;
    }
    if ( (Int) A->nrow != m)
    {
        return PARU_INVALID;
    }
 

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
        PRLEVEL(1, ("Paru: memory problem inside residual\n"));
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
 * resid = norm1(b-A*x) / norm1(A)
 *
 * */
ParU_Ret ParU_Residual(cholmod_sparse *A, double *X, double *B, Int m, Int nrhs,  
   double &resid, double &anorm, ParU_Control *Control)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% mRHS inside residual\n"));
    if (A == NULL || X == NULL || B == NULL)
    {
        return PARU_INVALID;
    }
    if ( (Int) A->nrow != m)
    {
        return PARU_INVALID;
    }
 

#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("%% mRHS before everything B is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, ("%%"));
        for (Int l = 0; l < nrhs; l++)
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
        PRLEVEL(1, ("Paru: memory problem inside mRHS residual\n"));
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
