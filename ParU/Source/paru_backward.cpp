////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_backward //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*! @brief     compute the backward error
 *
 *          get a factorized matrix A and a vector x1
 *          compute Ax1=b then solve for Ax2=b
 *          return ||x2-x1||
 *
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_Ret paru_backward(double *x1, double &resid, double &norm,
                       cholmod_sparse *A, ParU_Symbolic *Sym, ParU_Numeric *Num,
                       ParU_Control *Control)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside backward\n"));
    if(Sym == NULL || Num == NULL || x1 == NULL || A == NULL)
        return PARU_INVALID;

    Int m = Sym->m;
#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("%% before everything x1 is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(PR, (" %.2lf, ", x1[k]));
    }
    PRLEVEL(PR, (" \n"));
#endif
    double *b = (double *)paru_calloc(m, sizeof(double));
    if (b == NULL)
    {
        PRLEVEL(1, ("Paru: memory problem inside backward\n"));
        return PARU_OUT_OF_MEMORY;
    }
    paru_gaxpy(A, x1, b, 1);
#ifndef NDEBUG
    PRLEVEL(1, ("%% after gaxpy b is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", b[k]));
    }
    PRLEVEL(1, (" \n"));
#endif

    ParU_Ret info;
    info = ParU_Solve(Sym, Num, b, Control);
    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% A problem happend during factorization\n"));
        paru_free(m, sizeof(Int), b);
        return info;
    }

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("x2 = [ "));
    for (Int i = 0; i < MIN(m, 10); ++i) PRLEVEL(PR, ("%lf ", b[i]));
    PRLEVEL(PR, (" ...]\n"));
#endif

    for (Int k = 0; k < m; k++) b[k] -= x1[k];

    resid = paru_vec_1norm(b, m);
    PRLEVEL(1, ("%% resid =%lf\n", resid));
    norm = resid / (paru_spm_1norm(A) * paru_vec_1norm(x1, m));
    PRLEVEL(1, ("backward error is |%.2lf| and weigheted backward error is"
       "|%.2f|.\n",resid == 0 ? 0 : log10(resid), resid == 0 ? 0 :log10(norm)));
    paru_free(m, sizeof(Int), b);
    return PARU_SUCCESS;
}
