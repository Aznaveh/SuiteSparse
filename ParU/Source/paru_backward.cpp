/** =========================================================================  /
 * =======================  paru_backward ===================================  /
 * ==========================================================================  /
 * @brief    get a factorized matrix A and a vector x1
 *
 *          compute Ax1=b then solve for Ax2=b
 *          return ||x2-x1||
 *
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_ResultCode paru_backward(cholmod_sparse *A, paru_matrix *paruMatInfo,
                              double *x1,
                              double *Results)  // output
//  0 residual
//  1 weighted residual
//  2 time
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside backward\n"));
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int m = LUsym->m;
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
        printf("Memory problem inside residual\n");
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

    ParU_ResultCode info;
    info = paru_solve(paruMatInfo, b);
    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% A problem happend during factorization\n"));
        return info;
    }

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("x2 = [ "));
    for (Int i = 0; i < MIN(m, 10); ++i) PRLEVEL(PR, ("%lf ", b[i]));
    PRLEVEL(PR, (" ...]\n"));
#endif

    for (Int k = 0; k < m; k++) b[k] -= x1[k];

    double res = paru_vec_1norm(b, m);
    PRLEVEL(1, ("%% res=%lf\n", res));
    double weighted_res = res / (paru_spm_1norm(A) * paru_vec_1norm(x1, m));
    //    PRLEVEL(1, ("backward erroris |%.2lf| and weigheted backward error is
    //    |%.2f|.\n",
    //                res == 0 ? 0 : log10(res),
    //                res == 0 ? 0 :log10(weighted_res)));
    //
    printf("backward error is |%.2lf| and weigheted backward erroris |%.2f|.\n",
           res == 0 ? 0 : log10(res), res == 0 ? 0 : log10(weighted_res));

    paru_free(m, sizeof(Int), b);
    Results[0] = res;
    Results[1] = weighted_res;
    return PARU_SUCCESS;
}
