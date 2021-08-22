/** =========================================================================  /
 * =======================  paru_residual ===================================  /
 * ==========================================================================  /
 * @brief    get a factorized matrix A and a vector b
 *          compute ||Ax -b ||
 *
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

double paru_residual(cholmod_sparse *A, paru_matrix *paruMatInfo, double *b)
{
    DEBUGLEVEL(1);
    PRLEVEL(1, ("%% inside residual\n"));
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int m = LUsym->m;
#ifndef NDEBUG
    PRLEVEL(1, ("%% before everything b is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", b[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    double *x = (double *)paru_alloc(m, sizeof(double));
    if (x == NULL)
    {
        printf("Memory problem inside residual\n");
        return -1;
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

    ParU_ResultCode info;
    info = paru_solve(paruMatInfo, x);
    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% A problem happend during factorization\n"));
        return -1;
    }

#ifndef NDEBUG
    Int PR = -1;
    PRLEVEL(PR, ("x = [ "));
    for (Int i = 0; i < MIN(m, 10); ++i) PRLEVEL(PR, ("%lf ", x[i]));
    PRLEVEL(PR, (" ...]\n"));
#endif

    PRLEVEL(1, ("%% gaxpy\n"));
    paru_gaxpy(A, x, b, -1);
    double res = paru_vec_1norm(b, m);
    double weighted_res = res / (paru_spm_1norm(A) * paru_vec_1norm(x, m));
    PRLEVEL(1, ("Residual is |%.2lf| and weigheted residual is |%.2f|.\n",
                log10(res), log10(weighted_res)));

    paru_free(m, sizeof(Int), x);
    return res;
}
