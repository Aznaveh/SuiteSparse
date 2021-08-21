/** =========================================================================  /
 * =======================  paru_solve ======================================  /
 * ==========================================================================  /
 * @brief    get a matrix and factorize it
 *      Allocate space for paruMatInfo
 *      the user should free the space
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_ResultCode paru_solve(Int f, paru_matrix *paruMatInfo)
// This routine call paru_front from first(f)...f including f
// This routine is called recursively to make tasks
{
    // TODO: temporary to test perm
    paru_perm(paruMatInfo);
    Int m = LUsym->m;
    double b[m];
    double x[m];
    double xt[m];
    for (Int i = 0; i < m; ++i) b[i] = i + 1;
    paru_apply_perm(LUsym->Pfin, b, x, m);  // x = p (b)
    if (paruMatInfo->scale_row)
        paru_apply_scale (paruMatInfo->scale_row, LUsym->Ps, x, m, LUsym->n1);
    paru_lsolve(paruMatInfo, x);
    paru_usolve(paruMatInfo, x);
    paru_apply_inv_perm(LUsym->Qfill, x, xt, m);  // xt = qinv (x)


    printf ("x = [ ");
    for (Int i = 0; i < MIN(m,10); ++i) 
        printf ("%lf ",xt[i]);
    printf (" ...]\n");

    for (Int i = 0; i < m; ++i) b[i] *= -1;
    //    b[i] = 0;
    paru_gaxpy(A, xt, b);
    double res = paru_vec_1norm(b, m);
    double weighted_res = res / (paru_spm_1norm(A) * paru_vec_1norm(xt, m));
    printf("Residual is |%.2lf| and weigheted residual is |%.2f|.\n", 
            log10(res), log10(weighted_res) );

}
