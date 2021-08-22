/** =========================================================================  /
 * =======================  paru_solve ======================================  /
 * ==========================================================================  /
 * @brief  sovle Ax = b
 *      get a factorized matrix and a right hand side
 *      returns x
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_ResultCode paru_solve(paru_matrix *paruMatInfo, double *b)
{
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int m = LUsym->m;
    double *x = (double *)paru_alloc(m, sizeof(double));
    if (x == NULL)
    {
        printf("Memory problem inside solve\n");
        return PARU_OUT_OF_MEMORY;
    }

    paru_apply_perm(LUsym->Pfin, b, x, m);  // x = b (p)
    if (paruMatInfo->scale_row)             // x = s.Ps[x]
        paru_apply_scale(paruMatInfo->scale_row, LUsym->Ps, x, m, LUsym->n1);

    paru_lsolve(paruMatInfo, x);                 // x = L\x
    paru_usolve(paruMatInfo, x);                 // x = U\x
    paru_apply_inv_perm(LUsym->Qfill, x, b, m);  // b(q) = x

    return PARU_SUCCESS;
}
