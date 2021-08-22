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
    DEBUGLEVEL(1);
    PRLEVEL(1, ("%% inside solve\n"));
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

    PRLEVEL(1, ("%% lsolve\n"));
    paru_lsolve(paruMatInfo, x);                 // x = L\x
    PRLEVEL(1, ("%% usolve\n"));
    paru_usolve(paruMatInfo, x);                 // x = U\x
    paru_apply_inv_perm(LUsym->Qfill, x, b, m);  // b(q) = x

    paru_free(m, sizeof(Int), x);
    return PARU_SUCCESS;
}
