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
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside solve\n"));
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int m = LUsym->m;
    if (paruMatInfo->res == PARU_SINGULAR)
    {
        printf("Paru: the matrix is singular; cannot be solved.\n");
        return PARU_SINGULAR;
    }
    double *x = (double *)paru_alloc(m, sizeof(double));
    if (x == NULL)
    {
        printf("Paru: memory problem inside solve\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(x, b, m * sizeof(double));

    // if (LUsym->scale_row)             // x = s.x
    //    paru_apply_scale(LUsym->scale_row, x, m);
    // paru_memcpy(b, x, m * sizeof(double));
    // paru_apply_perm(LUsym->Pfin, b, x, m);  // x = b (p)
    paru_apply_perm_scale(LUsym->Pfin, LUsym->scale_row, b, x, m);

    PRLEVEL(1, ("%% lsolve\n"));
    paru_lsolve(paruMatInfo, x);  // x = L\x
    PRLEVEL(1, ("%% usolve\n"));
    paru_usolve(paruMatInfo, x);                 // x = U\x
    paru_apply_inv_perm(LUsym->Qfill, x, b, m);  // b(q) = x

    paru_free(m, sizeof(Int), x);
    return PARU_SUCCESS;
}
