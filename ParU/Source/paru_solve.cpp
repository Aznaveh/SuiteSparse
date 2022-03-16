////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_solve /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!  @brief  sovle Ax = b
 *      get a factorized matrix and a right hand side
 *      returns x; overwrites it on b
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_Ret paru_solve(double *b, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside solve\n"));
    ParU_Symbolic *Sym = Num->Sym;
    Int m = Sym->m;
    if (Num->res == PARU_SINGULAR)
    {
        printf("Paru: the matrix is singular; cannot be solved.\n");
        return PARU_SINGULAR;
    }
#ifndef NTIME
    double start_time = PARU_OPENMP_GET_WTIME;
#endif

    double *x = (double *)paru_alloc(m, sizeof(double));
    if (x == NULL)
    {
        printf("Paru: memory problem inside solve\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(x, b, m * sizeof(double));

    // if (Sym->scale_row)             // x = s.x
    //    paru_apply_scale(Sym->scale_row, x, m);
    // paru_memcpy(b, x, m * sizeof(double));
    // paru_apply_perm(Sym->Pfin, b, x, m);  // x = b (p)
    paru_apply_perm_scale(Sym->Pfin, Sym->scale_row, b, x, m);

    PRLEVEL(1, ("%% lsolve\n"));
    paru_lsolve(x, Num);  // x = L\x
    PRLEVEL(1, ("%% usolve\n"));
    paru_usolve(x, Num);                       // x = U\x
    paru_apply_inv_perm(Sym->Qfill, x, b, m);  // b(q) = x

    paru_free(m, sizeof(Int), x);
#ifndef NTIME
    double time = PARU_OPENMP_GET_WTIME;
    time -= start_time;
    PRLEVEL(-1, ("%%solve has been finished in %lf seconds\n", time));
#endif
    return PARU_SUCCESS;
}
//////////////////////////  paru_solve ////////////// mRHS /////////////////////
/*!  @brief  sovle AX = B
 *      get a factorized matrix and several right hand sides
 *      returns X; overwrites it on B
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_Ret paru_solve(double *B, Int n, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% mRHS inside Solve\n"));
    ParU_Symbolic *Sym = Num->Sym;
    Int m = Sym->m;
    if (Num->res == PARU_SINGULAR)
    {
        printf("Paru: the matrix is singular; cannot be solved.\n");
        return PARU_SINGULAR;
    }
#ifndef NTIME
    double start_time = PARU_OPENMP_GET_WTIME;
#endif
    double *X = (double *)paru_alloc(m * n, sizeof(double));
    if (X == NULL)
    {
        printf("Paru: memory problem inside Solve\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(X, B, m * n * sizeof(double));

    paru_apply_perm_scale(Sym->Pfin, Sym->scale_row, B, X, m, n);

    PRLEVEL(1, ("%%mRHS lsolve\n"));
    paru_lsolve(X, n, Num);  // X = L\X
    PRLEVEL(1, ("%%mRHS usolve\n"));
    paru_usolve(X, n, Num);                       // X = U\X
    paru_apply_inv_perm(Sym->Qfill, X, B, m, n);  // B(q) = X

    paru_free(m * n, sizeof(Int), X);

#ifndef NTIME
    double time = PARU_OPENMP_GET_WTIME;
    time -= start_time;
    PRLEVEL(-1, ("%% mRHS solve has been finished in %lf seconds\n", time));
#endif
    return PARU_SUCCESS;
}
