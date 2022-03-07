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
#ifndef NDEBUG
    double start_time = omp_get_wtime();
#endif

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
#ifndef NDEBUG
    double time = omp_get_wtime() - start_time;  
    PRLEVEL(1, ("%%solve has been finished in %lf seconds\n", time));
    printf ("%%solve has been finished in %lf seconds\n", time);
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

ParU_ResultCode paru_solve(paru_matrix *paruMatInfo, double *B, Int n)
{
    DEBUGLEVEL(1);
    PRLEVEL(1, ("%% mRHS inside Solve\n"));
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int m = LUsym->m;
    if (paruMatInfo->res == PARU_SINGULAR)
    {
        printf("Paru: the matrix is singular; cannot be solved.\n");
        return PARU_SINGULAR;
    }
#ifndef NDEBUG
    double start_time = omp_get_wtime();
#endif


    double *X = (double *)paru_alloc(m*n, sizeof(double));
    if (X == NULL)
    {
        printf("Paru: memory problem inside Solve\n");
        return PARU_OUT_OF_MEMORY;
    }
    paru_memcpy(X, B, m*n * sizeof(double));

    paru_apply_perm_scale(LUsym->Pfin, LUsym->scale_row, B, X, m, n);

    PRLEVEL(1, ("%%mRHS lsolve\n"));
    paru_lsolve(paruMatInfo, X, n);  // X = L\X
    PRLEVEL(1, ("%%mRHS usolve\n"));
    paru_usolve(paruMatInfo, X, n);                 // X = U\X
    paru_apply_inv_perm(LUsym->Qfill, X, B, m, n);     // B(q) = X

    paru_free(m*n, sizeof(Int), X);
#ifndef NDEBUG
    double time = omp_get_wtime() - start_time;  
    PRLEVEL(-1, ("%% mRHS solve has been finished in %lf seconds\n", time));
    printf ("%%RHS solve has been finished in %lf seconds\n", time);
#endif
    return PARU_SUCCESS;
}
