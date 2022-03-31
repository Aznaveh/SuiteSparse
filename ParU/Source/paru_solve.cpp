////////////////////////////////////////////////////////////////////////////////
//////////////////////////  ParU_Solve /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!  @brief  sovle Ax = b
 *      get a factorized matrix and a right hand side
 *      returns x; overwrites it on b
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_Ret ParU_Solve(ParU_Symbolic *Sym, ParU_Numeric *Num, double *b,
                    ParU_Control *user_Control)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% inside solve\n"));
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
    // making copy of user input and check it
    ParU_Control my_Control = *user_Control;
    {
        Int mem_chunk = my_Control.mem_chunk;
        if (mem_chunk < 1024) my_Control.mem_chunk = 1024 * 1024;
        Int max_threads = PARU_OPENMP_MAX_THREADS;
        if (my_Control.paru_max_threads > 0)
            my_Control.paru_max_threads =
                MIN(max_threads, my_Control.paru_max_threads);
        else
            my_Control.paru_max_threads = max_threads;
    }
    ParU_Control *Control = &my_Control;

    paru_memcpy(x, b, m * sizeof(double), Control);
    paru_apply_perm_scale(Sym->Pfin, Num->Rs, b, x, m);

    PRLEVEL(1, ("%% lsolve\n"));
    paru_lsolve(x, Sym, Num, Control);  // x = L\x
    PRLEVEL(1, ("%% usolve\n"));
    paru_usolve(x, Sym, Num, Control);              // x = U\x
    paru_apply_inv_perm(Sym->Qfill, x, b, m);  // b(q) = x

    paru_free(m, sizeof(Int), x);
#ifndef NTIME
    double time = PARU_OPENMP_GET_WTIME;
    time -= start_time;
    PRLEVEL(-1, ("%%solve has been finished in %lf seconds\n", time));
#endif
    return PARU_SUCCESS;
}
////////////////////////Keeps a copy of b //////////////////////////////////////
ParU_Ret ParU_Solve(ParU_Symbolic *Sym, ParU_Numeric *Num, double *b, double *x,
                    ParU_Control *user_Control)
{
    if (Num->res == PARU_SINGULAR)
    {
        printf("Paru: the matrix is singular; cannot be solved.\n");
        return PARU_SINGULAR;
    }
    ParU_Control my_Control = *user_Control;
    {
        Int mem_chunk = my_Control.mem_chunk;
        if (mem_chunk < 1024) my_Control.mem_chunk = 1024 * 1024;
        Int max_threads = PARU_OPENMP_MAX_THREADS;
        if (my_Control.paru_max_threads > 0)
            my_Control.paru_max_threads =
                MIN(max_threads, my_Control.paru_max_threads);
        else
            my_Control.paru_max_threads = max_threads;
    }
    ParU_Control *Control = &my_Control;

    Int m = Sym->m;
    paru_memcpy(x, b, m * sizeof(double), Control);
    return ParU_Solve(Sym, Num, x, user_Control);
}
//////////////////////////  ParU_Solve ////////////// mRHS /////////////////////
/*!  @brief  sovle AX = B
 *      get a factorized matrix and several right hand sides
 *      returns X; overwrites it on B
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_Ret ParU_Solve(ParU_Symbolic *Sym, ParU_Numeric *Num, Int nrhs, double *B,
                    ParU_Control *user_Control)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% mRHS inside Solve\n"));
    Int m = Sym->m;
    if (Num->res == PARU_SINGULAR)
    {
        printf("Paru: the matrix is singular; cannot be solved.\n");
        return PARU_SINGULAR;
    }
#ifndef NTIME
    double start_time = PARU_OPENMP_GET_WTIME;
#endif
    double *X = (double *)paru_alloc(m * nrhs, sizeof(double));
    if (X == NULL)
    {
        printf("Paru: memory problem inside Solve\n");
        return PARU_OUT_OF_MEMORY;
    }

    // making copy of user input and check it
    ParU_Control my_Control = *user_Control;
    {
        Int mem_chunk = my_Control.mem_chunk;
        if (mem_chunk < 1024) my_Control.mem_chunk = 1024 * 1024;
        Int max_threads = PARU_OPENMP_MAX_THREADS;
        if (my_Control.paru_max_threads > 0)
            my_Control.paru_max_threads =
                MIN(max_threads, my_Control.paru_max_threads);
        else
            my_Control.paru_max_threads = max_threads;
    }
    ParU_Control *Control = &my_Control;

    paru_memcpy(X, B, m * nrhs * sizeof(double), Control);
    paru_apply_perm_scale(Sym->Pfin, Num->Rs, B, X, m, nrhs);

    PRLEVEL(1, ("%%mRHS lsolve\n"));
    paru_lsolve(X, nrhs, Sym, Num, Control);  // X = L\X
    PRLEVEL(1, ("%%mRHS usolve\n"));
    paru_usolve(X, nrhs, Sym, Num, Control);              // X = U\X
    paru_apply_inv_perm(Sym->Qfill, X, B, m, nrhs);  // B(q) = X

    paru_free(m * nrhs, sizeof(Int), X);

#ifndef NTIME
    double time = PARU_OPENMP_GET_WTIME;
    time -= start_time;
    PRLEVEL(-1, ("%% mRHS solve has been finished in %lf seconds\n", time));
#endif
    return PARU_SUCCESS;
}
