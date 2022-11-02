//  =========================================================================  /
// =======================  paru_simple =====================================  /
// ==========================================================================  /
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*
 * @brief   A very simple example of using ParU by a C main program
 *
 * @author Aznaveh
 * */

#include "ParU_C.h"

int main(int argc, char **argv)
{
    cholmod_common Common, *cc;
    cholmod_sparse *A;
    ParU_C_Symbolic *Sym;
    /* start CHOLMOD */
    cc = &Common;
    int mtype;
    cholmod_l_start(cc);
    /* load A */
    A = (cholmod_sparse *)cholmod_l_read_matrix(stdin, 1, &mtype, cc);

    ParU_C_Control Control;
    init_control (&Control);
    ParU_Ret info;
    info = ParU_C_Analyze(A, &Sym, &Control);
    ParU_C_Numeric *Num;
    info = ParU_C_Factorize(A, Sym, &Num, &Control);
    //FIXME: I am not sure how to deal with it yet
    //~~~~~~~~~~~~~~~~~~~Test the results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Int m = Sym->m;
    double my_time, my_solve_time;
    if (info == PARU_SUCCESS)
    {
        double *b = (double *)malloc(m * sizeof(double));
        double *xx = (double *)malloc(m * sizeof(double));
        for (Int i = 0; i < m; ++i) b[i] = i + 1;
        double my_solve_time_start = omp_get_wtime();
        info = ParU_Solve(Sym, Num, b, xx, &Control);
        if (info != PARU_SUCCESS)
        {
            printf("ParU: Solve has a problem.\n");
            free(b);
            free(xx);
            cholmod_l_free_sparse(&A, cc);
            cholmod_l_finish(cc);
            ParU_Freesym(&Sym, &Control);
            return info;
        }
        my_solve_time = omp_get_wtime() - my_solve_time_start;
        my_time = omp_get_wtime() - my_start_time;
        printf("Solve time is %lf seconds.\n", my_solve_time);
        double resid, anorm;
        info = ParU_Residual(A, xx, b, m, resid, anorm, &Control);
        if (info != PARU_SUCCESS)
        {
            printf("ParU: Residual has a problem.\n");
            free(b);
            free(xx);
            cholmod_l_free_sparse(&A, cc);
            cholmod_l_finish(cc);
            ParU_Freesym(&Sym, &Control);
            return info;
        }

        printf("Residual is |%.2lf| and anorm is %.2e and rcond is %.2e.\n",
                resid == 0 ? 0 : log10(resid), anorm, Num->rcond);

        free(b);
        free(xx);
        const Int nrhs = 16;  // number of right handsides
        double *B = (double *)malloc(m * nrhs * sizeof(double));
        double *X = (double *)malloc(m * nrhs * sizeof(double));
        for (Int i = 0; i < m; ++i)
            for (Int j = 0; j < nrhs; ++j) B[j * m + i] = (double)(i + j + 1);

        info = ParU_Solve(Sym, Num, nrhs, B, X, &Control);
        if (info != PARU_SUCCESS)
        {
            printf("ParU: mRhs Solve has a problem.\n");
            free(B);
            free(X);
            cholmod_l_free_sparse(&A, cc);
            cholmod_l_finish(cc);
            ParU_Freesym(&Sym, &Control);
            return info;
        }
        info = ParU_Residual(A, X, B, m, nrhs, resid, anorm, &Control);
        if (info != PARU_SUCCESS)
        {
            printf("ParU: mRhs Residual has a problem.\n");
            free(B);
            free(X);
            cholmod_l_free_sparse(&A, cc);
            cholmod_l_finish(cc);
            ParU_Freesym(&Sym, &Control);
            return info;
        }

        printf("mRhs Residual is |%.2lf|\n", resid == 0 ? 0 : log10(resid));

        free(B);
        free(X);
    }

    ParU_C_Freenum(&Num, &Control);
    ParU_C_Freesym(&Sym, &Control);
    cholmod_l_free_sparse(&A, cc);
    cholmod_l_finish(cc);
}
