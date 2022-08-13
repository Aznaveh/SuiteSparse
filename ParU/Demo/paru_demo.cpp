//  =========================================================================  /
// =======================  paru_demo =======================================  /
// ==========================================================================  /
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*
 * @brief    test to see how to call umfpack symbolic analysis
 *
 * @author Aznaveh
 * */
#include <math.h>
#include <omp.h>

#include "ParU.hpp"

int main(int argc, char **argv)
{
    cholmod_common Common, *cc;
    cholmod_sparse *A;
    ParU_Symbolic *Sym = NULL;

    //~~~~~~~~~Reading the input matrix and test if the format is OK~~~~~~~~~~~~
    // start CHOLMOD
    cc = &Common;
    int mtype;
    cholmod_l_start(cc);

    // A = mread (stdin) ; read in the sparse matrix A
    A = (cholmod_sparse *)cholmod_l_read_matrix(stdin, 1, &mtype, cc);
    if (A == NULL)
    {
        printf("ParU: input matrix is invalid\n");
        exit(1);
    }

    if (mtype != CHOLMOD_SPARSE)
    {
        printf("ParU: input matrix must be sparse\n");
        exit(1);
    }

    if (A->xtype != CHOLMOD_REAL)
    {
        printf("ParU: input matrix must be real\n");
        exit(1);
    }

    //~~~~~~~~~~~~~~~~~~~Starting computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int ver[3]; char date[128];
    ParU_Version (ver ,date);
    printf("ParU %d.%d.%d",ver[0],ver[1],ver[2]);
    printf(" %s\n",date);

    double my_start_time = omp_get_wtime();

    ParU_Control Control;
    ParU_Ret info;

    //Control.umfpack_ordering = UMFPACK_ORDERING_AMD;
    //Control.umfpack_strategy = UMFPACK_STRATEGY_UNSYMMETRIC;
    //Control.umfpack_strategy = UMFPACK_STRATEGY_SYMMETRIC;
    //Control.umfpack_default_singleton = 0;
    Control.paru_max_threads = 6;
    info = ParU_Analyze(A, &Sym, &Control);
    double my_time_analyze = omp_get_wtime() - my_start_time;
    if (info != PARU_SUCCESS)
    {
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    printf("In: %ldx%ld nnz = %ld \n", Sym->m, Sym->n, Sym->anz);
    printf("ParU: Symbolic factorization is done in %lfs!\n", my_time_analyze);
    ParU_Numeric *Num;
    double my_start_time_fac = omp_get_wtime();
    info = ParU_Factorize(A, Sym, &Num, &Control);
    double my_time_fac = omp_get_wtime() - my_start_time_fac;
    if (info != PARU_SUCCESS)
    {
        printf("ParU: factorization was NOT successfull in %lf seconds!",
               my_time_fac);
        if (info == PARU_OUT_OF_MEMORY)
            printf("\nOut of memory\n");
        if (info == PARU_INVALID)
            printf("\nInvalid!\n");
        if (info == PARU_SINGULAR)
            printf("\nSingular!\n");
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        ParU_Freesym(&Sym, &Control);
        return info;
    }
    else
    {
        printf("ParU: factorization was successfull in %lf seconds.\n",
               my_time_fac);
    }

    //~~~~~~~~~~~~~~~~~~~Test the results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Int m = Sym->m;
    double my_time, my_solve_time;
#if 1
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
#endif

    //~~~~~~~~~~~~~~~~~~~End computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Int max_threads = omp_get_max_threads();
    BLAS_set_num_threads(max_threads);

    //~~~~~~~~~~~~~~~~~~~Calling umfpack~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double umf_time = 0;

#if 1
    double umf_start_time = omp_get_wtime();
    double status,           // Info [UMFPACK_STATUS]
           Info[UMFPACK_INFO],  // Contains statistics about the symbolic analysis

           umf_Control[UMFPACK_CONTROL];  // it is set in umfpack_dl_defaults and
    // is used in umfpack_dl_symbolic; if
    // passed NULL it will use the defaults
    umfpack_dl_defaults(umf_Control);
    //umf_Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_AMD;
    //umf_Control [UMFPACK_STRATEGY] =   UMFPACK_STRATEGY_UNSYMMETRIC;
    //umf_Control [UMFPACK_STRATEGY] =   UMFPACK_STRATEGY_SYMMETRIC;
    umf_Control[UMFPACK_SINGLETONS] = Control.umfpack_default_singleton;

    Int *Ap = (Int *)A->p;
    Int *Ai = (Int *)A->i;
    double *Ax = (double *)A->x;
    // Int m = A->nrow;
    Int n = A->ncol;
    void *Symbolic, *Numeric;  // Output argument in umf_dl_symbolc;

    status =
        umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, umf_Control, Info);
    //umf_Control[UMFPACK_PRL] = 0;
    //umfpack_dl_report_info(umf_Control, Info);
    if (status < 0)
    {
        umfpack_dl_report_info(umf_Control, Info);
        umfpack_dl_report_status(umf_Control, status);
        printf("umfpack_dl_symbolic failed\n");
        exit(0);
    }
    double umf_symbolic = omp_get_wtime() - umf_start_time;
    double umf_fac_start= omp_get_wtime();
    status =
        umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, umf_Control, Info);
    //umf_Control[UMFPACK_PRL] = 2;
    //umfpack_dl_report_info(umf_Control, Info);
    //umfpack_dl_report_status(umf_Control, status);
    if (status < 0)
    {
        umfpack_dl_report_info(umf_Control, Info);
        umfpack_dl_report_status(umf_Control, status);
        printf("umfpack_dl_numeric failed\n");
    }

    double umf_time_fac = omp_get_wtime() - umf_fac_start;

    double *b = (double *)malloc(m * sizeof(double));
    double *x = (double *)malloc(m * sizeof(double));
    for (Int i = 0; i < m; ++i) b[i] = i + 1;

    double solve_start = omp_get_wtime();
    status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, umf_Control,
                              Info);
    double umf_solve_time = omp_get_wtime() - solve_start;
    umf_time = omp_get_wtime() - umf_start_time;
    double umf_resid, umf_anorm;
    info = ParU_Residual(A, x, b, m, umf_resid, umf_anorm, &Control);
    printf("UMFPACK Residual is |%.2lf| and anorm is %.2e and rcond is %.2e.\n",
            umf_resid == 0 ? 0 : log10(umf_resid), umf_anorm, Num->rcond);


    free(x);
    free(b);

    umfpack_dl_free_symbolic(&Symbolic);
    umfpack_dl_free_numeric(&Numeric);

    // Writing results to a file
    if (info == PARU_SUCCESS)
    {
        FILE *res_file;
        char res_name[] = "../Demo/Res/res.txt";
        res_file = fopen(res_name, "a");
        if (res_file == NULL)
        {
            printf("Par: error in making %s to write the results!\n", res_name);
        }
        fprintf(res_file, "%ld %ld %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                Sym->m, Sym->anz, 
                my_time_analyze, my_time_fac, my_solve_time, my_time,
                umf_symbolic, umf_time_fac, umf_solve_time, umf_time);
        fclose(res_file);
    }
    printf("my_time = %lf umf_time=%lf umf_solv_t = %lf ratio = %lf\n", my_time,
            umf_time, umf_solve_time, my_time / umf_time);

#endif
    //~~~~~~~~~~~~~~~~~~~Free Everything~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ParU_Freenum(&Num, &Control);
    ParU_Freesym(&Sym, &Control);

    cholmod_l_free_sparse(&A, cc);
    cholmod_l_finish(cc);
}
