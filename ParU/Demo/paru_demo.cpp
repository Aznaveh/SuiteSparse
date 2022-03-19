/** =========================================================================  /
 * =======================  paru_demo =======================================  /
 * ==========================================================================  /
 * @brief    test to see how to call umfpack symbolic analysis
 *
 * @author Aznaveh
 * */
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
        printf("Paru: input matrix is invalid\n");
        exit(1);
    }

    if (mtype != CHOLMOD_SPARSE)
    {
        printf("Paru: input matrix must be sparse\n");
        exit(1);
    }

    if (A->xtype != CHOLMOD_REAL)
    {
        printf("Paru: input matrix must be real\n");
        exit(1);
    }

    //~~~~~~~~~~~~~~~~~~~Starting computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double my_start_time = omp_get_wtime();

    ParU_Control Control;
    ParU_Ret info;

    info = ParU_Analyze(A, &Sym, &Control);
    if (info != PARU_SUCCESS)
    {
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    printf("Paru: Symbolic factorization is done!\n");
    ParU_Numeric *Num;
    info = ParU_Factorize(A, Sym, &Num, &Control);
    double my_time = omp_get_wtime() - my_start_time;
    if (info != PARU_SUCCESS)
    {
        ParU_Freenum(&Num, &Control);
        ParU_Freesym(&Sym, &Control);
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    printf("Paru: factorization was successfull in %lf seconds.\n", my_time);

    //~~~~~~~~~~~~~~~~~~~Test the results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if 1
    Int m = Sym->m;
    double *b = (double *)malloc(m * sizeof(double));
    for (Int i = 0; i < m; ++i) b[i] = i + 1;
    double resid, norm;
    ParU_Residual(b, resid, norm, A, Num, &Control);
    for (Int i = 0; i < m; ++i) b[i] = i + 1;
    ParU_Backward(b, resid, norm, A, Num, &Control);
    free(b);

    const Int nn = 16;  // number of right handsides
    double *B = (double *)malloc(m * nn * sizeof(double));
    double Res[4];
    for (Int i = 0; i < m; ++i)
        for (Int j = 0; j < nn; ++j) B[j * m + i] = (double)(i + j + 1);
    ParU_Residual(A, Num, B, Res, nn, &Control);
    free(B);
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
    umf_Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    // umf_Control [UMFPACK_STRATEGY] =   UMFPACK_STRATEGY_UNSYMMETRIC;

    Int *Ap = (Int *)A->p;
    Int *Ai = (Int *)A->i;
    double *Ax = (double *)A->x;
    // Int m = A->nrow;
    Int n = A->ncol;
    void *Symbolic, *Numeric;  // Output argument in umf_dl_symbolc;

    status =
        umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, umf_Control, Info);
    if (status < 0)
    {
        umfpack_dl_report_info(umf_Control, Info);
        umfpack_dl_report_status(umf_Control, status);
        printf("umfpack_dl_symbolic failed\n");
        exit(0);
    }
    status =
        umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, umf_Control, Info);
    if (status < 0)
    {
        umfpack_dl_report_info(umf_Control, Info);
        umfpack_dl_report_status(umf_Control, status);
        printf("umfpack_dl_numeric failed\n");
    }

    umf_time = omp_get_wtime() - umf_start_time;

    b = (double *)malloc(m * sizeof(double));
    double *x = (double *)malloc(m * sizeof(double));
    for (Int i = 0; i < m; ++i) b[i] = i + 1;

    double solve_start = omp_get_wtime();
    status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, umf_Control,
                              Info);
    double solve_time = omp_get_wtime() - solve_start;
    free(x);
    free(b);

    umfpack_dl_free_symbolic(&Symbolic);
    umfpack_dl_free_numeric(&Numeric);

    // Writing results to a file
    if (info == PARU_SUCCESS)
    {
        // Num->umf_time = umf_time;
        // Writing LU factors into a file, can be time consuming
        // paru_write(Num, scale, argv[1]);
        FILE *res_file;
        char res_name[] = "../Demo/Res/res.txt";
        res_file = fopen(res_name, "a");
        if (res_file == NULL)
        {
            printf("Par: error in making %s to write the results!\n", res_name);
        }
        fprintf(res_file, "%ld %ld %lf %lf %lf\n", Sym->m, Sym->anz, my_time,
                umf_time, my_time / umf_time);
        fclose(res_file);
    }
    printf("my_time = %lf umf_time=%lf umf_solv_t = %lf ratio = %lf\n", my_time,
           umf_time, solve_time, my_time / umf_time);

#endif
    //~~~~~~~~~~~~~~~~~~~Free Everything~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ParU_Freenum(&Num, &Control);
    ParU_Freesym(&Sym, &Control);

    cholmod_l_free_sparse(&A, cc);
    cholmod_l_finish(cc);
}
