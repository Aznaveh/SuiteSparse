/** =========================================================================  /
 * =======================  paru_demo =======================================  /
 * ==========================================================================  /
 * @brief    test to see how to call umfpack symbolic analysis
 *
 * @author Aznaveh
 * */
#include "ParU.hpp"
#include <omp.h>

int main(int argc, char **argv)
{
    cholmod_common Common, *cc;
    cholmod_sparse *A;
    paru_symbolic *LUsym = NULL;

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

    BLAS_set_num_threads(1);

    double my_start_time = omp_get_wtime();

    ParU_ResultCode info;

    info = paru_analyze(A, &LUsym);
    if (info != PARU_SUCCESS)
    {
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    printf ("Paru: Symbolic factorization is done!\n");
    paru_matrix *paruMatInfo;
    info = paru_factorize(A, LUsym, &paruMatInfo);
    double my_time = omp_get_wtime() - my_start_time;
    if (info != PARU_SUCCESS)
    {
        paru_freemat(&paruMatInfo);
        paru_freesym(&LUsym);
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    printf ("Paru: factorization was successfull in %lf seconds.\n",my_time);

    //~~~~~~~~~~~~~~~~~~~Test the results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if 1
    Int m = LUsym->m;
    double *b = (double *)malloc(m * sizeof(double));
    for (Int i = 0; i < m; ++i) b[i] = i + 1;
    double resid, norm;
    paru_residual(b, resid, norm, A, paruMatInfo);
    //for (Int i = 0; i < m; ++i) b[i] = i + 1;
    //paru_backward(A, paruMatInfo, b, Res);
    free(b);

    const Int nn = 16; //number of right handsides
    double *B = (double *)malloc(m*nn * sizeof(double));
    double Res[4];
    for (Int i = 0; i < m; ++i) 
        for (Int j = 0; j < nn; ++j) B[j*m+i] = (double) (i+j + 1);
    paru_residual(A, paruMatInfo, B, Res, nn);
    free(B);
#endif

    //~~~~~~~~~~~~~~~~~~~End computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Int max_threads = paruMatInfo->paru_max_threads;
    BLAS_set_num_threads(max_threads);

    //~~~~~~~~~~~~~~~~~~~Calling umfpack~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double umf_time = 0;

#if 1
    double umf_start_time = omp_get_wtime();
    double status,           // Info [UMFPACK_STATUS]
        Info[UMFPACK_INFO],  // Contains statistics about the symbolic analysis

        Control[UMFPACK_CONTROL];  // it is set in umfpack_dl_defaults and
    // is used in umfpack_dl_symbolic; if
    // passed NULL it will use the defaults
    umfpack_dl_defaults(Control);
    Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    // Control [UMFPACK_STRATEGY] =   UMFPACK_STRATEGY_UNSYMMETRIC;

    Int *Ap = (Int *)A->p;
    Int *Ai = (Int *)A->i;
    double *Ax = (double *)A->x;
    // Int m = A->nrow;
    Int n = A->ncol;
    void *Symbolic, *Numeric;  // Output argument in umf_dl_symbolc;

    status = umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control, Info);
    if (status < 0)
    {
        umfpack_dl_report_info(Control, Info);
        umfpack_dl_report_status(Control, status);
        printf("umfpack_dl_symbolic failed\n");
        exit(0);
    }
    status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
    if (status < 0)
    {
        umfpack_dl_report_info(Control, Info);
        umfpack_dl_report_status(Control, status);
        printf("umfpack_dl_numeric failed\n");
    }

    umf_time = omp_get_wtime() - umf_start_time;
    umfpack_dl_free_symbolic(&Symbolic);
    umfpack_dl_free_numeric(&Numeric);

    // Writing results to a file
    if (info == PARU_SUCCESS)
    {
        // paruMatInfo->umf_time = umf_time;
        // Writing LU factors into a file, can be time consuming
        // paru_write(paruMatInfo, scale, argv[1]);
        FILE *res_file;
        char res_name[] = "../Demo/Res/res.txt";
        res_file = fopen(res_name, "a");
        if (res_file == NULL)
        {
            printf("Par: error in making %s to write the results!\n",
                   res_name);
        }
        fprintf(res_file, "%ld %ld %lf %lf %lf\n", LUsym->m, LUsym->anz,
                my_time, umf_time, my_time / umf_time);
        fclose(res_file);
    }
    printf("my_time = %lf umf_time=%lf ratio = %lf\n", my_time, umf_time,
           my_time / umf_time);

#endif
    //~~~~~~~~~~~~~~~~~~~Free Everything~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    paru_freemat(&paruMatInfo);
    paru_freesym(&LUsym);

    cholmod_l_free_sparse(&A, cc);
    cholmod_l_finish(cc);
}
