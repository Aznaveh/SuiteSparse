////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_tasked_dgemm //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief      a wrapper around  BLAS_DGEMM for tasked base dgemmed 
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define L 512
#define SMALL 8
void paru_tasked_dgemm(Int f,  BLAS_INT M, BLAS_INT N, BLAS_INT K, 
        double *A, BLAS_INT lda, double *B, BLAS_INT ldb, 
        double beta, double *C, BLAS_INT ldc, paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(1);
    //alpha is always -1  in my DGEMMs
    Int num_active_tasks = paruMatInfo->num_active_tasks;
    Int max_threads = omp_get_max_threads();
    PRLEVEL(1, ("%% XXX DGEMM (%d,%d,%d)%1.1f in %ld {%ld}\n", 
                M, N, K, beta, f, num_active_tasks));
    if (M < SMALL && N < SMALL && K < SMALL)
    //if(0)
    {
        PRLEVEL(1, ("%% SMALL DGEMM (%d,%d,%d) in %ld\n", M, N, K, f));
        for (Int i = 0 ; i < M; i++)
            for (Int j = 0 ; j < N; j++)
            {
                if (beta == 0) 
                    C[i+j*ldc]  = 0; 
                for (Int k = 0 ; k < K; k++)
                {
                    C[i+j*ldc] -= A[i+k*lda]*B[k+j*ldb];
                }
            }

    }
    else if (M < L && N < L)
        //if(1)
    { 
        PRLEVEL(1, ("%% No tasking for DGEMM (%dx%d) in %ld\n", M, N, f));
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
                M, N, K, -1, A, lda, B, ldb, beta, C, ldc);
    }
    else if (num_active_tasks == 1)
    { //using all the threads
        BLAS_set_num_threads(max_threads);
        PRLEVEL(1, ("%% A single task DGEMM (%dx%d) in %ld\n", M, N, f));
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
                M, N, K, -1, A, lda, B, ldb, beta, C, ldc);
      //  BLAS_set_num_threads(1);
 
    }
    else
    {
       #ifdef MKLROOT
        Int my_share = max_threads / (num_active_tasks+1);
        if (my_share == 0 ) my_share = 1;
        PRLEVEL(1, ("%% MKL local threads for DGEMM (%dx%d) in %ld [[%ld]]\n", 
                    M, N, f, my_share));
        //using my share of threads
        mkl_set_num_threads_local(my_share);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
                M, N, K, -1, A, lda, B, ldb, beta, C, ldc);
        mkl_set_num_threads_local(0);
        #else
        PRLEVEL(1, ("%%YES tasking for DGEMM (%dx%d) in %ld \n", M, N, f));
        Int num_col_blocks =  N / L + 1 ;
        Int num_row_blocks =  M / L + 1 ;

        Int len_col = N / num_col_blocks;
        Int len_row = M / num_row_blocks;

        PRLEVEL(1, ("%% col-blocks=%ld,row-blocks=%ld [%ld]\n", 
                    num_col_blocks,num_row_blocks,
                    num_col_blocks*num_row_blocks));
        #pragma omp parallel proc_bind(close)
        #pragma omp single nowait
        {
            for (Int I = 0; I < num_row_blocks; I++)
            {
                BLAS_INT m = (I + 1) == num_row_blocks ? 
                    (M - I * len_row) : len_row;

                for (Int J = 0; J < num_col_blocks; J++)
                {
                    BLAS_INT n = (J + 1) == num_col_blocks ? 
                        (N - J * len_col) : len_col;
                    PRLEVEL(1, ("%% I=%ld J=%ld m=%d n=%d in %ld\n", I, J,
                                m, n, f));
                    #pragma omp task
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
                            m, n,K, -1, 
                            A + (I * len_row), lda,
                            B + (J * len_col * ldb), ldb, 
                            beta, C+ (J * ldc * len_col + I * len_row), ldc);
                }
            }
        }  
        #endif
    }

#ifdef COUNT_FLOPS
    #pragma omp atomic update
    paruMatInfo->flp_cnt_dgemm += (double)2 * M * N * K;
#endif
}
