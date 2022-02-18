////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_tasked_trsm  //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief      a wrapper around  BLAS_TRSM for tasking
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define L 4096
void paru_tasked_trsm(Int f, int m, int n, double alpha, double *a, int lda,
                      double *b, int ldb, paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    Int naft;
    #pragma omp atomic read
    naft = paruMatInfo->naft;
    const Int max_threads = paruMatInfo->paru_max_threads;
    if (naft == 1)
        BLAS_set_num_threads(max_threads);
    else
        BLAS_set_num_threads(1);
    if ( n < L || (naft == 1) || (naft >= max_threads)) 
    //if(1)
    {
#ifndef NDEBUG
        if (n < L)
            PRLEVEL(1, ("%% Small TRSM (%dx%d) in %ld\n", m, n, f));
        if (naft == 1)
            PRLEVEL(1, ("%% All threads for trsm(%dx%d) in %ld\n", m, n, f));
#endif
        cblas_dtrsm (CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, 
                CblasUnit, m, n, alpha, a, lda, b, ldb);
    }
    else
    {
       #ifdef MKLROOT
        Int my_share = max_threads / naft;
        if (my_share == 0 ) my_share = 1;
        PRLEVEL(1, ("%% MKL local threads for trsm(%dx%d) in %ld [[%ld]]\n", 
                    m, n, f, my_share));
        mkl_set_num_threads_local(my_share);
        cblas_dtrsm (CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, 
                CblasUnit, m, n, alpha, a, lda, b, ldb);
        mkl_set_num_threads_local(0);
        #else
        PRLEVEL(1, ("%%YES tasksingt for trsm(%dx%d) in %ld \n", m, n, f));
        Int num_blocks = n / L + 1;
        Int len_bloc = n / num_blocks;
        PRLEVEL(1, ("%%  num_blocks = %ld", num_blocks));
        #pragma omp parallel proc_bind(close)
        #pragma omp single nowait
        {
            for (Int J = 0; J < num_blocks; J++)
            {
                BLAS_INT n_b = (J + 1) == num_blocks  ? 
                    (n - J * len_bloc) : len_bloc;
                PRLEVEL(1, ("%%  n_b= %d\n", n_b));
                #pragma omp task 
                cblas_dtrsm (CblasColMajor, CblasLeft, CblasLower, 
                        CblasNoTrans, CblasUnit, m, n_b, alpha, a, lda, 
                        (b + J * len_bloc* ldb), ldb);
            }
        }  
        #endif
    }

#ifdef COUNT_FLOPS
    #pragma omp atomic update
    paruMatInfo->flp_cnt_trsm += (double)(m + 1) * m * n;
#endif
}
