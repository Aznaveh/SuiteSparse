//  =========================================================================  /
// =======================  paru_test =======================================  /
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

#include "paru_cov.hpp"

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
    int ver[3]; char date[128];
    ParU_Version (ver ,date);
    printf("Paru %d.%d.%d",ver[0],ver[1],ver[2]);
    printf(" %s\n",date);
 
    double my_start_time = omp_get_wtime();

    ParU_Control Control;
    ParU_Ret info;

    //RUTAL_ALLOC_TEST(info, ParU_Analyze(A, &Sym, &Control) );
    //info = ParU_Analyze(A, &Sym, &Control);
    {                                           
        paru_set_malloc_tracking(true);         
        for (Int nmalloc = 0;; nmalloc++)       
        {                                       
            printf("#####nmalloc=%ld\n",nmalloc);
            paru_set_nmalloc(nmalloc);          
            info = ParU_Analyze(A, &Sym, &Control);
            if (info != PARU_OUT_OF_MEMORY)     
            {                                   
                printf("nmalloc=%ld\n",nmalloc);
                break;                          
            }                                   
            if (nmalloc > 100)              
            {                                   
                printf("ParU: too much malloc\n"); 
                break;                          
            }                                   
        }                                       
        printf("ParU: test failure\n");         
        paru_set_malloc_tracking(false);        
    }
    if (info != PARU_SUCCESS)
    {
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    printf("In: %ldx%ld nnz = %ld \n", Sym->m, Sym->n, Sym->anz);
    printf("Paru: Symbolic factorization is done!\n");
    ParU_Numeric *Num;
    info = ParU_Factorize(A, Sym, &Num, &Control);
    double my_time = omp_get_wtime() - my_start_time;
    if (info != PARU_SUCCESS)
    {
        printf("Paru: factorization was NOT successfull in %lf seconds.\n",
                my_time);
    }
    else
        printf("Paru: factorization was successfull in %lf seconds.\n",
                my_time);

    //~~~~~~~~~~~~~~~~~~~Test the results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Int m = Sym->m;
#if 1
    if (info == PARU_SUCCESS)
    {
        double *b = (double *)malloc(m * sizeof(double));
        double *xx = (double *)malloc(m * sizeof(double));
        for (Int i = 0; i < m; ++i) b[i] = i + 1;
        double my_solve_time_start = omp_get_wtime();
        info = ParU_Solve(Sym, Num, b, xx, &Control);
        double my_solve_time = omp_get_wtime() - my_solve_time_start;
        printf("Solve time is %lf seconds.\n", my_solve_time);
        my_start_time = omp_get_wtime();
        double resid, anorm;
        ParU_Residual(A, xx, b, m, resid, anorm, &Control);

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
        ParU_Residual(A, X, B, m, nrhs, resid, anorm, &Control);
        printf("mRhs Residual is |%.2lf|\n", resid == 0 ? 0 : log10(resid));

        free(B);
        free(X);
    }
#endif

    //~~~~~~~~~~~~~~~~~~~End computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //~~~~~~~~~~~~~~~~~~~Free Everything~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ParU_Freenum(&Num, &Control);
    ParU_Freesym(&Sym, &Control);

    cholmod_l_free_sparse(&A, cc);
    cholmod_l_finish(cc);
}
