//  =========================================================================  /
// =======================  paru_quicker_test.cpp  ==========================  /
// ==========================================================================  /
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*
 * @brief    for coverage test of bigger matrices
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

    if (mtype != CHOLMOD_SPARSE)
    {
        printf("Paru: input matrix must be sparse\n");
        exit(1);
    }

    /////This part is for covering the codes that cannot be covered through
    ///// factorizing 
    // covering alloc lines
    Int *t = NULL;
    t = (Int *)paru_alloc(1, sizeof(Int) * 0);
    t = (Int *)paru_alloc(Size_max, sizeof(Int));
    size_t size = 0;
    t = (Int *)paru_realloc(10, sizeof(Int) * 0, t, &size);
    t = (Int *)paru_realloc(10, sizeof(Int), t, &size);
    paru_free(10, sizeof(Int), t);
    Int *test_new = new Int[0];
    delete[] test_new;
    // covering elementList
    paru_element* elementList[] = {NULL};
    Int lac_0 = lac_el(elementList, 0);
    if (lac_0 != LONG_MAX)
        printf ("Some problem happend\n");

    //~~~~~~~~~~~~~~~~~~~Starting computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int ver[3];
    char date[128];
    ParU_Version(ver, date);
    printf("Paru %d.%d.%d", ver[0], ver[1], ver[2]);
    printf(" %s\n", date);

    ParU_Control Control;
    // puting Control lines to work
    Control.mem_chunk = 1024;
    Control.umfpack_ordering = 23;
    Control.umfpack_strategy = 23;
    Control.paru_max_threads = 4;
    Control.relaxed_amalgamation_threshold = -1;
    Control.paru_strategy = 23;
    Control.scale = -1;
    Control.panel_width = -1;
    Control.piv_toler = -1;
    Control.diag_toler = -1;
    Control.trivial = -1;
    Control.worthwhile_dgemm = -2;
    Control.worthwhile_trsm = -1;
    Control.paru_strategy = PARU_STRATEGY_SYMMETRIC;

    ParU_Ret info;

    info = ParU_Analyze(A, &Sym, &Control);
    if (info != PARU_SUCCESS)
    {
        printf("Paru: some problem detected during symbolic analysis\n");
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    Control.paru_strategy = PARU_STRATEGY_AUTO;
    info = ParU_Analyze(A, &Sym, &Control);
    if (info != PARU_SUCCESS)
    {
        printf("Paru: some problem detected during symbolic analysis\n");
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    printf("In: %ldx%ld nnz = %ld \n", Sym->m, Sym->n, Sym->anz);
    printf("Paru: Symbolic factorization is done!\n");
    ParU_Numeric *Num;

    info = ParU_Factorize(NULL, Sym, &Num, &Control);
    info = ParU_Factorize(A, NULL, &Num, &Control);
    info = ParU_Factorize(A, Sym, &Num, &Control);
    if (info != PARU_SUCCESS)
    {
        printf("Paru: factorization was NOT succssfull.\n");
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        ParU_Freesym(&Sym, &Control);
        return info;
    }
    else
        printf("Paru: factorization was successfull.\n");

    //~~~~~~~~~~~~~~~~~~~Test the results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Int m = Sym->m;
#if 1
    if (info == PARU_SUCCESS)
    {
        double *b = (double *)malloc(m * sizeof(double));
        double *xx = (double *)malloc(m * sizeof(double));
        for (Int i = 0; i < m; ++i) b[i] = i + 1;
        info = ParU_Solve(Sym, NULL, b, xx, &Control);  // coverage
        info = ParU_Solve(Sym, NULL, b, &Control);      // coverage
        info = ParU_Solve(Sym, Num, b, xx, &Control);
        if (info != PARU_SUCCESS)
        {
            free(b);
            free(xx);
            printf("Paru: Solve has a problem.\n");
            cholmod_l_free_sparse(&A, cc);
            cholmod_l_finish(cc);
            ParU_Freesym(&Sym, &Control);
            return info;
        }

        double resid, anorm;
        info =
            ParU_Residual(A, xx, NULL, m, resid, anorm, &Control);  // coverage
        info = ParU_Residual(A, xx, b, m, resid, anorm, &Control);
        if (info != PARU_SUCCESS)
        {
            free(b);
            free(xx);
            printf("Paru: Residual has a problem.\n");
            cholmod_l_free_sparse(&A, cc);
            cholmod_l_finish(cc);
            ParU_Freesym(&Sym, &Control);
            return info;
        }

        printf("Residual is |%.2lf| and anorm is %.2e and rcond is %.2e.\n",
               resid == 0 ? 0 : log10(resid), anorm, Num->rcond);

        for (Int i = 0; i < m; ++i) b[i] = i + 1;
        info = paru_backward(b, resid, anorm, NULL, Sym, Num, &Control);
        free(b);
        free(xx);
        const Int nrhs = 16;  // number of right handsides
        double *B = (double *)malloc(m * nrhs * sizeof(double));
        double *X = (double *)malloc(m * nrhs * sizeof(double));
        for (Int i = 0; i < m; ++i)
            for (Int j = 0; j < nrhs; ++j) B[j * m + i] = (double)(i + j + 1);

        info = ParU_Solve(Sym, NULL, nrhs, B, X, &Control);  // for coverage
        info = ParU_Solve(Sym, NULL, nrhs, B, &Control);     // for coverage
        info = ParU_Solve(Sym, Num, nrhs, B, X, &Control);
        if (info != PARU_SUCCESS)
        {
            printf("Paru: mRhs Solve has a problem.\n");
            free(B);
            free(X);
            cholmod_l_free_sparse(&A, cc);
            cholmod_l_finish(cc);
            ParU_Freesym(&Sym, &Control);
            return info;
        }
        // This one is just for more coverage
        info = ParU_Residual(A, X, NULL, m, nrhs, resid, anorm, &Control);

        info = ParU_Residual(A, X, B, m, nrhs, resid, anorm, &Control);
        if (info != PARU_SUCCESS)
        {
            printf("Paru: mRhs Residual has a problem.\n");
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

    //~~~~~~~~~~~~~~~~~~~Free Everything~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ParU_Freenum(&Num, &Control);
    ParU_Freesym(&Sym, &Control);

    cholmod_l_free_sparse(&A, cc);
    cholmod_l_finish(cc);
}
