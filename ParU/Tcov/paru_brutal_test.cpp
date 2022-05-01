//  =========================================================================  /
// =======================  paru_brutal_test.cpp  ===========================  /
// ==========================================================================  /
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*
 * @brief    to test all the allocations: malloc, calloc, realloc
 *              They fail one by one; very slow
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
    int ver[3];
    char date[128];
    ParU_Version(ver, date);
    printf("Paru %d.%d.%d", ver[0], ver[1], ver[2]);
    printf(" %s\n", date);

    ParU_Control Control;
    // puting Control lines to work
    Control.mem_chunk = 1;
    Control.umfpack_ordering = 23;
    Control.umfpack_strategy = 23;
    Control.paru_max_threads = -1;
    Control.relaxed_amalgamation_threshold = -1;
    Control.paru_strategy = 23;
    Control.scale = -1;
    Control.panel_width = -1;
    Control.piv_toler = -1;
    Control.diag_toler = -1;
    Control.trivial = -1;
    Control.worthwhile_dgemm = -2;
    Control.worthwhile_trsm = -1;

    ParU_Ret info;

    // info = ParU_Analyze(A, &Sym, &Control);
    BRUTAL_ALLOC_TEST(info, ParU_Analyze(A, &Sym, &Control));
    if (info == PARU_OUT_OF_MEMORY)
    {
        printf("Paru: some problem detected during symbolic analysis\n");
        cholmod_l_free_sparse(&A, cc);
        cholmod_l_finish(cc);
        return info;
    }
    if (Sym != NULL)
        printf("In: %ldx%ld nnz = %ld \n", Sym->m, Sym->n, Sym->anz);
    printf("Paru: Symbolic factorization is done!\n");
    ParU_Numeric *Num;

    // info = ParU_Factorize(A, Sym, &Num, &Control);
    BRUTAL_ALLOC_TEST(info, ParU_Factorize(A, Sym, &Num, &Control));
    if (info == PARU_OUT_OF_MEMORY)
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
#if 1
    if (info != PARU_OUT_OF_MEMORY)
    {
        Int m = 0;
        if (Sym) m = Sym->m;
        printf("Testing the resutls\n");
        double *b = (double *)malloc(m * sizeof(double));
        double *xx = (double *)malloc(m * sizeof(double));
        for (Int i = 0; i < m; ++i) b[i] = i + 1;
        // info = ParU_Solve(Sym, Num, b, xx, &Control);
        printf("Testing Solve\n");
        BRUTAL_ALLOC_TEST(info, ParU_Solve(Sym, Num, b, xx, &Control));
        double resid, anorm;
        // info = ParU_Residual(A, xx, b, m, resid, anorm, &Control);
        printf("Testing Residual\n");
        BRUTAL_ALLOC_TEST(info,
                          ParU_Residual(A, xx, b, m, resid, anorm, &Control));
        for (Int i = 0; i < m; ++i) b[i] = i + 1;
        // info = paru_backward(b, resid, anorm, A, Sym, Num, &Control);
        printf("Testing backward\n");
        BRUTAL_ALLOC_TEST(
            info, paru_backward(b, resid, anorm, A, Sym, Num, &Control));
        free(b);
        free(xx);
        const Int nrhs = 16;  // number of right handsides
        double *B = (double *)malloc(m * nrhs * sizeof(double));
        double *X = (double *)malloc(m * nrhs * sizeof(double));
        for (Int i = 0; i < m; ++i)
            for (Int j = 0; j < nrhs; ++j) B[j * m + i] = (double)(i + j + 1);

        // info = ParU_Solve(Sym, Num, nrhs, B, X, &Control);
        printf("Testing mRHS Solve\n");
        BRUTAL_ALLOC_TEST(info, ParU_Solve(Sym, Num, nrhs, B, X, &Control));
            // info = ParU_Residual(A, X, B, m, nrhs, resid, anorm, &Control);
        printf("Testing mRHS Residual\n");
        BRUTAL_ALLOC_TEST(
            info, ParU_Residual(A, X, B, m, nrhs, resid, anorm, &Control));
            printf("mRhs Residual is |%.2lf|\n", resid == 0 ? 0 : log10(resid));

        free(B);
        free(X);
    }
#endif

    //~~~~~~~~~~~~~~~~~~~End computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //~~~~~~~~~~~~~~~~~~~Free Everything~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    printf("freeing stuff\n");
    ParU_Freenum(&Num, &Control);
    ParU_Freesym(&Sym, &Control);

    cholmod_l_free_sparse(&A, cc);
    cholmod_l_finish(cc);
}
