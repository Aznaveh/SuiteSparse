/** =========================================================================  /
 * =======================  paru_factorize  =================================  /
 * ========================================================================== */
/*! @brief Doing the BLAS factorization in different panels and call degree
 * update it it is necessary
 *
 *
 * @author Aznaveh
 */

#include "paru_internal.hpp"
#define TOLER 0.1  // pivot tolerance

extern "C" void dgetrf_(BLAS_INT *dim1, BLAS_INT *dim2, double *a,
                        BLAS_INT *lda, BLAS_INT *ipiv, BLAS_INT *info);
// dger is already defined in ~/SuiteSparse/CHOLMOD/Include/cholmod_blas.h
void inline swap_int(Int *a, Int *b)
{
    Int tmp = *a;
    *a = *b;
    *b = tmp;
}

template <class T>
void inline swap(T &a, T &b)
{
    T c(a);
    a = b;
    b = c;
}

void swap_rows(double *F, Int *frowList, Int m, Int n, Int r1, Int r2)
{
    // This function also swap rows r1 and r2 wholly and indices
    if (r1 == r2) return;
    swap(frowList[r1], frowList[r2]);
    for (Int j = 0; j < n; j++)
        // each column
        swap(F[j * m + r1], F[j * m + r2]);
}

Int paru_panel_factorize(Int f, Int m, Int n, const Int panel_width,
                         Int panel_num, Int row_end, paru_matrix *paruMatInfo)
{
    // works like dgetf2f.f in netlib v3.0  here is a link:
    // https://github.com/xianyi/OpenBLAS/blob/develop/reference/dgetf2f.f
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% Inside panel factorization %ld \n", panel_num));

    Int *row_degree_bound = paruMatInfo->row_degree_bound;
    Int j1 = panel_num * panel_width;  // panel starting column

    //  j1 <= panel columns < j2
    //     last panel might be smaller
    Int j2 = (j1 + panel_width < n) ? j1 + panel_width : n;

    PRLEVEL(1, ("%% j1= %ld j2 =%ld \n", j1, j2));
    PRLEVEL(1, ("%% row_end= %ld\n", row_end));

    // ASSERT(row_end >= j2);

    Int *frowList = paruMatInfo->frowList[f];
    paru_fac *LUs = paruMatInfo->partial_LUs;
    double *F = LUs[f].p;

#ifndef NDEBUG  // Printing the panel
    Int num_col_panel = j2 - j1;
    Int p = 1;
    PRLEVEL(p, ("%% Starting the factorization\n"));
    PRLEVEL(p, ("%% This Panel:\n"));
    for (Int r = j1; r < row_end; r++)
    {
        PRLEVEL(p, ("%% %ld\t", frowList[r]));
        for (Int c = j1; c < j2; c++) PRLEVEL(p, (" %2.5lf\t", F[c * m + r]));
        PRLEVEL(p, ("\n"));
    }
#endif
    Int *Super = paruMatInfo->LUsym->Super;
    Int col1 = Super[f]; /* fornt F has columns col1:col2-1 */
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *Qfill = LUsym->Qfill;
    Int *Pinit = LUsym->Pinit;

    // column jth of the panel
    for (Int j = j1; j < j2; j++)
    {
        // for fat fronts
        if (j >= row_end) break;

        PRLEVEL(1, ("%% j = %ld\n", j));

        // Initializing maximum element in the column
        Int row_max = j;
#ifndef NDEBUG
        // FIXME a problem here
        Int row_deg_max = row_degree_bound[frowList[row_max]];
#endif
        double maxval = F[j * m + row_max];
        PRLEVEL(1, ("%% before search max value= %2.4lf row_deg = %ld\n",
                    maxval, row_deg_max));

        Int origCol = Qfill ? Qfill[j + col1] : j + col1;
        Int row_diag = (origCol == Pinit[frowList[j]]) ? j : -1;
        double diag_val = maxval;  // initialization
        PRLEVEL(1, ("%%curCol=%ld origCol= %ld row_diag=%ld\n", j + col1,
                    origCol, row_diag));

        PRLEVEL(1, ("%%##j=%ld value= %2.4lf\n", j, F[j * m + j]));
        // find max
        for (Int i = j + 1; i < row_end; i++)
        {
            PRLEVEL(1, ("%%i=%ld value= %2.4lf", i, F[j * m + i]));
            PRLEVEL(1, (" deg = %ld \n", row_degree_bound[frowList[i]]));
            if (fabs(maxval) < fabs(F[j * m + i]))
            {
                row_max = i;
                maxval = F[j * m + i];
            }
            Int origRow = Pinit[frowList[i]];
            PRLEVEL(1, ("%%curRow=%ld origRow= %ld\n", frowList[i], origRow));
            if (origRow == origCol)
            {
                PRLEVEL(1, ("%%Found it %2.4lf\n", F[j * m + i]));
                row_diag = i;
                diag_val = F[j * m + i];
            }
        }

#ifndef NDEBUG
        row_deg_max = row_degree_bound[frowList[row_max]];
#endif

        PRLEVEL(1, ("%% max value= %2.4lf\n", maxval));

        if (maxval == 0)
        {
            PRLEVEL(-1, ("%% NO pivot found in %ld\n", f));
            continue;
        }

        // initialzing pivot as max numeric value
        double piv = maxval;
        Int row_piv = row_max;
        Int chose_diag = 0;

        if (LUsym->strategy == UMFPACK_STRATEGY_SYMMETRIC)
        {
            if (row_diag != -1)
                if (fabs(TOLER * TOLER * maxval) < fabs(diag_val))
                {
                    piv = diag_val;
                    row_piv = row_diag;
                    PRLEVEL(1, ("%% symmetric pivot piv value= %2.4lf"
                                " row_piv=%ld\n",
                                piv, row_piv));
                    chose_diag = 1;
                }
        }

        // find sparsest between accepteble ones
        // if not symmetric or the diagonal is not good enough
        Int row_deg_sp = row_degree_bound[frowList[row_max]];
        if (chose_diag == 0)
        {
            Int row_sp = row_max;
            for (Int i = j; i < row_end; i++)
                if (fabs(TOLER * maxval) < fabs(F[j * m + i]) &&
                    row_degree_bound[frowList[i]] < row_deg_sp)
                {  // numerically acceptalbe and sparser
                    piv = F[j * m + i];
                    row_deg_sp = row_degree_bound[frowList[i]];
                    row_sp = i;
                }
            row_piv = row_sp;
        }

        PRLEVEL(1, ("%% piv value= %2.4lf row_deg=%ld\n", piv, row_deg_sp));

        // swap rows
        PRLEVEL(1, ("%% Swaping rows j=%ld, row_piv=%ld\n", j, row_piv));
        swap_rows(F, frowList, m, n, j, row_piv);

#ifndef NDEBUG  // Printing the pivotal front
        p = 1;
        if (row_piv != row_max) PRLEVEL(p, ("%% \n"));
        PRLEVEL(p, ("%% After Swaping\n"));
        PRLEVEL(p, (" \n"));
        for (Int r = 0; r < row_end; r++)
        {
            PRLEVEL(p, ("%% %ld\t", frowList[r]));
            for (Int c = 0; c < num_col_panel; c++)
                PRLEVEL(p, (" %2.5lf\t", F[c * m + r]));
            PRLEVEL(p, ("\n"));
        }
#endif

        // dscal //TODO?: loop unroll is also possible

        if (j < row_end - 1)
        {
            PRLEVEL(1, ("%% dscal\n"));
            for (Int i = j + 1; i < row_end; i++)
            {
                PRLEVEL(1, ("%%i=%ld value= %2.4lf", i, F[j * m + i]));
                F[j * m + i] /= piv;
                PRLEVEL(1, (" -> %2.4lf\n", F[j * m + i]));
            }
        }

        // dger
        /*               dgemm   A := alpha *x*y**T + A
         *
         *
         *                <----------fp------------------->
         *                        j1 current j2
         *                         ^  panel  ^
         *                         |         |
         *                         |   j     |
         *             F           |   ^     |
         *              \  ____..._|___|_____|__________...___
         * ^              |\      |         |                 |
         * |              | \     |<--panel | rest of         |
         * |              |  \    | width-> |  piv front      |
         * |              |___\...|_______ _|_________ ... ___|
         * |   ^    j1--> |       |\* *|    |                 |
         * | panel        |       |**\*|    |                 |
         * | width        |       |___|Pyyyy|                 |
         * |   v          |____...|___|xAAAA__________...____ |
         * |        j2--> |       |   |xAAAA|                 |
         * rowCount       |       |   |xAAAA|                 |
         * |              |       |   |xAAAA|                 |
         * |              |       |   |xAAAA|                 |
         * |              |       |row_end  |                 |
         * |              .       .         .                 .
         * |              .       .         .                 .
         * |              .       .         .                 .
         * v              |___....____________________..._____|
         *
         */

        if (j < j2 - 1)
        {
            BLAS_INT M = (BLAS_INT)row_end - 1 - j;
            BLAS_INT N = (BLAS_INT)j2 - 1 - j;
            double alpha = -1.0;
            double *X = F + j * m + j + 1;
            BLAS_INT Incx = (BLAS_INT)1;
            double *Y = F + j * m + j + m;
            BLAS_INT Incy = (BLAS_INT)m;
            double *A = F + j * m + j + m + 1;
            BLAS_INT lda = (BLAS_INT)m;

#ifndef NDEBUG  // Printing dger input
            Int p = 1;
            PRLEVEL(p, ("%% lda =%d ", lda));
            PRLEVEL(p, ("%% M =%d ", M));
            PRLEVEL(p, ("N =%d \n %%", N));
            PRLEVEL(p, ("%% x= (%d)", N));
            for (Int i = 0; i < M; i++) PRLEVEL(p, (" %lf ", X[i]));
            PRLEVEL(p, ("\n %% y= (%d)", N));
            for (Int j = 0; j < N; j++) PRLEVEL(p, (" %lf ", Y[j * m]));
            PRLEVEL(p, ("\n"));

#endif
            BLAS_DGER(&M, &N, &alpha, X, &Incx, Y, &Incy, A, &lda);
#ifdef COUNT_FLOPS
            paruMatInfo->flp_cnt_dger += (double)2 * M * N;
#ifndef NDEBUG
            PRLEVEL(p, ("\n%% FlopCount Dger fac %d %d ", M, N));
            PRLEVEL(p, ("cnt = %lf\n ", paruMatInfo->flp_cnt_dger));
#endif
#endif
        }

#ifndef NDEBUG  // Printing the pivotal front
        Int p = 1;
        PRLEVEL(p, ("%% After dger\n"));
        for (Int r = j; r < row_end; r++)
        {
            PRLEVEL(p, ("%% %ld\t", frowList[r]));
            for (Int c = j; c < j2; c++)
                PRLEVEL(p, (" %2.5lf\t", F[c * m + r]));
            PRLEVEL(p, ("\n"));
        }
#endif
    }
    return 1;
}

Int paru_factorize(Int f, Int start_fac, std::vector<Int> &panel_row,
                   std::set<Int> &stl_colSet,
                   std::vector<Int> &pivotal_elements, paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);

    Int *Super = paruMatInfo->LUsym->Super;
    Int col1 = Super[f]; /* fornt F has columns col1:col2-1 */
    Int col2 = Super[f + 1];
    Int fp = col2 - col1; /* first fp columns are pivotal */

    paru_fac *LUs = paruMatInfo->partial_LUs;
    Int rowCount = paruMatInfo->frowCount[f];
    double *F = LUs[f].p;

    Int panel_width = paruMatInfo->panel_width;
    for (Int panel_num = 0;; panel_num++)
    {
#ifndef NDEBUG  // Printing the pivotal front
        Int *frowList = paruMatInfo->frowList[f];
        Int p = 1;
        PRLEVEL(p, ("%%Pivotal Front Before %ld\n", panel_num));

        for (Int r = 0; r < rowCount; r++)
        {
            PRLEVEL(p, ("%% %ld\t", frowList[r]));
            for (Int c = 0; c < fp; c++)
            {
                PRLEVEL(p, (" %2.5lf\t", F[c * rowCount + r]));
            }
            PRLEVEL(p, ("\n"));
        }
#endif

        Int row_end = panel_row[panel_num];
        Int j1 = panel_num * panel_width;
        Int j2 = (panel_num + 1) * panel_width;
        // factorize current panel
        paru_panel_factorize(f, rowCount, fp, panel_width, panel_num, row_end,
                             paruMatInfo);

        // This can be done parallel to the  next part
        if (paruMatInfo->LUsym->Cm[f] != 0)  // if there is potential column
                                             // left
            paru_update_rowDeg(panel_num, row_end, f, start_fac, stl_colSet,
                               pivotal_elements, paruMatInfo);

        if (j2 >= fp)  // if it is the last panel
            break;

        /*               trsm
         *
         *        F = fully summed part of the pivotal front
         *           op( A ) * B = alpha*B
         *
         *                <----------fp------------------->
         *                        j1 current j2
         *    F                    ^  panel  ^
         *     \           ____..._|_________|__________...___
         * ^              |\      |         |                 |
         * |              | \     |<--panel | rest of         |
         * |              |  \    | width-> |  piv front      |
         * |              |___\...|_______ _|_________ ... ___|
         * |   ^    j1--> |       |\        |                 |
         * | panel        |       |**\ A    |   B(In out)     |
         * | width        |       |*L**\    |                 |
         * |   v          |____...|******\ _|_________...____ |
         * |        j2--> |       |         |                 |
         * rowCount       |       |         |                 |
         * |              .       .         .                 .
         * |              .       .         .                 .
         * |              .       .         .                 .
         * v              |___....____________________..._____|
         *
         */
        ASSERT(j2 < fp);
        {
            BLAS_INT M = (BLAS_INT)panel_width;
            BLAS_INT N = (BLAS_INT)fp - j2;
            double alpha = 1.0;
            double *A = F + j1 * rowCount + j1;
            BLAS_INT lda = (BLAS_INT)rowCount;
            double *B = F + j2 * rowCount + j1;
            BLAS_INT ldb = (BLAS_INT)rowCount;
#ifndef NDEBUG
            Int p = 1;
            PRLEVEL(p, ("%% M =%d N = %d alpha = %f \n", M, N, alpha));
            PRLEVEL(p, ("%% lda =%d ldb =%d\n", lda, ldb));
            PRLEVEL(
                p, ("%% Pivotal Front Before Trsm: %ld x %ld\n", fp, rowCount));
            for (Int r = 0; r < rowCount; r++)
            {
                PRLEVEL(p, ("%% %ld\t", frowList[r]));
                for (Int c = 0; c < fp; c++)
                    PRLEVEL(p, (" %2.5lf\t", F[c * rowCount + r]));
                PRLEVEL(p, ("\n"));
            }

#endif
            BLAS_DTRSM("L", "L", "N", "U", &M, &N, &alpha, A, &lda, B, &ldb);
#ifdef COUNT_FLOPS
            paruMatInfo->flp_cnt_trsm += (double)(M + 1) * M * N;
#ifndef NDEBUG
            p = 0;
            PRLEVEL(p, ("\n%% FlopCount Trsm factorize %d %d ", M, N));
            PRLEVEL(p, ("cnt = %lf\n ", paruMatInfo->flp_cnt_trsm));
#endif

#endif

#ifndef NDEBUG
            PRLEVEL(p, ("%% Pivotal Front After Trsm: %ld x %ld\n %%", fp,
                        rowCount));
            for (Int r = 0; r < rowCount; r++)
            {
                PRLEVEL(p, ("%% %ld\t", frowList[r]));
                for (Int c = 0; c < fp; c++)
                    PRLEVEL(p, (" %2.5lf\t", F[c * rowCount + r]));
                PRLEVEL(p, ("\n"));
            }
#endif
        }

        /*               dgemm   C := alpha*op(A)*op(B) + beta*C
         *
         *        F = fully summed part of the pivotal front
         *
         *                <----------fp------------------->
         *                        j1 current j2
         *    F                    ^  panel  ^
         *     \           ____..._|_________|__________...___
         * ^              |\      |         |                 |
         * |              | \     |<--panel | rest of         |
         * |              |  \    | width-> |  piv front      |
         * |              |___\...|_______ _|_________ ... ___|
         * |   ^    j1--> |       |\        |**************** |
         * | panel        |       |  \      |****In***B****** |
         * | width        |       |    \    |**************** |
         * |   v          |____...|_______\___________...____ |
         * |        j2--> |       |******** |ccccccccccccccccc|
         * rowCount       |       |******** |ccccccccccccccccc|
         * |              .       .******** .ccccccCcccccccccc.
         * |              .       .***In*** .ccccOutcccccccccc.
         * |              .       .***A**** .ccccccccccccccccc.
         * |              |       |******** |ccccccccccccccccc|
         * |              |       |******** |ccccccccccccccccc|
         * |              |       |row_end  |                 |
         * |              |       |         |                 |
         * v              |___....|_______ _|__________...____|
         *
         */

        {
            BLAS_INT M = (BLAS_INT)(row_end - j2);

            BLAS_INT N = (BLAS_INT)fp - j2;
            BLAS_INT K = (BLAS_INT)panel_width;
            double alpha = -1;
            double *A = F + j1 * rowCount + j2;
            BLAS_INT lda = (BLAS_INT)rowCount;
            double *B = F + j2 * rowCount + j1;
            BLAS_INT ldb = (BLAS_INT)rowCount;
            double beta = 1;  // keep current values
            double *C = F + j2 * rowCount + j2;
            BLAS_INT ldc = (BLAS_INT)rowCount;
#ifndef NDEBUG
            Int p = 1;
            PRLEVEL(p, ("%% DGEMM "));
            PRLEVEL(p,
                    ("%% M =%d K = %d N = %d alpha = %f \n", M, K, N, alpha));
            PRLEVEL(p, ("%% lda =%d ldb =%d\n", lda, ldb));
            PRLEVEL(p, ("%% j2 =%ld j1=%ld\n", j2, j1));
            PRLEVEL(p, ("\n %%"));
#endif

            // double start_time = omp_get_wtime();
            BLAS_DGEMM("N", "N", &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C,
                       &ldc);
            // double tot_time = omp_get_wtime() - start_time;
            // printf ("%ld  %lf ",f, tot_time);
            // printf ("%d %d %d ",M ,N, K);
            // printf ("%d %d %d\n ",lda ,ldb, ldc);
#ifdef COUNT_FLOPS
            paruMatInfo->flp_cnt_dgemm += (double)2 * M * N * K;
            paruMatInfo->flp_cnt_real_dgemm += (double)2 * M * N * K;
#ifndef NDEBUG
            PRLEVEL(p, ("\n%% FlopCount Dgemm factorize %d %d %d ", M, N, K));
            PRLEVEL(p, ("%d %d %d \n", M, N, K));

            PRLEVEL(p, ("cnt = %lf\n ", paruMatInfo->flp_cnt_dgemm));
#endif
#endif
        }

#ifndef NDEBUG
        PRLEVEL(p,
                ("%% Pivotal Front After Dgemm: %ld x %ld\n %%", fp, rowCount));
        for (Int r = 0; r < rowCount; r++)
        {
            PRLEVEL(p, ("%% %ld\t", frowList[r]));
            for (Int c = 0; c < fp; c++)
                PRLEVEL(p, (" %2.5lf\t", F[c * rowCount + r]));
            PRLEVEL(p, ("\n"));
        }
#endif
    }
    return 0;
}
