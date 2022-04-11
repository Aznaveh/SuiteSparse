////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_fs_factorize  /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*! @brief Doing the BLAS factorization in different panels and call degre////e
 * update when it is necessary.
 *
 * @author Aznaveh
 */

#include "paru_internal.hpp"

void swap_rows(double *F, Int *frowList, Int m, Int n, Int r1, Int r2,
               ParU_Numeric *Num)
{
    // This function also swap rows r1 and r2 wholly and indices
    if (r1 == r2) return;
    std::swap(frowList[r1], frowList[r2]);
    // So dissappointed in parallelism this part; SLOWDOWNN

    // Int naft; //number of active frontal tasks
    // pragma omp atomic read
    // naft = Num->naft;
    // const Int max_threads = Control->paru_max_threads;
    // if ( (naft == 1) && (n > 1024) )
    // printf ("naft=%ld, max_threads=%ld num_tasks=%ld n =%ld \n",
    //        naft, max_threads, max_threads/(naft), n);
    // pragma omp parallel if ( (naft == 1) && (n > 1024) )
    // pragma omp single
    // pragma omp taskloop num_tasks(max_threads/(naft+1))

    for (Int j = 0; j < n; j++)
        // each column
        std::swap(F[j * m + r1], F[j * m + r2]);
}

Int paru_panel_factorize(Int f, Int m, Int n, const Int panel_width, 
        Int panel_num, Int row_end, paru_work *Work, ParU_Numeric *Num)
{
    // works like dgetf2f.f in netlib v3.0  here is a link:
    // https://github.com/xianyi/OpenBLAS/blob/develop/reference/dgetf2f.f
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;
    PRLEVEL(1, ("%% Inside panel factorization %ld \n", panel_num));

    Int *row_degree_bound = Work->row_degree_bound;
    ParU_Control *Control = Num->Control;
    Int j1 = panel_num * panel_width;  // panel starting column

    //  j1 <= panel columns < j2
    //     last panel might be smaller
    Int j2 = (j1 + panel_width < n) ? j1 + panel_width : n;

    PRLEVEL(1, ("%% j1= %ld j2 =%ld \n", j1, j2));
    PRLEVEL(1, ("%% row_end= %ld\n", row_end));

    // ASSERT(row_end >= j2);

    Int *frowList = Num->frowList[f];
    ParU_Factors *LUs = Num->partial_LUs;
    double *F = LUs[f].p;

#ifndef NDEBUG  // Printing the panel
    Int num_col_panel = j2 - j1;
    PRLEVEL(PR, ("%% Starting the factorization\n"));
    PRLEVEL(PR, ("%% This Panel:\n"));
    for (Int r = j1; r < row_end; r++)
    {
        PRLEVEL(PR, ("%% %ld\t", frowList[r]));
        for (Int c = j1; c < j2; c++) PRLEVEL(PR, (" %2.5lf\t", F[c * m + r]));
        PRLEVEL(PR, ("\n"));
    }
#endif
    ParU_Symbolic *Sym = Work->Sym;
    Int *Super = Sym->Super;
    Int col1 = Super[f]; /* fornt F has columns col1:col2-1 */
    Int *Diag_map = Work->Diag_map;
    Int n1 = Sym->n1;

    // column jth of the panel
    for (Int j = j1; j < j2; j++)
    {
        // for fat fronts
        if (j >= row_end) break;

        PRLEVEL(1, ("%% j = %ld\n", j));

        // Initializing maximum element in the column
        Int row_max = j;

        double maxval = F[j * m + row_max];

#ifndef NDEBUG
        Int row_deg_max = row_degree_bound[frowList[row_max]];
        PRLEVEL(1, ("%% before search max value= %2.4lf row_deg = %ld\n",
                    maxval, row_deg_max));
#endif

        Int row_diag = (Diag_map) ? Diag_map[col1 + j + n1] - n1 : -1;
        double diag_val = maxval;  // initialization
        Int diag_found = frowList[j] == row_diag ? j : -1;
        PRLEVEL(1, ("%%curCol=%ld row_diag=%ld\n", j + col1 + n1, row_diag));
        PRLEVEL(1, ("%%##j=%ld value= %2.4lf\n", j, F[j * m + j]));

        for (Int i = j + 1; i < row_end; i++)
        {  // find max
            PRLEVEL(1, ("%%i=%ld value= %2.4lf", i, F[j * m + i]));
            PRLEVEL(1, (" deg = %ld \n", row_degree_bound[frowList[i]]));
            if (fabs(maxval) < fabs(F[j * m + i]))
            {
                row_max = i;
                maxval = F[j * m + i];
            }
            if (frowList[i] == row_diag)  // find diag
            {
                PRLEVEL(1, ("%%Found it %2.4lf\n", F[j * m + i]));
                // row_diag = i;
                diag_found = i;
                diag_val = F[j * m + i];
            }
        }

        /*** ATTENTION: IT FACES A NUMERICAL PROBLEM EVEN IF I USE REDUCTION**/
        /*  pragma omp declare reduction
        //  (maxfabs : double:
        //   omp_out = fabs(omp_in) > fabs(omp_out) ? omp_in : omp_out )
        */

        /* pragma omp taskloop default(none)
        //shared(maxval, F,row_max, row_end, j,
        //        row_diag, diag_val, diag_found, m, frowList) grainsize(512)
        //for (Int i = j + 1; i < row_end; i++)
        //{  // find max
        //    double value = F[j * m + i];
        //    if (fabs(maxval) < fabs(value))
        //    {
        //        pragma omp critical
        //        {
        //            maxval = value;
        //            row_max = i;
        //        }
        //    }
        //    if (frowList[i] == row_diag)  // find diag
        //    {
        //        pragma omp critical
        //        {//actually this will be seen at most one time
        //            diag_found = i;
        //            diag_val = value;
        //        }
        //    }
        //}
        */

#ifndef NDEBUG
        row_deg_max = row_degree_bound[frowList[row_max]];
#endif

        PRLEVEL(1, ("%% max value= %2.4lf\n", maxval));

        if (maxval == 0)
        {
            PRLEVEL(1, ("%% NO pivot found in %ld\n", n1 + col1 + j));
            #pragma omp atomic write
            Num->res = PARU_SINGULAR;
            continue;
        }
        //XXX if piv valuse is less than epsilon the matrix might be singular
        //double eps = 1e-15;
        //if (maxval < eps)
        //{
        //    #pragma omp atomic write
        //    Num->res = PARU_SINGULAR;
        //}

        // initialzing pivot as max numeric value
        double piv = maxval;
        Int row_piv = row_max;
        Int chose_diag = 0;

        if (Control->paru_strategy == PARU_STRATEGY_SYMMETRIC)
        {
            if (diag_found != -1)
            {
                if (fabs(Control->diag_toler * maxval) < fabs(diag_val))
                {
                    piv = diag_val;
                    row_piv = diag_found;
                    PRLEVEL(1, ("%% symmetric pivot piv value= %2.4lf"
                                " row_piv=%ld\n",
                                piv, row_piv));
                    chose_diag = 1;
                }
#ifndef NDEBUG
                else
                {
                    PRLEVEL(1, ("%% diag found but too small %ld"
                                " maxval=%2.4lf diag_val=%e \n",
                                row_piv, maxval, diag_val));
                }
#endif
            }
#ifndef NDEBUG
            else
            {
                PRLEVEL(1, ("%% diag not found %ld\n", row_piv));
            }
#endif
        }

        // find sparsest between accepteble ones
        // if not symmetric or the diagonal is not good enough
        Int row_deg_sp = row_degree_bound[frowList[row_max]];
        if (chose_diag == 0)
        {
            Int row_sp = row_max;
            // pragma omp taskloop  default(none) shared(maxval, F, row_sp, j,
            // row_end, m, piv, frowList, row_degree_bound, row_deg_sp)
            // grainsize(512)

            for (Int i = j; i < row_end; i++)
            {
                double value = F[j * m + i];
                if (fabs(Control->piv_toler * maxval) < fabs(value) &&
                    row_degree_bound[frowList[i]] < row_deg_sp)
                {  // numerically acceptalbe and sparser
                    // pragma omp critical
                    {
                        piv = value;
                        row_deg_sp = row_degree_bound[frowList[i]];
                        row_sp = i;
                    }
                }
            }
            row_piv = row_sp;
        }

        if (Control->paru_strategy == PARU_STRATEGY_SYMMETRIC 
                && chose_diag == 0)
        {
            Int pivcol = col1 + j + n1;      // S col index + n1
            Int pivrow = frowList[row_piv];  // S row index
            paru_Diag_update(pivcol, pivrow, Work);
            PRLEVEL(1, ("%% symmetric matrix but the diag didn't picked for "
                        "row_piv=%ld\n",
                        row_piv));
        }
        PRLEVEL(1, ("%% piv value= %2.4lf row_deg=%ld\n", piv, row_deg_sp));
        PRLEVEL(-1, ("%% piv value= %e \n", piv));
        //XXX if piv valuse is less than epsilon the matrix might be singular
        //double eps = 1e-15;
        //if (piv < eps)
        //{
        //    #pragma omp atomic write
        //    Num->res = PARU_SINGULAR;
        //}

        // swap rows
        PRLEVEL(1, ("%% Swaping rows j=%ld, row_piv=%ld\n", j, row_piv));
        swap_rows(F, frowList, m, n, j, row_piv, Num);

#ifndef NDEBUG  // Printing the pivotal front
        PR = 1;
        if (row_piv != row_max) PRLEVEL(PR, ("%% \n"));
        PRLEVEL(PR, ("%% After Swaping\n"));
        PRLEVEL(PR, (" \n"));
        for (Int r = 0; r < row_end; r++)
        {
            PRLEVEL(PR, ("%% %ld\t", frowList[r]));
            for (Int c = 0; c < num_col_panel; c++)
                PRLEVEL(PR, (" %2.5lf\t", F[c * m + r]));
            PRLEVEL(PR, ("\n"));
        }
#endif

        // dscal loop unroll is also possible

        if (j < row_end - 1)
        {
            PRLEVEL(1, ("%% dscal\n"));
            // pragma omp taskloop simd default(none)
            // shared(j, row_end, F, m, piv) if(row_end-j > 1024)
            #pragma omp simd
            for (Int i = j + 1; i < row_end; i++)
            {
                // printf("%%i=%ld value= %2.4lf", i, F[j * m + i]);
                F[j * m + i] /= piv;
                // printf(" -> %2.4lf\n", F[j * m + i]);
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
            Int PR = 1;
            PRLEVEL(PR, ("%% lda =%d ", lda));
            PRLEVEL(PR, ("%% M =%d ", M));
            PRLEVEL(PR, ("N =%d \n %%", N));
            PRLEVEL(PR, ("%% x= (%d)", N));
            for (Int i = 0; i < M; i++) PRLEVEL(PR, (" %lf ", X[i]));
            PRLEVEL(PR, ("\n %% y= (%d)", N));
            for (Int j = 0; j < N; j++) PRLEVEL(PR, (" %lf ", Y[j * m]));
            PRLEVEL(PR, ("\n"));

#endif
            // BLAS_DGER(&M, &N, &alpha, X, &Incx, Y, &Incy, A, &lda);
            cblas_dger(CblasColMajor, M, N, alpha, X, Incx, Y, Incy, A, lda);
#ifdef COUNT_FLOPS
// printf("dger adding to flop count %ld\n", M*N*2);
            #pragma omp atomic update
            Num->flp_cnt_dger += (double)2 * M * N;
#ifndef NDEBUG
            PRLEVEL(PR, ("\n%% FlopCount Dger fac %d %d ", M, N));
            PRLEVEL(PR, ("cnt = %lf\n ", Num->flp_cnt_dger));
#endif
#endif
        }

#ifndef NDEBUG  // Printing the pivotal front
        Int PR = 1;
        PRLEVEL(PR, ("%% After dger\n"));
        for (Int r = j; r < row_end; r++)
        {
            PRLEVEL(PR, ("%% %ld\t", frowList[r]));
            for (Int c = j; c < j2; c++)
                PRLEVEL(PR, (" %2.5lf\t", F[c * m + r]));
            PRLEVEL(PR, ("\n"));
        }
#endif
    }
    return 1;
}

Int paru_factorize_full_summed(Int f, Int start_fac,
                               std::vector<Int> &panel_row,
                               std::set<Int> &stl_colSet,
                               std::vector<Int> &pivotal_elements,
                               paru_work *Work, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;

    Int *Super = Work->Sym->Super;
    Int col1 = Super[f]; /* fornt F has columns col1:col2-1 */
    Int col2 = Super[f + 1];
    Int fp = col2 - col1; /* first fp columns are pivotal */

    ParU_Factors *LUs = Num->partial_LUs;
    Int rowCount = Num->frowCount[f];
    double *F = LUs[f].p;

    ParU_Control *Control = Num->Control;
    Int panel_width = Control->panel_width;
 
    Int num_panels =
        (fp % panel_width == 0) ? fp / panel_width : fp / panel_width + 1;
    for (Int panel_num = 0; panel_num < num_panels; panel_num++)
    {
#ifndef NDEBUG  // Printing the pivotal front
        Int *frowList = Num->frowList[f];
        PRLEVEL(PR, ("%%Pivotal Front Before %ld\n", panel_num));

        for (Int r = 0; r < rowCount; r++)
        {
            PRLEVEL(PR, ("%% %ld\t", frowList[r]));
            for (Int c = 0; c < fp; c++)
            {
                PRLEVEL(PR, (" %2.5lf\t", F[c * rowCount + r]));
            }
            PRLEVEL(PR, ("\n"));
        }
#endif

        Int row_end = panel_row[panel_num];
        Int j1 = panel_num * panel_width;
        Int j2 = (panel_num + 1) * panel_width;
        // factorize current panel
        paru_panel_factorize(f, rowCount, fp, panel_width, panel_num, row_end,
                             Work, Num);
        // Int naft; //number of active frontal tasks
        // pragma omp atomic read
        // naft = Work->naft;
        // pragma omp parallel  proc_bind(close) if(naft == 1)
        // pragma omp single
        {
            // update row degree and dgeem can be done in parallel
            // pragma omp task default(none) mergeable
            // shared(Num, pivotal_elements, stl_colSet)
            // shared(panel_num, row_end, f, start_fac)

            if (Work->Sym->Cm[f] != 0)
            {  // if there is potential column left
                paru_update_rowDeg(panel_num, row_end, f, start_fac, stl_colSet,
                                   pivotal_elements, Work, Num);
            }

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
            // pragma omp task  shared(F)
            // shared(panel_width, j1, j2, fp, f, rowCount)
            if (j2 < fp)  // if it is not the last
            {
                BLAS_INT M = (BLAS_INT)panel_width;
                BLAS_INT N = (BLAS_INT)fp - j2;
                double alpha = 1.0;
                double *A = F + j1 * rowCount + j1;
                BLAS_INT lda = (BLAS_INT)rowCount;
                double *B = F + j2 * rowCount + j1;
                BLAS_INT ldb = (BLAS_INT)rowCount;
#ifndef NDEBUG
                Int PR = 1;
                PRLEVEL(PR, ("%% M =%d N = %d alpha = %f \n", M, N, alpha));
                PRLEVEL(PR, ("%% lda =%d ldb =%d\n", lda, ldb));
                PRLEVEL(PR, ("%% Pivotal Front Before Trsm: %ld x %ld\n", fp,
                             rowCount));
                for (Int r = 0; r < rowCount; r++)
                {
                    PRLEVEL(PR, ("%% %ld\t", frowList[r]));
                    for (Int c = 0; c < fp; c++)
                        PRLEVEL(PR, (" %2.5lf\t", F[c * rowCount + r]));
                    PRLEVEL(PR, ("\n"));
                }

#endif
                paru_tasked_trsm(f, M, N, alpha, A, lda, B, ldb, Work, Num);
#ifndef NDEBUG
                PRLEVEL(PR, ("%% Pivotal Front After Trsm: %ld x %ld\n %%", fp,
                             rowCount));
                for (Int r = 0; r < rowCount; r++)
                {
                    PRLEVEL(PR, ("%% %ld\t", frowList[r]));
                    for (Int c = 0; c < fp; c++)
                        PRLEVEL(PR, (" %2.5lf\t", F[c * rowCount + r]));
                    PRLEVEL(PR, ("\n"));
                }
#endif
            }
        }  // end of parallel region; it doesn't show good performance

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

        if (j2 < fp)
        {
            BLAS_INT M = (BLAS_INT)(row_end - j2);

            BLAS_INT N = (BLAS_INT)fp - j2;
            BLAS_INT K = (BLAS_INT)panel_width;
            // alpha = -1;
            double *A = F + j1 * rowCount + j2;
            BLAS_INT lda = (BLAS_INT)rowCount;
            double *B = F + j2 * rowCount + j1;
            BLAS_INT ldb = (BLAS_INT)rowCount;
            // double beta = 1;  // keep current values
            double *C = F + j2 * rowCount + j2;
            BLAS_INT ldc = (BLAS_INT)rowCount;
#ifndef NDEBUG
            Int PR = 1;
            PRLEVEL(PR, ("%% DGEMM "));
            PRLEVEL(PR, ("%% M =%d K = %d N = %d \n", M, K, N));
            PRLEVEL(PR, ("%% lda =%d ldb =%d\n", lda, ldb));
            PRLEVEL(PR, ("%% j2 =%ld j1=%ld\n", j2, j1));
            PRLEVEL(PR, ("\n %%"));
#endif
            paru_tasked_dgemm(f, M, N, K, A, lda, B, ldb, 1, C, ldc, Work, Num);
            // printf ("%d %d %d ",M ,N, K);
            // printf ("%d %d %d\n ",lda ,ldb, ldc);
#ifdef COUNT_FLOPS
            // printf("dgemm adding to flop count %ld\n", M*N*2);
            //#pragma omp atomic
            // Num->flp_cnt_real_dgemm += (double)2 * M * N * K;
#ifndef NDEBUG
            PRLEVEL(PR, ("\n%% FlopCount Dgemm factorize %d %d %d ", M, N, K));
            PRLEVEL(PR, ("%d %d %d \n", M, N, K));
#endif
#endif
        }

#ifndef NDEBUG
        if (j2 < fp)
        {
            PRLEVEL(PR, ("%% Pivotal Front After Dgemm: %ld x %ld\n %%", fp,
                         rowCount));
            for (Int r = 0; r < rowCount; r++)
            {
                PRLEVEL(PR, ("%% %ld\t", frowList[r]));
                for (Int c = 0; c < fp; c++)
                    PRLEVEL(PR, (" %2.5lf\t", F[c * rowCount + r]));
                PRLEVEL(PR, ("\n"));
            }
        }
#endif
    }
    return 0;
}
