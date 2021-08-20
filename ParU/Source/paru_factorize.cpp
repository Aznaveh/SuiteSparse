/** =========================================================================  /
 * =======================  paru_factorize ==================================  /
 * ==========================================================================  /
 * @brief    get a matrix and factorize it
 *      Allocate space for paruMatInfo
 *      the user should free the space
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"
#define TASK_FL_THRESHOLD 1024 * 1024

ParU_ResultCode paru_do_fronts(Int f, paru_matrix *paruMatInfo)
// This routine call paru_front from first(f)...f including f
// This routine is called recursively to make tasks
{
    DEBUGLEVEL(1);
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    ParU_ResultCode info;

    // double *front_flop_bound = LUsym->front_flop_bound;
    double *stree_flop_bound = LUsym->stree_flop_bound;

    info = PARU_SUCCESS;
    if (stree_flop_bound[f] < TASK_FL_THRESHOLD)
    {
        Int *first = LUsym->first;
        PRLEVEL(1, ("%% Sequential %ld - %ld is small\n", first[f], f));
        ASSERT(first[f] >= 0);
        for (Int i = first[f]; i <= f; i++)
        {
            PRLEVEL(1, ("%% Wroking on front %ld\n", i));
            info = paru_front(i, paruMatInfo);
            if (info != PARU_SUCCESS)
            {
                PRLEVEL(1, ("%% A problem happend in %ld\n", i));
                return info;
            }
        }
    }
    else
    {
        Int *Childp = LUsym->Childp;
        Int *Child = LUsym->Child;

        PRLEVEL(1, ("%% tasks are generating for children of %ld\n", f));
        if (Childp[f + 1] - Childp[f] > 100)
            printf("%% lots of children here\n");
        //#pragma omp taskgroup
        for (Int i = Childp[f]; i <= Childp[f + 1] - 1; i++)
        {
            //   #pragma omp task default(none) shared(paruMatInfo, Child, info)
            //   firstprivate(i)
            {
                ParU_ResultCode myInfo = paru_do_fronts(Child[i], paruMatInfo);
                if (myInfo != PARU_SUCCESS)
                {
                    PRLEVEL(1, ("%% A problem happend in %ld\n", i));
                    info = myInfo;
                    //#pragma omp cancel taskgroup
                    // return info;
                }
            }
        }
        // I could also use it but it doesnt work with cancel
        //#pragma omp taskwait
        info = paru_front(f, paruMatInfo);
        if (info != PARU_SUCCESS)
        {
            PRLEVEL(1, ("%% A problem happend in %ld\n", f));
            // return info;
        }
    }
    return info;
}

ParU_ResultCode paru_factorize(cholmod_sparse *A, paru_symbolic *LUsym,
                               paru_matrix **paruMatInfo_handle)
{
    DEBUGLEVEL(0);
    double my_start_time = omp_get_wtime();
    if (A == NULL)
    {
        printf("Paru: input matrix is invalid\n");
        return PARU_INVALID;
    }

    if (A->xtype != CHOLMOD_REAL)
    {
        printf("Paru: input matrix must be real\n");
        return PARU_INVALID;
    }

    if (LUsym == NULL)
    {
        return PARU_INVALID;
    }

    int scale = 0;

    paru_matrix *paruMatInfo;
    paruMatInfo = *paruMatInfo_handle;

    ParU_ResultCode info;
    info = paru_init_rowFronts(&paruMatInfo, A, scale, LUsym);
    *paruMatInfo_handle = paruMatInfo;

    PRLEVEL(1, ("%% init_row is done\n"));
    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% init_row has a problem\n"));
        return info;
    }

    Int nf = paruMatInfo->LUsym->nf;

    //// do_fronts generate a task parallel region
    // Int *Parent = LUsym->Parent;
    // #pragma omp taskgroup
    // for (Int i = 0; i < nf; i++)
    // {
    //     #pragma omp single nowait
    //     if (Parent[i] == -1)
    //     {
    //  #pragma omp task default(none) shared(paruMatInfo, info) firstprivate(i)
    //         {
    //             ParU_ResultCode myInfo = paru_do_fronts(i, paruMatInfo);
    //             if (myInfo != PARU_SUCCESS)
    //             {
    //                 PRLEVEL(1, ("%% A problem happend in %ld\n", i));
    //                 info = myInfo;
    //                 #pragma omp cancel taskgroup
    //                 // return info;
    //             }
    //         }
    //     }
    // }

    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% factorization has some problem\n"));
        return info;
    }

    // The following code can be substituted in a sequential case
    for (Int i = 0; i < nf; i++)
    {
        PRLEVEL(1, ("%% Wroking on front %ld\n", i));
        info = paru_front(i, paruMatInfo);
        if (info != PARU_SUCCESS)
        {
            PRLEVEL(1, ("%% A problem happend in %ld\n", i));
            return info;
        }
    }
    // TODO: temporary to test perm
    paru_perm(paruMatInfo);
    Int m = LUsym->m;
    double b[m];
    double x[m];
    double xt[m];
    for (Int i = 0; i < m; ++i) b[i] = i + 1;
    paru_apply_perm(LUsym->Pfin, b, x, m);  // x = p (b)
    if (paruMatInfo->scale_row)
        paru_apply_scale (paruMatInfo->scale_row, LUsym->Ps, x, m, LUsym->n1);
    paru_lsolve(paruMatInfo, x);
    paru_usolve(paruMatInfo, x);
    paru_apply_inv_perm(LUsym->Qfill, x, xt, m);  // xt = qinv (x)


    printf ("x = [ ");
    for (Int i = 0; i < MIN(m,10); ++i) 
        printf ("%lf ",xt[i]);
    printf (" ...]\n");

    for (Int i = 0; i < m; ++i) b[i] *= -1;
    //    b[i] = 0;
    paru_gaxpy(A, xt, b);
    double res = paru_vec_1norm(b, m);
    double weighted_res = res / (paru_spm_1norm(A) * paru_vec_1norm(xt, m));
    printf("Residual is |%.2lf| and weigheted residual is |%.2f|.\n", 
            log10(res), log10(weighted_res) );

    paruMatInfo->my_time = omp_get_wtime() - my_start_time;
    return PARU_SUCCESS;
}
