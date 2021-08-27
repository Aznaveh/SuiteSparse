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


    paru_matrix *paruMatInfo;
    paruMatInfo = *paruMatInfo_handle;

    ParU_ResultCode info;
    info = paru_init_rowFronts(&paruMatInfo, A, LUsym);
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
    paru_perm(paruMatInfo); // to form the final permutation
    Int m = LUsym->m;
    paruMatInfo->my_time = omp_get_wtime() - my_start_time;
    printf("time = %lf\n",paruMatInfo->my_time);

    double *b = (double *)paru_alloc(m, sizeof(double));
    
    //double b[m];
    //double x[m];
    //double xt[m];
    for (Int i = 0; i < m; ++i) b[i] = i + 1;

    //paru_apply_perm(LUsym->Pfin, b, x, m);  // x = p (b)
    //if (paruMatInfo->scale_row)
    //    paru_apply_scale (paruMatInfo->scale_row, LUsym->Ps, x, m, LUsym->n1);
    //paru_lsolve(paruMatInfo, x);
    //paru_usolve(paruMatInfo, x);
    //paru_apply_inv_perm(LUsym->Qfill, x, xt, m);  // xt = qinv (x)


    //printf ("x = [ ");
    //for (Int i = 0; i < MIN(m,10); ++i) 
    //    printf ("%lf ",xt[i]);
    //printf (" ...]\n");

    //paru_gaxpy(A, xt, b, -1);
    //double res = paru_vec_1norm(b, m);
     
    double Res[4];
    paru_residual(A, paruMatInfo, b, Res);
    for (Int i = 0; i < m; ++i) b[i] = i + 1;
    paru_backward(A, paruMatInfo, b, Res);
    // printf("Hooh %lf\n", log10(Res[1]));
    //double res = paru_residual(A, paruMatInfo, b);
    //double weighted_res = res / (paru_spm_1norm(A) * paru_vec_1norm(xt, m));
    //printf("Residual is |%.2lf| and weigheted residual is |%.2f|.\n", 
    //        log10(res), log10(weighted_res) );

    
    //forward error
    //for (Int i = 0; i < m; ++i) b[i] = i + 1;
    //paru_gaxpy(A, b, x, 1);
    //paru_apply_perm(LUsym->Pfin, b, x, m);  // x = p (b)
    //if (paruMatInfo->scale_row)
    //    paru_apply_scale (paruMatInfo->scale_row, LUsym->Ps, x, m, LUsym->n1);
    //paru_lsolve(paruMatInfo, x);
    //paru_usolve(paruMatInfo, x);
    //paru_apply_inv_perm(LUsym->Qfill, x, xt, m);  // xt = qinv (x)
    //paru_gaxpy(A, xt, b, -1);

    //res = paru_vec_1norm(b, m);
    //weighted_res = res / (paru_spm_1norm(A) * paru_vec_1norm(xt, m));
    //printf("forward eorror is |%.2lf| and weigheted forward error is |%.2f|.\n", 
    //        log10(res), log10(weighted_res) );


    paru_free (m, sizeof(double), b);
    return PARU_SUCCESS;
}
