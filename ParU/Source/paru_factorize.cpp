////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_factorize /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief    get a matrix and factorize it
 *      specify the order of eliminating fronts
 *      Allocate space for paruMatInfo
 *      the user should free the space
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

ParU_ResultCode paru_exec_tasks(Int t, Int *task_num_child,
                                paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(1);
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *task_parent = LUsym->task_parent;
    Int daddy = task_parent[t];
    Int *task_map = LUsym->task_map;

    Int num_original_children = 0;
    if (daddy != -1) num_original_children = LUsym->task_num_child[daddy];
    PRLEVEL(1, ("executing task %ld fronts %ld-%ld (%ld children)\n", t,
                task_map[t] + 1, task_map[t + 1], num_original_children));
    ParU_ResultCode myInfo;
    for (Int f = task_map[t] + 1; f <= task_map[t + 1]; f++)
    {
        myInfo = paru_front(f, paruMatInfo);
        if (myInfo != PARU_SUCCESS)
        {
            return myInfo;
        }
    }
    Int num_rem_children;

#ifndef NDEBUG
    if (daddy == -1) PRLEVEL(1, ("%% finished task root(%ld)\n", t));
#endif

    if (daddy != -1)  // if it is not a root
    {
        if (num_original_children != 1)
        {
            #pragma omp atomic capture
            {
                task_num_child[daddy]--;
                num_rem_children = task_num_child[daddy];
            }

            PRLEVEL(1,
                    ("%% finished task %ld(%ld,%ld)  Parent has %ld left\n", t,
                     task_map[t] + 1, task_map[t + 1], task_num_child[daddy]));
            if (num_rem_children == 0)
            {
                PRLEVEL(1,
                        ("%% task %ld executing its parent %ld\n", t, daddy));
                return myInfo =
                           paru_exec_tasks(daddy, task_num_child, paruMatInfo);
            }
        }
        else  // I was the only spoiled kid in the family;
        {
            PRLEVEL(1, ("%% task %ld only child executing its parent %ld\n", t,
                        daddy));
            return myInfo = paru_exec_tasks(daddy, task_num_child, paruMatInfo);
        }
    }
    return myInfo;
}
ParU_ResultCode paru_factorize(cholmod_sparse *A, paru_symbolic *LUsym,
                               paru_matrix **paruMatInfo_handle)
{
    DEBUGLEVEL(1);
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
    paruMatInfo->num_active_tasks = 0;
    Int nf = LUsym->nf;
    //////////////// Using task tree //////////////////////////////////////////
    Int ntasks = LUsym->ntasks;
    Int *task_depth = LUsym->task_depth;
    std::vector<Int> task_Q;
    // This vector changes during factorization
    // std::vector<Int> task_num_child(ntasks);
    // paru_memcpy ( &task_num_child[0] , LUsym->task_num_child,
    //        ntasks * sizeof(Int));

    Int task_num_child[ntasks];
    paru_memcpy(task_num_child, LUsym->task_num_child, ntasks * sizeof(Int));

    for (Int t = 0; t < ntasks; t++)
    {
        if (task_num_child[t] == 0) task_Q.push_back(t);
    }

    std::sort(task_Q.begin(), task_Q.end(),
              [&task_depth](const Int &t1, const Int &t2) -> bool {
                  return task_depth[t1] > task_depth[t2];
              });

    // Int *Depth = LUsym->Depth;
    //    std::sort(task_Q.begin(), task_Q.end(),
    //            [&Depth, &task_map](const Int &t1, const Int &t2)-> bool
    //            {return Depth[task_map[t1]+1] > Depth[task_map[t2]+1];});

    Int max_chain = LUsym->max_chain;
    double chainess = 2;
    double maxchain_ratio = 2;
    printf("ntasks=%ld task_Q.size=%ld\n", ntasks, task_Q.size());
    if (ntasks > 0)
    {
        chainess = (task_depth[task_Q[0]] + 1) / (double)nf;
        maxchain_ratio = (((double)max_chain+1)/nf);
        printf("nf = %ld, deepest = %ld, chainess = %lf max_chain=%ld" 
                " maxchain_ratio =%lf\n", 
                nf,
               task_depth[task_Q[0]],
               chainess,
               max_chain , maxchain_ratio);
    } 
#ifndef NDEBUG
    Int PR = -1;
    Int *task_map = LUsym->task_map;
    PRLEVEL(PR, ("\n%% task_Q:\n"));
    for (Int i = 0; i < (Int)task_Q.size(); i++)
    {
        Int t = task_Q[i];
        PRLEVEL(PR, ("%ld[%ld-%ld](%ld) ", t, task_map[t] + 1, task_map[t + 1],
                     task_depth[t]));
    }
    PRLEVEL(PR, ("\n"));
#endif

    //if (chainess < .6 && maxchain_ratio < .25)
    if (1)
    {
        printf("Parallel\n");
        const Int size = (Int)task_Q.size();
        const Int steps = size == 0 ? 1 : size;
        const Int stages = size / steps + 1;
        Int start = 0;
        PRLEVEL(
            1, ("%% size=%ld, steps =%ld, stages =%ld\n", size, steps, stages));

        for (Int ii = 0; ii < stages; ii++)
        {
            if (start >= size) break;
            Int end = start + steps > size ? size : start + steps;
            PRLEVEL(-1, ("%% doing Queue tasks <%ld,%ld>\n", start, end));
            #pragma omp parallel for
            //#pragma omp single nowait
            //#pragma omp task //untied  //clang seg fault on untied
            for (Int i = start; i < end; i++)
            // for (Int i = 0; i < (Int)task_Q.size(); i++)
            {
                Int t = task_Q[i];
                // printf("poping %ld \n", f);
                // Int d = Depth[task_map[t]+1];
                Int d = task_depth[t];
                // #pragma omp task priority(d)
                {
                    #pragma omp atomic update
                    paruMatInfo->num_active_tasks++;

                    ParU_ResultCode myInfo =
                        // paru_exec_tasks(t, &task_num_child[0], paruMatInfo);
                        paru_exec_tasks(t, task_num_child, paruMatInfo);
                    if (myInfo != PARU_SUCCESS)
                    {
                        #pragma omp atomic write
                        info = myInfo;
                    }
                    #pragma omp atomic update
                    paruMatInfo->num_active_tasks--;
                }
            }
            start += steps;
        }
        if (info != PARU_SUCCESS)
        {
            PRLEVEL(1, ("%% factorization has some problem\n"));
            if (info == PARU_OUT_OF_MEMORY)
                printf("Paru: out of memory during factorization\n");
            else if (info == PARU_SINGULAR)
                printf("Paru: Input matrix is singular\n");
            return info;
        }
    }

    //////////////// Making task tree ///End////////////////////////////////////
    // The following code can be substituted in a sequential case
    else
    {
        printf("Sequential\n");
        paruMatInfo->num_active_tasks = 1;
        for (Int i = 0; i < nf; i++)
        {
            // if (i %1000 == 0) PRLEVEL(1, ("%% Wroking on front %ld\n", i));

            info = paru_front(i, paruMatInfo);
            if (info != PARU_SUCCESS)
            {
                PRLEVEL(1, ("%% A problem happend in %ld\n", i));
                return info;
            }
        }
    }

    info = paru_perm(paruMatInfo);  // to form the final permutation

    if (info == PARU_OUT_OF_MEMORY)
    {
        printf("Paru: memory problem after factorizaiton, in perumutaion.\n");
        return info;
    }

    #ifdef COUNT_FLOPS
    double flop_count = paruMatInfo->flp_cnt_dgemm + paruMatInfo->flp_cnt_dger +
        paruMatInfo->flp_cnt_trsm;
    PRLEVEL(-1, ("Flop count = %.17g\n", flop_count));
    #endif
    paruMatInfo->my_time = omp_get_wtime() - my_start_time;
    return PARU_SUCCESS;
}
