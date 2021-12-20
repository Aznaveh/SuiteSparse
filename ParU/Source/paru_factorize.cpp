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
#define TASK_FL_THRESHOLD (double(1024 * 1024))


ParU_ResultCode paru_exec(Int f, 
        Int* num_active_children, 
        paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("executing front %ld\n", f));

    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *Parent = LUsym->Parent;
    Int *Childp = LUsym->Childp;
    ParU_ResultCode myInfo = paru_front(f, paruMatInfo);
    if (myInfo != PARU_SUCCESS)
    {
        return myInfo;
    }
    Int num_rem_children;
    Int daddy = Parent[f];
    Int num_original_children = Childp[daddy+1] - Childp[daddy];
    if (daddy != -1) //if it is not a root
    {
        if (num_original_children != 1)
        {
            #pragma omp atomic capture
            { 
                num_active_children[daddy]--;
                num_rem_children = num_active_children[daddy];
            }

            PRLEVEL(1, ("%% finished %ld  Parent has %ld left\n", 
                        f, num_active_children[daddy]));
            if (num_rem_children == 0)
            {
                return myInfo = 
                    paru_exec(Parent[f], num_active_children, paruMatInfo);
            }
        }
        else //I was the only spoiled kid in the family; no need for critical
        {
            return myInfo = 
                paru_exec(Parent[f], num_active_children, paruMatInfo);
        }
    }
    return myInfo;
}

ParU_ResultCode paru_exec_tasks(Int t, 
        Int* task_num_child, 
        Int* task_map, 
        Int* task_parent,
        Int* task_num_child_nochange, 
        paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    Int daddy = task_parent[t];
    Int num_original_children = task_num_child_nochange[daddy];
    PRLEVEL(1, ("executing task %ld fronts %ld-%ld (%ld children)\n",t+1, 
                task_map[t]+1+1,task_map[t+1]+1, num_original_children));
    ParU_ResultCode myInfo;
    for (Int f = task_map[t]+1; f <= task_map[t+1]; f++)
    {
        myInfo = paru_front(f, paruMatInfo);
        if (myInfo != PARU_SUCCESS)
        {
            return myInfo;
        }
    }
    Int num_rem_children;

#ifndef NDEBUG
    if (daddy == -1)
        PRLEVEL(1, ("%% finished task root(%ld)\n", t));
#endif
 
    if (daddy != -1) //if it is not a root
    {
        if (num_original_children != 1)
        {
            #pragma omp atomic capture
            { 
                task_num_child[daddy]--;
                num_rem_children = task_num_child[daddy];
            }

            PRLEVEL(1, ("%% finished task %ld(%ld,%ld)  Parent has %ld left\n", 
                        t+1, task_map[t]+1+1,task_map[t+1]+1,task_num_child[daddy]));
            if (num_rem_children == 0)
            {
                PRLEVEL(1, ("%% task %ld executing its parent %ld\n", t+1, task_parent[t]+1));
                return myInfo = 
                    paru_exec_tasks(task_parent[t], task_num_child, task_map, task_parent, 
                            task_num_child_nochange, paruMatInfo);
            }
        }
        else //I was the only spoiled kid in the family; no need for critical
        {
            PRLEVEL(1, ("%% task %ld only child executing its parent %ld\n", t+1, task_parent[t]+1));
            return myInfo = 
                    paru_exec_tasks(task_parent[t], task_num_child, task_map, task_parent, 
                            task_num_child_nochange, paruMatInfo);
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
    // printf ("Starting init row\n");
    info = paru_init_rowFronts(&paruMatInfo, A, LUsym);
    // printf ("Finishing init row\n");
    *paruMatInfo_handle = paruMatInfo;

    PRLEVEL(1, ("%% init_row is done\n"));
    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% init_row has a problem\n"));
        return info;
    }

    Int nf = LUsym->nf;
    Int *Childp = LUsym->Childp;
    Int *Depth = LUsym->Depth;
 
    //////////////// Making task tree //////////////////////////////////////////
    //TODO fix jungles
    Int *Parent = LUsym->Parent;
    std::vector<Int> task_helper(nf, 0);
    std::vector<Int> task_map;
    task_map.push_back(-1);
    //double *front_flop_bound = LUsym->front_flop_bound;
    double *stree_flop_bound = LUsym->stree_flop_bound;
    const double flop_thresh = 1024;
    Int num_task_nodes = 0;
 
    
    for (Int f = 0; f < nf; f++)
    {
        Int daddy = Parent[f];
        if (daddy == -1) 
        {
            task_helper[f] = num_task_nodes++;
            task_map.push_back(f);
            continue;
        }
        Int num_sibl = Childp[daddy+1] - Childp[daddy] - 1;
        if (num_sibl == 0 ) //getting rid of chains
        {
            task_helper[f] = -1;
            continue;
        }
        //Int num_child = Childp[f+1] - Childp[f];
        if (stree_flop_bound[f] < flop_thresh)
        {//getting rid of  small tasks
            task_helper[f] = -1;
        }
        else
        {
            task_helper[f] = num_task_nodes++;
            task_map.push_back(f);
        }
    }
    std::vector<Int> task_parent(num_task_nodes,-1);
    for (Int i = 0; i < num_task_nodes; i++) 
    {
        Int node = task_map[i+1];
        Int parent = Parent[node];
        //printf("i=%ld node=%ld parent=%ld\n",i,node,parent);
        if (parent == -1)  
        {
            task_parent[i] = -1;
            continue;
        }
        while (task_helper[parent] < 0 )
            parent = Parent[parent];
        task_parent[i] = task_helper[parent];
    }
    std::vector<Int> task_num_child(num_task_nodes,0);
    for (Int t = 0; t < num_task_nodes; t++)
    {
        Int par = task_parent[t];
        if (par != -1)
            task_num_child[par]++;
    }
    //std::vector<Int> task_num_child_nochange(task_num_child);
    std::vector<Int> task_num_child_nochange(num_task_nodes);
    for (Int i = 0; i < (Int)task_num_child.size(); i++) 
        task_num_child_nochange[i] = task_num_child[i];

    std::vector<Int> task_Q;
    for (Int t = 0; t < num_task_nodes; t++)
    {
        if (task_num_child[t] == 0) task_Q.push_back(t);
    }
    
    std::sort(task_Q.begin(), task_Q.end(), 
            [&Depth, &task_map](const Int &t1, const Int &t2)-> bool 
            {return Depth[task_map[t1]+1] > Depth[task_map[t2]+1];});

//    printf("%% Task tree helper:\n");
//    for (Int i = 0; i < (Int)task_helper.size(); i++) 
//        printf("%ld)%ld ", i, task_helper[i]);
//    printf("\n%% tasknodes map (%ld):\n",num_task_nodes);
//    for (Int i = 0; i < (Int)task_map.size(); i++) 
//        printf("%ld ", task_map[i]);
//    printf("\n%% tasktree :\n");
//    for (Int i = 0; i < (Int)task_parent.size(); i++) 
//        printf("%ld ", task_parent[i]);
//    printf("\n%% task_num_child:\n");
//    for (Int i = 0; i < (Int)task_num_child.size(); i++) 
//        printf("%ld ", task_num_child[i]);
//    printf("\n%% tasks\n");
//    for (Int t = 0; t < num_task_nodes; t++)
//    {
//        printf("%ld[%ld-%ld] ",t, task_map[t]+1, task_map[t+1]);
//    }
//    printf("\n%% task_Q:\n");
//    for (Int i = 0; i < (Int)task_Q.size(); i++) 
//    {
//        Int t = task_Q[i];
//        printf("%ld[%ld-%ld] ",t, task_map[t]+1, task_map[t+1]);
//    }
//    printf("\n");

    #pragma omp parallel
    #pragma omp single nowait
    #pragma omp task untied
    for (Int i = 0; i < (Int)task_Q.size(); i++)
    {
        Int t = task_Q[i];
        //printf("poping %ld \n", f);
        Int d = Depth[task_map[t]+1];
        #pragma omp task priority(d) 
        {
            ParU_ResultCode myInfo =
                paru_exec_tasks(t, &task_num_child[0], &task_map[0],
                        &task_parent[0], &task_num_child_nochange[0],
                        paruMatInfo);
            if (myInfo != PARU_SUCCESS)
            {
                #pragma omp critical
                info = myInfo;
            }
        }
    }


    //////////////// Making task tree ///End////////////////////////////////////

    //////////////// Queue based ///////////////////////////////////////////////
//    std::vector<Int> num_active_children (nf, 0);
//    #pragma omp parallel for
//    for (Int f = 0; f < nf; f++)
//        num_active_children[f] = Childp[f+1] - Childp[f];
//#ifndef NDEBUG
//    Int PR = 1;
//    PRLEVEL(PR, ("Number of children:\n"));
//    for (Int f = 0; f < nf; f++)
//        PRLEVEL(PR, ("%ld ",num_active_children[f]));
//    PRLEVEL(PR, ("\n"));
//    PR = 1;
//#endif
//    std::vector<Int> Q;
//    for (Int f = 0; f < nf; f++)
//        if (num_active_children[f] == 0) Q.push_back(f);
//    //Int *Depth = LUsym->Depth;
//    std::sort(Q.begin(), Q.end(), [&Depth](const Int &a, const Int &b)-> bool 
//            {return Depth[a] > Depth[b];});
//
//#ifndef NDEBUG
//    PR = 0;
//    PRLEVEL(PR, ("%ld Leaves with their depth:\n", Q.size()));
//    for (auto l: Q) 
//        PRLEVEL(PR, ("%ld(%ld) ",l, Depth[l]));
//    PRLEVEL(PR, ("\n"));
//    PR = 1;
//#endif
//    
//    #pragma omp parallel
//    #pragma omp single nowait
//    #pragma omp task untied
//    for (Int i = 0; i < (Int)Q.size(); i++)
//    {
//        Int f = Q[i];
//        //printf("poping %ld \n", f);
//        Int d = Depth[f];
//        #pragma omp task priority(d) 
//        {
//            ParU_ResultCode myInfo =
//                paru_exec(f, &num_active_children[0], paruMatInfo);
//            if (myInfo != PARU_SUCCESS)
//            {
//                #pragma omp critical
//                info = myInfo;
//            }
//        }
//    }
//    //////////////// Queue based ///END/////////////////////////////////////////

    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% factorization has some problem\n"));
        if (info == PARU_OUT_OF_MEMORY)
            printf("Paru: out of memory during factorization\n");
        else if (info == PARU_SINGULAR)
            printf("Paru: Input matrix is singular\n");
        return info;
    }

#ifdef COUNT_FLOPS
    double flop_count = paruMatInfo->flp_cnt_dgemm + paruMatInfo->flp_cnt_dger +
        paruMatInfo->flp_cnt_trsm;
    PRLEVEL (-1, ("Flop count = %.17g\n",flop_count));
#endif
    // The following code can be substituted in a sequential case
    //for (Int i = 0; i < nf; i++)
    //{
    //    if (i %1000 == 0) PRLEVEL(1, ("%% Wroking on front %ld\n", i));

    //    info = paru_front(i, paruMatInfo);
    //    if (info != PARU_SUCCESS)
    //    {
    //        PRLEVEL(1, ("%% A problem happend in %ld\n", i));
    //        return info;
    //    }
    //}

    info = paru_perm(paruMatInfo);  // to form the final permutation

    if (info == PARU_OUT_OF_MEMORY)
    {
        printf("Paru: memory problem after factorizaiton, in perumutaion.\n");
        return info;
    }

    paruMatInfo->my_time = omp_get_wtime() - my_start_time;

    return PARU_SUCCESS;
}
