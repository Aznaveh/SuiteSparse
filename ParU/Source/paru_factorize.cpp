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
    DEBUGLEVEL(1);
    PRLEVEL(1, ("executing front %ld\n", f));

    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *Parent = LUsym->Parent;
    ParU_ResultCode myInfo = paru_front(f, paruMatInfo);
    if (myInfo != PARU_SUCCESS)
    {
        return myInfo;
    }
    Int num_rem_children;
    Int daddy = Parent[f];
    if (daddy != -1) //if it is not a root
    {
        #pragma omp critical
        num_rem_children = --num_active_children[daddy];
        
        //These two operations are possible with atomic but race can happen
        //#pragma omp atomic 
        //num_active_children[daddy]--;
        //#pragma omp atomic read
        //num_rem_children = num_active_children[daddy];

        PRLEVEL(1, ("%% finished %ld  Parent has %ld left\n", 
        f, num_active_children[daddy]));
        if (num_rem_children == 0)
        {
            return myInfo = 
                paru_exec(Parent[f], num_active_children, paruMatInfo);
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

    //////////////// Queue based ///////////////////////////////////////////////
    Int nf = LUsym->nf;
    Int *Parent = LUsym->Parent;
    std::vector<Int> num_active_children (nf, 0);
    #pragma omp parallel for
    for (Int f = 0; f < nf; f++)
        if (Parent[f] != -1) 
        {
            #pragma omp atomic
            num_active_children[Parent[f]]++;
        }
#ifndef NDEBUG
    Int PR = 1;
    PRLEVEL(PR, ("Number of children:\n"));
    for (Int f = 0; f < nf; f++)
    PRLEVEL(PR, ("%ld ",num_active_children[f]));
    PR = 1;
#endif
    std::vector<Int> Q;
    for (Int f = 0; f < nf; f++)
        if (num_active_children[f] == 0) Q.push_back(f);
    Int *Depth = LUsym->Depth;
    std::sort(Q.begin(), Q.end(), [&Depth](const Int &a, const Int &b)-> bool 
            {return Depth[a] > Depth[b];});

#ifndef NDEBUG
    Int PR = 0;
    PRLEVEL(PR, ("Leaves with their depth:\n"));
    for (auto l: Q) 
        PRLEVEL(PR, ("%ld(%ld) ",l, Depth[l]));
    PRLEVEL(PR, ("\n"));
    PR = 1;
#endif
    
    #pragma omp parallel
    #pragma omp single nowait
    #pragma omp task untied
    for (Int i = 0; i < (Int)Q.size(); i++)
    {
        Int f = Q[i];
        //printf("poping %ld \n", f);
        Int d = Depth[f];
        #pragma omp task priority(d) 
        {
            ParU_ResultCode myInfo =
                paru_exec(f, &num_active_children[0], paruMatInfo);
            if (myInfo != PARU_SUCCESS)
            {
                #pragma omp critical
                info = myInfo;
            }
        }
    }
    //////////////// Queue based ///END/////////////////////////////////////////
    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% factorization has some problem\n"));
        if (info == PARU_OUT_OF_MEMORY)
            printf("Paru: out of memory during factorization\n");
        else if (info == PARU_SINGULAR)
            printf("Paru: Input matrix is singular\n");
        return info;
    }

    // The following code can be substituted in a sequential case
    // Int nf = LUsym->nf;
    // for (Int i = 0; i < nf; i++)
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
