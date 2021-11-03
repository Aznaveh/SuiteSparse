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

#ifndef NDEBUG
Int ntasks = 0 ;
#endif

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
    paruMatInfo->res = PARU_SUCCESS;
    if (stree_flop_bound[f] < TASK_FL_THRESHOLD)
    {
        Int *first = LUsym->first;
        PRLEVEL(1, ("%% Sequential %ld - %ld is small (%lf)\n", first[f], f,
                    stree_flop_bound[f]));
        ASSERT(first[f] >= 0);
        for (Int i = first[f]; i <= f; i++)
        {
            PRLEVEL(2, ("%% Wroking on front %ld\n", i));
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

#ifndef NDEBUG
        PRLEVEL(1, ("%% tasks are generating for children of %ld(%lf)\n", f,
                    stree_flop_bound[f]));
        if (Childp[f + 1] - Childp[f] > 100)
            PRLEVEL(1, ("%% lots of children here\n"));
#endif

        Int nchild = (Childp[f + 1] - Childp[f]) ;

        if (nchild == 1)
        {
            PRLEVEL(1, ("%% Just one child for %ld\n",f));

            Int i = Childp[f] ;
            {
                {
                    ParU_ResultCode myInfo = paru_do_fronts(Child[i], paruMatInfo);
                    if (myInfo != PARU_SUCCESS)
                    {
                        // PRLEVEL(1, ("%% A problem happend in %ld\n", i));
                        info = myInfo;
                        return info;
                    }
                }
            }
            // I could also use it but it doesnt work with cancel
            info = paru_front(f, paruMatInfo);
            if (info != PARU_SUCCESS)
            {
                PRLEVEL(1, ("%% A problem happend in %ld\n", f));
                return info;
            }

        }
        else
        {

#if 1

            // at least 2 children

#ifndef NDEBUG
            #pragma omp atomic
            ntasks += nchild ;
#endif

            #pragma omp taskloop grainsize(1) shared(info)
            for (Int i = Childp[f]; i <= Childp[f + 1] - 1; i++)
            {
                {
                    ParU_ResultCode myInfo = 
                        paru_do_fronts(Child[i], paruMatInfo);
                    if (myInfo != PARU_SUCCESS)
                    {
                        // PRLEVEL(1, ("%% A problem happend in %ld\n", i));
                        #pragma omp critical
                        info = myInfo;
                        //#pragma omp cancel taskgroup
                        // return info;
                    }
                }
            }

#else

            // at least 2 children

#ifndef NDEBUG
            #pragma omp atomic
            ntasks += (nchild - 1) ;
#endif

            #pragma omp taskloop nogroup
            for (Int i = Childp[f] + 1; i <= Childp[f + 1] - 1; i++)
            {
                {
                    ParU_ResultCode myInfo = paru_do_fronts(Child[i], paruMatInfo);
                    if (myInfo != PARU_SUCCESS)
                    {
                        // PRLEVEL(1, ("%% A problem happend in %ld\n", i));
                        info = myInfo;
                        #pragma omp cancel taskgroup
                        // return info;
                    }
                }
            }

            Int i = Childp[f];
            {
                {
                    ParU_ResultCode myInfo = paru_do_fronts(Child[i], paruMatInfo);
                    if (myInfo != PARU_SUCCESS)
                    {
                        // PRLEVEL(1, ("%% A problem happend in %ld\n", i));
                        info = myInfo;
                        // return info;
                    }
                }
            }

            #pragma omp taskwait
#endif

            if (info != PARU_SUCCESS)
            {
                return info;
            }
            // I could also use it but it doesnt work with cancel
            info = paru_front(f, paruMatInfo);
            if (info != PARU_SUCCESS)
            {
                PRLEVEL(1, ("%% A problem happend in %ld\n", f));
                // return info;
            }

        }
    }
    return info;
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

    // do_fronts generate a task parallel region
#ifndef NDEBUG
    Int *Parent = LUsym->Parent;
#endif
    if (LUsym->num_roots == 1)
        info = paru_do_fronts(LUsym->roots[0], paruMatInfo);
    else
    {

        #pragma omp taskloop nogroup 
        for (Int i = 0; i < LUsym->num_roots; i++)
        {
            Int r = LUsym->roots[i];
            ASSERT(Parent[r] == -1);
            ParU_ResultCode myInfo = paru_do_fronts(r, paruMatInfo);
            if (myInfo != PARU_SUCCESS)
            {
                // PRLEVEL(1, ("%% A problem happend in %ld\n", i));
                info = myInfo;
                //#pragma omp cancel taskgroup
                // return info;
            }
        }
    }

    if (info != PARU_SUCCESS)
    {
        PRLEVEL(1, ("%% factorization has some problem\n"));
        return info;
    }
#ifndef NDEBUG
    else 
    {
        PRLEVEL(1, ("%% factorization is done with %ld tasks\n",ntasks));
    }
#endif

    // The following code can be substituted in a sequential case
    //Int nf = LUsym->nf;
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
