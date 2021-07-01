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

ParU_ResultCode paru_do_fronts(Int f, paru_matrix *paruMatInfo)
// This routine call paru_front from first(f)...f including f
{
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *first = LUsym->first;

    ParU_ResultCode info;
    ASSERT(first[f] >= 0);
    //#pragma omp task shared(paruMatInfo)
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
    //#pragma omp task shared(paruMatInfo)
    //info = paru_front(f, paruMatInfo);
    return PARU_SUCCESS;
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

    //TODO my plan is to make do_fronts a task parallel region
    //#pragma omp parallel shared(paruMatInfo)
    {
        //#pragma omp single
        info = paru_do_fronts(nf-1, paruMatInfo);
    }

    // for (Int i = 0; i < nf; i++)
    // {
    //     PRLEVEL(1, ("%% Wroking on front %ld\n", i));
    //     info = paru_front(i, paruMatInfo);
    //     if (info != PARU_SUCCESS)
    //     {
    //         PRLEVEL(1, ("%% A problem happend in %ld\n", i));
    //         return info;
    //     }
    // }

    paruMatInfo->my_time = omp_get_wtime() - my_start_time;
    return PARU_SUCCESS;
}
