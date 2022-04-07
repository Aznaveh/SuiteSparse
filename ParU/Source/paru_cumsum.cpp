////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_cumsum ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*!
 * @brief   Overwrite a vector of length n with its cumulitive sum of length
 *          n+1.
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

Int paru_cumsum(Int n, Int *X, ParU_Control *Control)
{  // n is size, X is size n and in/out
    Int tot = 0;
    if (X == NULL) return tot;
    
    Int mem_chunk = Control->mem_chunk;
    if (n < mem_chunk)
    {
        for (Int k = 0; k < n; k++)
        {
            X[k] += tot;  // tot = sum (X[0:k-1])
            tot = X[k];
        }
        return tot;
    }
    Int mid = n/2;
    Int sum = 0;
    #pragma omp parallel shared(sum, n, X, Control) firstprivate(mid)
    {
        #pragma omp single
        {
            #pragma omp task 
            sum = paru_cumsum(mid, X, Control);
            #pragma omp task 
            paru_cumsum(n - mid, X+mid, Control);
            #pragma omp taskwait 
            #pragma omp taskloop 
            for (int i = mid; i < n; i ++)
            {
                X[i] += sum;
            }
        }
    }
    return X[n-1];
}
