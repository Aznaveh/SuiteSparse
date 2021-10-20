////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_cumsum ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!
 * @brief   Overwrite a vector of length n with its cumulitive sum of length
 *          n+1.
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define base (1024 * 1024 * 256)
Int paru_cumsum(Int n, Int *X)
{  // n is size, X is size n and in/out
    Int tot = 0;
    if (X == NULL) return tot;
    if (n < base)
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
    #pragma omp parallel
    {
        #pragma omp single
        {
            #pragma omp task shared(sum)
            sum = paru_cumsum(mid, X);
            #pragma omp task 
            paru_cumsum(n - mid, X+mid);
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
