/** =========================================================================  /
 * =======================  paru_dgemm ======================================  /
 * ========================================================================== */


#include "Parallel_LU.hpp"
                                                        
/*!
 * @brief   Overwrite a vector of length n with its cumulitive sum of length
 *          n+1.  This is very similar to cumsum in SPQR
 * 
 * @author Aznaveh
 */
Int paru_cumsum (Int n, Int *X)
{   // n is size, X is size n and in/out
    Int tot = 0;
    if (X != NULL)
    {
        for(Int k = 0; k < n ; k++)
        {
           X[k] += tot;     // tot = sum (X[0:k-1])
           tot = X[k];
        }
    }
    return tot;
}


