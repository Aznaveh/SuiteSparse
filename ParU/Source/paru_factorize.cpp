#include "Parallel_LU.hpp"
Int paru_factorize (double *F, Int m, Int n)
{
    DEBUGLEVEL(0);
#ifndef NDEBUG  // Printing the pivotal front
    Int p = 1;
   for (int r = 0; r < m; r++){
        for (Int c = 0; c < n; c++){
            PRLEVEL (p, (" %3.1lf\t", F[c*m+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif


    return 0;
}
