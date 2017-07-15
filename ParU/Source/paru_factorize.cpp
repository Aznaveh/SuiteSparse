#include "Parallel_LU.hpp"

extern "C" void dgetrf_(Int* dim1, Int* dim2, double* a, Int* lda,
        Int* ipiv, Int* info);


Int paru_factorize (double *F, Int m, Int n, Int *ipiv)
{
    DEBUGLEVEL(1);
#ifndef NDEBUG  // Printing the pivotal front
    Int p = 0;
    for (Int r = 0; r < m; r++){
        for (Int c = 0; c < n; c++){
            PRLEVEL (p, (" %3.1lf\t", F[c*m+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
    Int info;
    Int lda = m;

    PRLEVEL (1, ("ipiv =%p\n", ipiv));
    PRLEVEL (1, ("F=%p\n", F));
    
    
#ifndef NDEBUG  // Initializing permutation for debug
    for (Int i = 0; i < lda ; i++){
        ipiv [i] = -1;
    }
#endif


    dgetrf_(&m, &n, F, &lda, ipiv, &info);

    for (Int i = 0; i < lda ; i++){
        PRLEVEL (1, ("ipiv[%ld] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (1, ("\n"));

    PRLEVEL (1, ("info =%ld\n", info));

#ifndef NDEBUG  // Printing the LU decomposition
    p = 1;
    for (Int r = 0; r < m; r++){
        for (Int c = 0; c < n; c++){
            PRLEVEL (p, (" %3.1lf\t", F[c*m+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
 
    return 0;
}
