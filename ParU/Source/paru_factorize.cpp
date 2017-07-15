#include "Parallel_LU.hpp"

extern "C" void dgetrf_(Int* dim1, Int* dim2, double* a, Int* lda,
        int* ipiv, Int* info);


Int paru_factorize (double *F, Int m, Int n, int *ipiv)
{
    DEBUGLEVEL(0);
    PRLEVEL (1, (" %ld x %ld\n", m, n));
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
    for (int i = 0; i < lda ; i++){
        ipiv [i] = -1;
    }
#endif


    dgetrf_(&m, &n, F, &lda, ipiv, &info);

#ifndef NDEBUG  // Printing the permutation
    p = 0;
    // ATTENTION: ipiv is 1 based
    for (int i = 0; i < n; i++){
        PRLEVEL (p, ("ipiv[%d] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (p, ("\n"));

    // Printing the LU decomposition
    p = 0;
    for (Int r = 0; r < m; r++){
        for (Int c = 0; c < n; c++){
            PRLEVEL (p, (" %3.1lf\t", F[c*m+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif

    PRLEVEL (1, ("info = %ld\n", info));
    if (info !=0 ){
        printf("Some problem in factorization\n");
        return info;
    }
    return 0;
}
