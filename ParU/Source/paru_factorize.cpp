/** =========================================================================  /
 * =======================  paru_factorize  =================================  /
 * ========================================================================== */


#include "Parallel_LU.hpp"

extern "C" void dgetrf_(Int* dim1, Int* dim2, double* a, Int* lda,
        int* ipiv, Int* info);


Int paru_factorize (double *F, Int m, Int n,
        int *ipiv, cholmod_common *cc)
{
    DEBUGLEVEL(0);
    PRLEVEL (0, (" %ld x %ld\n", m, n));
#ifndef NDEBUG  // Printing the pivotal front before computation
    Int p = 1;
    PRLEVEL (1, ("Befor factorization:\n"));
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


#ifndef NDEBUG  // Initializing permutation; just for debug
    for (int i = 0; i < lda ; i++){
        ipiv [i] = -1;
    }
#endif


    dgetrf_(&m, &n, F, &lda, ipiv, &info);

//    ASSERT (m >= n);
    if (m < n) {
        ipiv [0] = -1;
        return info;
    }
    /* changing swap permutation to a real permutation */
    int* tmpPinv = (int*) paru_alloc (m, sizeof (int), cc);
#ifndef NDEBUG  // Printing the swap permutation
    p = 1;
    // ATTENTION: ipiv is 1 based
    PRLEVEL (p, ("swap permutation:\n"));
    for (int i = 0; i < m; i++){
        PRLEVEL (p, ("ipiv[%d] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (p, ("\n"));
#endif

    PRLEVEL (1, ("\n"));
    for (int i = 0; i < m; i++) tmpPinv[i] = i;
    for (int i = 0; i < n; i++){
        int tmp;
        // swap (tmpPinv [ipiv [i]], tmpPinv[i] ) and it is off by one
        ASSERT (ipiv [i] <= m);
        tmp = tmpPinv [ipiv [i]-1];
        PRLEVEL (1, ("tmp =%d\n", tmp));
        ASSERT (tmp < m);
        tmpPinv [ipiv [i]-1] = tmpPinv [i];
        tmpPinv [i] = tmp;
    }

    for (int i = 0; i < n; i++) 
        ipiv [i] = tmpPinv[i]; //copying back the important chunck

    paru_free (m, sizeof (int), tmpPinv, cc);



#ifndef NDEBUG  // Printing the permutation
    p = 1;
    // ATTENTION: ipiv is 1 based
    PRLEVEL (p, ("Real permutation:\n"));
    for (int i = 0; i < n; i++){
        PRLEVEL (p, ("ipiv[%d] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (p, ("\n"));

    // Printing the LU decomposition
    p = 1;
    PRLEVEL (p, ("After factorization:\n"));
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
