/** =========================================================================  /
 * =======================  paru_factorize  =================================  /
 * ========================================================================== */


#include "Parallel_LU.hpp"

extern "C" void dgetrf_( BLAS_INT *dim1, BLAS_INT *dim2, double *a, 
        BLAS_INT *lda, BLAS_INT *ipiv, BLAS_INT *info);


Int paru_factorize (double *F, Int lm, Int ln,
        BLAS_INT *ipiv)
{
    DEBUGLEVEL(1);

    BLAS_INT m = lm;
    BLAS_INT n = ln;


    PRLEVEL (1, (" %d x %d\n", m, n));
#ifndef NDEBUG  // Printing the pivotal front before computation
    Int p = 1;
    PRLEVEL (1, ("Befor factorization:\n"));
    for (Int r = 0; r < m; r++){
        for (Int c = 0; c < n; c++){
            PRLEVEL (p, (" %3.4lf\t", F[c*m+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
    BLAS_INT info;
    BLAS_INT lda = m;

    PRLEVEL (1, ("ipiv =%p\n", ipiv));
    PRLEVEL (1, ("F=%p\n", F));


#ifndef NDEBUG  // Initializing permutation; just for debug
    for (int i = 0; i < lda ; i++){
        ipiv [i] = -1;
    }
#endif


    dgetrf_(&m, &n, F, &lda, ipiv, &info);



     ASSERT (m >= n);

     /* changing swap permutation to a real permutation */

    BLAS_INT *tmpPinv = ipiv + n; // using the rest of scratch memory
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
    PRLEVEL (1, (" m=%d n=%d\n", m, n));
    for (int i = 0; i < m; i++){
        int tmp;
        // swap (tmpPinv [ipiv [i]], tmpPinv[i] ) and it is off by one
        PRLEVEL (1, ("ipiv[%d] =%d\n",i, ipiv[i]));
        ASSERT (ipiv [i] <= m);
        tmp = tmpPinv [ipiv [i]-1];
        PRLEVEL (1, ("tmp =%d\n", tmp));
        ASSERT (tmp < m);
        tmpPinv [ipiv [i]-1] = tmpPinv [i];
        tmpPinv [i] = tmp;
    }

    for (int i = 0; i < m; i++) 
        ipiv [i] = tmpPinv[i]; //copying back the important chunck



#ifndef NDEBUG  // Printing the permutation
    p = 1;
    // ATTENTION: ipiv is 1 based
    PRLEVEL (p, ("Real permutation:\n"));
    for (int i = 0; i < m; i++){
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
    if (info != 0 ){
        printf("Some problem in factorization\n");
        return info;
    }
    return 0;
}
