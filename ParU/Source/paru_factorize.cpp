/** =========================================================================  /
 * =======================  paru_factorize  =================================  /
 * ========================================================================== */


#include "Parallel_LU.hpp"

extern "C" void dgetrf_( BLAS_INT *dim1, BLAS_INT *dim2, double *a, 
        BLAS_INT *lda, BLAS_INT *ipiv, BLAS_INT *info);


Int paru_factorize (double *F, Int *fsRowList, Int lm, Int ln,
        BLAS_INT *ipiv)
{
    DEBUGLEVEL(0);

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


#ifndef NDEBUG  // Printing the swap permutation
    p = 1;
    // ATTENTION: ipiv is 1 based
    PRLEVEL (p, ("swap permutation:\n"));
    for (int i = 0; i < m; i++){
        PRLEVEL (p, ("ipiv[%d] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (p, ("\n"));
#endif 

    PRLEVEL (1, (" m=%d n=%d\n", m, n));

    
    for (int i = 0; i < n; i++){
        Int tmp;
        // swap (fsRowList[ipiv [i]], fsRowList[i] ) and it is off by one
        PRLEVEL (1, ("ipiv[%d] =%d\n",i, ipiv[i]));
        ASSERT (ipiv [i] <= m);
        ASSERT (ipiv [i] > 0);
        fsRowList[ipiv [i]-1] = fsRowList [i];
        fsRowList [i] = tmp;
    }


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
