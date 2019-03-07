/** =========================================================================  /
 * =======================  paru_factorize  =================================  /
 * ========================================================================== 
 * @brief Doing the BLAS factorization 
 *         
 * @author Aznaveh
 **/
 

#include "Parallel_LU.hpp"

extern "C" void dgetrf_( BLAS_INT *dim1, BLAS_INT *dim2, double *a, 
        BLAS_INT *lda, BLAS_INT *ipiv, BLAS_INT *info);

Int paru_dgetrf (double *F, Int *fsRowList, Int lm, Int ln,
        BLAS_INT *ipiv)
{
    DEBUGLEVEL(0);

    BLAS_INT m = (BLAS_INT) lm;
    BLAS_INT n = (BLAS_INT) ln;


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
    if (m < n ){
        PRLEVEL (0, ("%%!!!!! FAIL m= %d  n= %d\n", m, n));
        return -1;
    }
#endif

#ifndef NDEBUG  // Printing the list of rows
    p = 1;
    PRLEVEL (p, ("Befor factorization (inside factorize): \n"));
    for (int i = 0; i < m ; i++){
        PRLEVEL (p, ("fsRowList [%d] =%d\n",i, fsRowList [i]));
    }
    PRLEVEL (p, ("\n"));
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

    
    // swap (fsRowList[ipiv [i]], fsRowList[i] ) and it is off by one
    for (Int i = 0; i < n; i++){
        PRLEVEL (1, ("ipiv[%d] =%d\n", i, ipiv[i]));
        ASSERT (ipiv [i] > 0);
        ASSERT (ipiv [i] <= m);
        Int tmp =  fsRowList[ipiv [i]-1];
        PRLEVEL (1, ("tmp =%ld\n", tmp));
        ASSERT (tmp >= 0);

        fsRowList[ipiv [i]-1] = fsRowList [i];
        fsRowList [i] = tmp;
    }

#ifndef NDEBUG   // Printing the LU decomposition
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
        printf("%%Some problem in factorization\n");
        return info;
    }
    return 0;
}

Int paru_factorize(double *F, Int *fsRowList, Int rowCount, Int fp, 
        paru_matrix *paruMatInfo)
        
{

    work_struct *Work =  paruMatInfo->Work;
    Int *row_degree_bound = paruMatInfo->row_degree_bound;
    BLAS_INT *ipiv = (BLAS_INT*) (Work->scratch+rowCount);
    return paru_dgetrf (F , fsRowList, rowCount, fp, ipiv);
    return 0;
}


