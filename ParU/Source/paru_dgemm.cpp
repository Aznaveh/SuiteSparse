/** =========================================================================  /
 * =======================  paru_dgemm ======================================  /
 * ========================================================================== */


#include "Parallel_LU.hpp"
                                                        
/*!
 *
 *        C := alpha*op( A )*op( B ) + beta*C,   
 * 
 */



/// Already defined in /SuiteSparse/CHOLMOD/Include/cholmod_blas.h
//          use BLAS_dgemm_
//
// void BLAS_DGEMM (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
//	BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
//	BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc) ;
//

Int paru_dgemm(double *pF, double *uPart, 
        double *el, Int fp, Int rowCount, Int colCount){

    BLAS_INT mA = rowCount-fp;

    BLAS_INT nB = colCount;

    BLAS_INT nA = fp;

    double alpha = -1;
    BLAS_INT lda = rowCount;
    BLAS_INT ldb = colCount;
    BLAS_INT ldc = rowCount-fp;

    double beta= 1;
    
    BLAS_DGEMM ("N" ,"N" , &mA, &nB, &nA, &alpha, pF, &lda, uPart, &ldb, &beta,
            el, &ldc);
    return 0;

}
