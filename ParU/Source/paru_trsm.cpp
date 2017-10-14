/** =========================================================================  /
 * =======================  paru_trsm =======================================  /
 * ========================================================================== */


#include "Parallel_LU.hpp"
                                                        
/*!
 * l11*u12=a12 and u12 is unkonwn  so we need this:
 *          op( A ) * X = alpha*B --> SIDE = 'L' or 'l'
 *         part(pF) * X = 1 * upart
 * UPLO = 'L' or 'l'; A is a lower triangular matrix.   
 *  TRANSA = 'N' or 'n'   op( A ) = A.
 *  DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *  M     M specifies the number of rows of B.
 *  N     N specifies the number of columns of B.
 *  TRANSA = 'N' or 'n'   op( A ) = A.
 *  DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *  LDA  leading dimension of A. 
 */



/// Already defined in /SuiteSparse/CHOLMOD/Include/cholmod_blas.h
//          use BLAS_DTRSM insted
//extern "C"  int dtrsm_(char *side, char *uplo, char *transa, char *diag, 
//                           int *m, int *n, double *alpha, double *a, int *lda, 
//                                                         double *b, int *ldb);
 
Int paru_trsm(double *pF, double *uPart, Int fp, Int rowCount, Int colCount){
    
    /*! TODO: Long int to int conversion     */
    BLAS_INT mB = fp;
    BLAS_INT nB = colCount;
    double alpha = 1;
    BLAS_INT lda = rowCount;
    BLAS_INT ldb = colCount;

    BLAS_DTRSM ("L" ,"L" ,"N" ,"U", &mB, &nB, &alpha, pF, &lda, 
            uPart, &ldb);
    return 0;

}
