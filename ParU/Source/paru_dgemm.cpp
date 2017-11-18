/** =========================================================================  /
 * =======================  paru_dgemm ======================================  /
 * ========================================================================== */


#include "Parallel_LU.hpp"
                                                        
/*!
 *
 *      This does the outer product
 *        C := alpha*op( A )*op( B ) + beta*C,   
 *        el:= -1*part(pF) * uPart + zero 
 *        TRANSA = 'N'   op(A) = A
 *        TRANSB = 'N'   op(B) = B
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

    DEBUGLEVEL(1);
    PRLEVEL (1, ("rowCount =%ld  ", rowCount));
    PRLEVEL (1, ("colCount =%ld  ", colCount));
    PRLEVEL (1, ("fp =%ld\n", fp));

    BLAS_INT mA = (BLAS_INT) (rowCount - fp); 
    BLAS_INT nB = (BLAS_INT)colCount;
    BLAS_INT nA = (BLAS_INT)fp;

    PRLEVEL (1, ("mA =%d  ", mA));
    PRLEVEL (1, ("nB =%d  ", nB));
    PRLEVEL (1, ("nA =%d\n", nA));

    double alpha = -1;
    BLAS_INT lda = (BLAS_INT) rowCount;
    BLAS_INT ldb = (BLAS_INT) fp;
    BLAS_INT ldc = (BLAS_INT) (rowCount - fp);
    PRLEVEL (1, ("alpha =%lf  ", alpha));
    PRLEVEL (1, ("lda =%ld  ", lda));
    PRLEVEL (1, ("ldb =%ld  ", ldb));
    PRLEVEL (1, ("ldc =%ld\n", ldc));


    // NOTE: beta can be zero, since C has just been calloc'd.
    // Also, C can be malloc'd not calloc'd.
    double beta= 1;
    
    BLAS_DGEMM ("N" ,"N" , &mA, &nB, &nA, &alpha, pF+fp, &lda, uPart, &ldb, &beta,
            el, &ldc);
    return 0;

}
