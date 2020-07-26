/** =========================================================================  /
 * =======================  paru_dgemm ======================================  /
 * ========================================================================== */
/*!
 *
 * @brief       This does the outer product
 *              C := alpha*op( A )*op( B ) + beta*C,   
 *              el:= -1*part(pF) * uPart + zero 
 *              TRANSA = 'N'   op(A) = A
 *              TRANSB = 'N'   op(B) = B
 * 
 * @author Aznaveh
 */
#include "Parallel_LU.hpp"



/// Already defined in /SuiteSparse/CHOLMOD/Include/cholmod_blas.h
//          use BLAS_dgemm_
//
// void BLAS_DGEMM (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
//	BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
//	BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc) ;
//

Int paru_dgemm(double *pF, double *uPart, 
        double *el, Int fp, Int rowCount, Int colCount)
{

    DEBUGLEVEL(0);
    PRLEVEL (1, ("%% rowCount =%ld  ", rowCount));
    PRLEVEL (1, ("%% colCount =%ld  ", colCount));
    PRLEVEL (1, ("%% fp =%ld\n", fp));

    BLAS_INT mA = (BLAS_INT) (rowCount - fp); 
    BLAS_INT nB = (BLAS_INT)colCount;
    BLAS_INT nA = (BLAS_INT)fp;

    PRLEVEL (1, ("%% mA =%d  ", mA));
    PRLEVEL (1, ("%% nB =%d  ", nB));
    PRLEVEL (1, ("%% nA =%d\n", nA));

#ifndef NDEBUG 
    double *Ap= pF+fp;
    PRLEVEL (1, ("%% A =\n"));
    for(Int i = 0; i < mA; i++)
    {
        PRLEVEL (1, ("%% "));
        for(Int j = 0; j < nA; j++)
            PRLEVEL (1, ("%2.4lf\t",Ap[j*rowCount+i]));
        PRLEVEL (1, ("\n"));
    }

    Int mB = nA;
    double *Bp = uPart;
    PRLEVEL (1, ("%% B =\n"));
    for(Int i = 0; i < mB; i++)
    {
        PRLEVEL (1, ("%% "));
        for(Int j = 0; j < nB; j++)
            PRLEVEL (1, ("%2.4lf\t",Bp[j*mB+i]));
        PRLEVEL (1, ("\n"));
    }

    double *Cp = el;
    PRLEVEL (1, ("%%Before DGEMM C =\n"));
    for(Int i = 0; i < mA; i++)
    {
        PRLEVEL (1, ("%% "));
        for(Int j = 0; j < nB; j++)
        {
            PRLEVEL (1, ("%2.4lf\t",Cp[j*mA+i]));
        }
        PRLEVEL (1, ("\n"));
    }
#endif 

    double alpha = -1;
    BLAS_INT lda = (BLAS_INT) rowCount;
    BLAS_INT ldb = (BLAS_INT) fp;
    BLAS_INT ldc = (BLAS_INT) (rowCount - fp);

    PRLEVEL (1, ("%% alpha =%lf  ", alpha));
    PRLEVEL (1, ("%% lda =%d  ", lda));
    PRLEVEL (1, ("%% ldb =%d  ", ldb));
    PRLEVEL (1, ("%% ldc =%d\n", ldc));


    double beta= 0; //U part is not initialized

    BLAS_DGEMM ("N" ,"N" , &mA, &nB, &nA, &alpha, pF+fp, &lda, uPart, &ldb, 
            &beta, el, &ldc);





#ifndef NDEBUG 
    PRLEVEL (1, ("%%After DGEMM C =\n"));
    for(Int i = 0; i < mA; i++)
    {
        PRLEVEL (1, ("%% "));
        for(Int j = 0; j < nB; j++)
            PRLEVEL (1, ("%2.4lf\t",Cp[j*mA+i]));
        PRLEVEL (1, ("\n"));
    }
#endif 


    return 0;

}
