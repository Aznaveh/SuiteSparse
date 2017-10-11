/** =========================================================================  /
 * =======================  paru_trsm =======================================  /
 * ========================================================================== */


#include "Parallel_LU.hpp"

//extern "C"  int dtrsm_(char *side, char *uplo, char *transa, char *diag, 
//                           int *m, int *n, double *alpha, double *a, int *lda, 
//                                                         double *b, int *ldb);
Int paru_trsm(double *F, Int m, Int n,
        int *ipiv, cholmod_common *cc)
{
 
}
