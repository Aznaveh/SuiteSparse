/** =========================================================================  /
 * =======================  paru_solve ======================================  /
 * ==========================================================================  /
 * @brief  sovle Ax = b
 *      get a factorized matrix and a right hand side 
 *      returns x
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"

ParU_ResultCode paru_solve(paru_matrix *paruMatInfo, double *b)
// This routine call paru_front from first(f)...f including f
// This routine is called recursively to make tasks
{
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int m = LUsym->m;
    // TODO:  this goes to the test
    //if (LUsym->Pfin != NULL ) paru_perm(paruMatInfo);
    //
    //double b[m];
    //double x[m];
    //double xt[m];
    //for (Int i = 0; i < m; ++i) b[i] = i + 1;
      
    double *x = (double *)paru_alloc(m, sizeof(double));
    if ( x == NULL)
    {
        printf("Memory problem inside solve\n");
        return PARU_OUT_OF_MEMORY;
    }
      
      
    paru_apply_perm(LUsym->Pfin, b, x, m);  // x = b (p)
    if (paruMatInfo->scale_row) // x = s.Ps[x]
        paru_apply_scale (paruMatInfo->scale_row, LUsym->Ps, x, m, LUsym->n1);

    paru_lsolve(paruMatInfo, x);    // x = L\x 
    paru_usolve(paruMatInfo, x);    // x = U\x 
    paru_apply_inv_perm(LUsym->Qfill, x, b, m);  // b(q) = x


    // This can go to test
    //printf ("x = [ ");
    //for (Int i = 0; i < MIN(m,10); ++i) 
    //    printf ("%lf ",xt[i]);
    //printf (" ...]\n");

    //for (Int i = 0; i < m; ++i) b[i] *= -1;
    ////    b[i] = 0;
    //paru_gaxpy(A, xt, b);
    //double res = paru_vec_1norm(b, m);
    //double weighted_res = res / (paru_spm_1norm(A) * paru_vec_1norm(xt, m));
    //printf("Residual is |%.2lf| and weigheted residual is |%.2f|.\n", 
    //        log10(res), log10(weighted_res) );

}
