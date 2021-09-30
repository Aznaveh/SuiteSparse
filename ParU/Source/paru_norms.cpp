////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_norms /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  computing norms: 1-norm for vectors and sparse matrix
 *  and matrix
 *  @author Aznaveh
 */
#include "paru_internal.hpp"
double paru_spm_1norm(cholmod_sparse *A)
{
    // 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum
    // CSparse
    DEBUGLEVEL(0);
    if (!(A) || !A->x) return (-1);
    Int n = A->ncol;
    Int *Ap = (Int *)A->p;
    double *Ax = (double *)A->x;

    double norm = 0;
    for (Int j = 0; j < n; j++)
    {
        double s = 0;
        for (Int p = Ap[j]; p < Ap[j + 1]; p++)
        {
            PRLEVEL(3, ("Ax[%ld] = %.2lf\n", p, Ax[p]));
            s += fabs(Ax[p]);
        }
        PRLEVEL(2, ("s = %.2lf\n", s));
        norm = MAX(norm, s);
    }
    PRLEVEL(1, ("norm = %.8lf\n", norm));
    return (norm);
}

double paru_vec_1norm(const double *x, Int n)
{
    DEBUGLEVEL(0);
    double norm = 0.0;
    for (Int i = 0; i < n; i++)
    {
        PRLEVEL(1, ("so far norm = %lf + %lf\n", norm, fabs(x[i])));
        norm += fabs(x[i]);
    }
    PRLEVEL(1, ("vec 1norm = %.8lf\n", norm));
    return (norm);
}
