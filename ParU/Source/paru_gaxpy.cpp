/** =========================================================================  /
 *  ======================  paru_gaxpy ======================================  /
 *  ========================================================================= */
/*! @brief  computing y = A*x+y
 *  doesn't need performance improvement mostly borrowed from CSparse
 *  just is used for tests
 *
 *
 *  @author Aznaveh
 */
#include "paru_internal.hpp"
Int paru_gaxpy (cholmod_sparse *A, const double *x, double *y)
{
    DEBUGLEVEL(0);

    if (!A || !x || !y) return (0);
    Int *Ap = (Int *)A->p;
    Int *Ai = (Int *)A->i;
    double *Ax = (double *)A->x;
    Int n = A->ncol;

    for (Int j = 0; j < n; j++)
    {
        for (Int p = Ap[j]; p < Ap[j + 1]; p++)
        {
            y[Ai[p]] += Ax[p] * x[j];
        }
    }

#ifndef NDEBUG
    Int m = A->nrow;
    PRLEVEL(1, ("%% after gaxpy y is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.8lf, ", y[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
