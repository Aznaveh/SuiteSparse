/** =========================================================================  /
 * =======================  paru_numeric_assemble ===========================  /
 * ========================================================================== */
/*  assemble numbers in columns, rows or fully  
 *  it is a gather scatter operation                                          */

#include "Parallel_LU.hpp"


void assemble_col (double *sC, double *dC,   //source col and destination col
                    Int m, Int *relRowInd)
{
    for (Int i = 0; i < m; i++) {
        if ( relRowInd[i] > 0 )  // If still valid
            dC [relRowInd[i] ] += sC[i];

    }
}

/*! TODO: design a different algorithm     */
void assemble_row (double *sR, double *dR,   //source Row and destination row
                    Int m, Int n, Int *relColInd)
{
    for (Int i = 0; i < m; i++) {
        if ( relColInd[i] > 0 )  // If still valid
            dR [relColInd[i] ] += sR[i];

    }
}
