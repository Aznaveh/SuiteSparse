/** =========================================================================  /
 * =======================  paru_numeric_assemble ===========================  /
 * ========================================================================== */
/*  assemble numbers in columns, rows or fully                                */

#include "Parallel_LU.hpp"


void assemble_cols (double *sR, double *dR,   //source Row and destination row
                    Int m, Int *relRowInd)
{
    for (Int i = 0; i < m; i++) {
        if ( relRowInd[i] > 0 )  // If still valid
            dR [relRowInd[i] ] += sR[i];

    }
}
