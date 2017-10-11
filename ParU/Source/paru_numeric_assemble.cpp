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

void assemble_row (double *sM, double *dM,//source Matrix and destination matrix 
                    Int sm, Int sn,    // dimension of source matrix
                    Int dm,     // dimension of destination matrix
                    Int sR, Int dR,     //source row and destination row
                    Int *relColInd)
//Source and destination are stored column based
{
    for (Int i = 0; i < sn; i++) {
        if ( relColInd[i] > 0 ){  // If still valid
            dM [relColInd [i]*dm + dR] += sM [sm*i + sR];
            PRLEVEL (0, ("adding %lf\n", sM [sm*i + sR]));
        }
    }
}


void assemble_all (double *s, double *d,   //source and destination
                    Int sm, Int sn,    // dimension of source matrix
                    Int dm,    // dimension of destination matrix
                    Int *relRowInd, Int *relColInd)
//Source and destination are stored column based
{
    for (Int i = 0; i < sn; i++) {
        if ( relColInd[i] > 0 ){  // If column is valid
            for (Int j = 0; j < sm; j++) {
                if ( relRowInd[j] > 0 )  // If row is valid
                    d [relColInd[i]*dm + relRowInd[j] ] += s[sm*i + j];

            }
        }
    }
}
