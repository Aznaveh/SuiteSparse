/** =========================================================================  /
 * =======================  paru_numeric_assemble ===========================  /
 * ========================================================================== */
/*  assemble numbers in columns, rows or fully  
 *  it is a gather scatter operation                                          */

#include "Parallel_LU.hpp"

void assemble_col (const double *sC, double *dC,   //source col and destination col
                    Int m, const Int *relRowInd)
{
    DEBUGLEVEL (0);
    for (Int i = 0; i < m; i++) {
        PRLEVEL (1, ("%% relRowInd [%ld] =%ld\n",i ,relRowInd [i] ));
        Int ri = relRowInd[i] ;
        if ( ri >= 0 ) { // If still valid
            PRLEVEL (1, ("%% sC [%ld] =%2.5lf \n", i, sC [i]));
            PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", i, dC [ri]));
            dC [ri ] += sC[i];
            PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", i, dC [ri]));
        }
    }
    PRLEVEL (1, ("\n"));
}

void assemble_row (const double *sM, double *dM,//source and destination matrix 
                    Int sm, Int sn,    // dimension of source matrix
                    Int dm,     // dimension of destination matrix
                    Int sR, Int dR,     //source row and destination row
                    const Int *relColInd)
//Source and destination are stored column based
{
    DEBUGLEVEL (1);
    for (Int j = 0; j < sn; j++) {
        Int rj =relColInd[j] ;
        if ( rj >= 0 ){  // If still valid
            PRLEVEL (1, ("%% sM [%ld] =%2.5lf \n", sm*j+sR, sM [sm*j+sR] ));
            PRLEVEL (1, ("%% dM [%ld] =%2.5lf \n", rj*dm+dR, dM [rj*dm+dR]));
            dM [rj*dm + dR] += sM [sm*j + sR];
            PRLEVEL (1, ("%% dM [%ld] =%2.5lf \n", rj*dm+dR, dM [rj*dm+dR]));
        }
    }
}


void assemble_all (double *s, double *d,   //source and destination
        Int sm, Int sn,    // dimension of source matrix
        Int dm,    // dimension of destination matrix
        Int *relRowInd, Int *relColInd)
    //Source and destination are stored column based
{
    DEBUGLEVEL (0);
    for (Int j = 0; j < sn; j++) {
        Int rj =relColInd[j] ;
        if (rj  >= 0 ){  // If column is valid
            for (Int i = 0; i < sm; i++) {
                Int ri =relRowInd[i] ;
                if (ri >= 0 )  // If row is valid
                    d [rj*dm + ri ] += s[sm*j + i];

            }
        }
    }
}
