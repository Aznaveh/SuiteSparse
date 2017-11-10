/** =========================================================================  /
 * =======================  paru_numeric_assemble ===========================  /
 * ========================================================================== */
/*  assemble numbers in columns, rows or fully  
 *  it is a gather scatter operation                                          */

#include "Parallel_LU.hpp"

/*! TODO: index validity must be checked somehow	 */
void assemble_col (const double *sC, double *dC,   //source col and destination col
                    Int m, const Int *relRowInd)
{
    DEBUGLEVEL (1);
    for (Int i = 0; i < m; i++) {
        PRLEVEL (1, ("relRowInd [%ld] =%ld\n",i ,relRowInd [i] ));
        Int ri = relRowInd[i] ;
        if ( ri >= 0 ) { // If still valid
            PRLEVEL (1, ("sC [%ld] =%2.5lf\n", i, sC [i]));
            dC [ri ] += sC[i];
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
    DEBUGLEVEL (0);
    for (Int j = 0; j < sn; j++) {
        Int rj =relColInd[j] ;
        if ( rj >= 0 ){  // If still valid
            dM [rj*dm + dR] += sM [sm*j + sR];
            PRLEVEL (1, ("adding %lf\n", sM [sm*j + sR]));
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
