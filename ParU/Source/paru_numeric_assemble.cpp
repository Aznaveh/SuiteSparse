/** =========================================================================  /
 * =======================  paru_numeric_assemble ===========================  /
 * ========================================================================== */
/*! @brief  assemble numbers in columns, rows or fully  
 *          it is a gather scatter operation
 *
 * @author Aznaveh
 * 
 */

#include "Parallel_LU.hpp"

void assemble_col (const double *sC,  //source col and destination col
        double *dC,
        Int m, const Int *relRowInd)
{
    DEBUGLEVEL (0);
    for (Int i = 0; i < m; i++) 
    {
        PRLEVEL (1, ("%% relRowInd [%ld] =%ld\n",i ,relRowInd [i] ));
        Int ri = relRowInd[i] ;
        if ( ri >= 0 ) 
        { // If still valid
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
    //TODO get some info from el not to start from first col
    DEBUGLEVEL (0);
    for (Int j = 0; j < sn; j++) 
    {
        Int rj =relColInd[j] ;
        if ( rj >= 0 )
        {  // If still valid
            PRLEVEL (1, ("%% sM [%ld] =%2.5lf \n", sm*j+sR, sM [sm*j+sR] ));
            PRLEVEL (1, ("%% dM [%ld] =%2.5lf \n", rj*dm+dR, dM [rj*dm+dR]));
            dM [rj*dm + dR] += sM [sm*j + sR];
            PRLEVEL (1, ("%% dM [%ld] =%2.5lf \n", rj*dm+dR, dM [rj*dm+dR]));
        }
    }
}


void assemble_row_toU (Int e, Int f, Int sR, Int dR, 
        std::vector <Int> colHash, 
        paru_matrix *paruMatInfo) 
{

    DEBUGLEVEL (0);

    paru_Element **elementList = paruMatInfo->elementList;
    paru_Element *el = elementList[e];

    if (el->cValid != paruMatInfo->time_stamp[f] )
        //if not updatated
        paru_update_rel_ind_col ( e, f, colHash, paruMatInfo) ;

    paru_fac *Us =  paruMatInfo->partial_Us;
    double *uPart = Us[f].p ; //uPart

    Int nEl = el->ncols;
    Int mEl = el->nrows;


    paru_fac *LUs =  paruMatInfo->partial_LUs;
    Int fp = LUs[f].n; 

    //Int *el_colIndex = colIndex_pointer (curEl);
    Int *el_colIndex = (Int*)(el+1);

    // Int *colRelIndex = relColInd (paru_Element *el);
    Int *colRelIndex = (Int*)(el+1) + mEl+ nEl;

    //double *el_Num = numeric_pointer (el);
    double *sM = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);

    Int ncolsSeen = el->ncolsleft;
    for (Int j = el->lac; j < nEl; j++) 
    {
        Int rj =colRelIndex[j] ;
        if ( el_colIndex[j] >= 0 )
        {  // If still valid
            ncolsSeen--;
            PRLEVEL (1, ("%% sM [%ld] =%2.5lf \n", mEl*j+sR, sM [mEl*j+sR] ));
            PRLEVEL (1, ("%% uPart [%ld] =%2.5lf \n", 
                        rj*fp+dR, uPart [rj*fp+dR]));
            uPart [rj*fp+ dR] += sM [mEl*j + sR];
            PRLEVEL (1, ("%% uPart [%ld] =%2.5lf \n", 
                        rj*fp+dR, uPart [rj*fp+dR]));
            if (ncolsSeen == 0) break;
        }
    }

}
