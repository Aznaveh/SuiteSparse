/** =========================================================================  /
 * =======================  paru_update_rel_ind =============================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"
/*! @brief updating element's relative indices in regard to another element 
 *      using my hash to find the columns and update relative indices
 *
 * @author Aznaveh
 * */


void paru_update_rel_ind_col ( Int e, Int f,  
        std::vector <Int> &colHash, 
        paru_matrix *paruMatInfo) 
{
    // updating relative column index 
    // it might be for curent element or for the Upart therefore we might even
    // dont have the curEl
    DEBUGLEVEL(0);
    PRLEVEL (1, ("%%update relative in %ld\n", f));

    paru_Element **elementList = paruMatInfo->elementList;
    paru_Element *el = elementList[e];

    //Int *el_Index = colIndex_pointer (el); //col global index of destination
    Int *el_Index = (Int*)(el+1); //col global index of destination

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    // Int *colRelIndex = relColInd (paru_Element *el);
    Int *colRelIndex = (Int*)(el+1) + mEl+ nEl;

    Int *fcolList = paruMatInfo->fcolList[f] ;

    for (Int i = el->lac; i < nEl; i++)
    {
        Int colInd = el_Index[i];
        if ( colInd < 0)
        {
            colRelIndex [i] = -1;
            continue;
        }
        PRLEVEL (1, ("%% searching for: cb_index[%ld]=%ld\n",
                    i,  colInd));
        Int found = paru_find_hash (colInd, colHash, fcolList);
        colRelIndex [i] = found;
        ASSERT (found != -1);
    }

    PRLEVEL (1, ("%%update relative in %ld finished\n", f));

    //update the cVal of el
    el->cValid = paruMatInfo->time_stamp[f];
}
