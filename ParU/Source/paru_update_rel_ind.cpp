////////////////////////////////////////////////////////////////////////////////
/////////////////////////  paru_update_rel_ind ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*! @brief updating element's relative indices in regard to another element
 *      using my hash to find the columns and update relative indices
 *
 * @author Aznaveh
 * */

#include "paru_internal.hpp"
void paru_update_rel_ind_col(Int e, Int f, std::vector<Int> &colHash,
                               paru_work *Work, ParU_Numeric *Num)
{
    // updating relative column index
    // it might be for curent element or for the Upart therefore we might even
    // dont have the curEl
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%%update relative in %ld\n", f));

    paru_element **elementList = Work->elementList;
    paru_element *el = elementList[e];

    // Int *el_Index = colIndex_pointer (el); //col global index of destination
    Int *el_Index = (Int *)(el + 1);  // col global index of destination

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    // Int *colRelIndex = relColInd (paru_element *el);
    Int *colRelIndex = (Int *)(el + 1) + mEl + nEl;

    Int *fcolList = Num->fcolList[f];

    for (Int i = el->lac; i < nEl; i++)
    {
        Int colInd = el_Index[i];
        if (colInd < 0)
        {
            colRelIndex[i] = -1;
            continue;
        }
        PRLEVEL(1, ("%% searching for: cb_index[%ld]=%ld\n", i, colInd));
        Int found = paru_find_hash(colInd, colHash, fcolList);
        colRelIndex[i] = found;
        ASSERT(found != -1);
    }

    PRLEVEL(1, ("%%update relative in %ld finished\n", f));

    // update the cVal of el
    el->cValid = Work->time_stamp[f];
}
