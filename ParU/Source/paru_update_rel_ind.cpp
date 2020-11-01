/** =========================================================================  /
 * =======================  paru_update_rel_in ==============================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"
/*! @brief updating element's relative indices in regard to another element 
 * @param paru_Element *el: current front which all the rleative indices are valid
 *        paru_Element *cb_el: contributing element that the relative indices need to
 *           be updated. There might be invalid rows/cols in contribution block
 *        char rc: it is either 'r' or 'c' row or column update
 *        cholmod_common *cc: Memory allcoation is needed here so cc must be
 *          informed about the amount of memory
 * @author Aznaveh
 * */

//TODO: I can romove this file totally
void inline swap_key(Int *srt_lst,Int *ind_lst, Int i, Int j )
{
    Int tmp = srt_lst[i]; srt_lst[i] = srt_lst[j]; srt_lst[j] = tmp;
    tmp = ind_lst[i]; ind_lst[i] = ind_lst[j]; ind_lst[j] = tmp;
}

Int partition (Int *srt_lst, Int *ind_lst, Int low, Int high)
{
    DEBUGLEVEL(0);
    Int mid= (low + high) /2 ; // pivot element 
    swap_key(srt_lst, ind_lst, mid, high);
    PRLEVEL (1, ("%% pivot = %ld\n",srt_lst[high])); 

    Int piv = srt_lst [high] ; // pivot element 

    Int j = low-1; 
    for(Int i = low; i < high; i++)
    {
        if (srt_lst[i] <= piv)
        {
            j++; 
            swap_key(srt_lst, ind_lst, i, j);
        }
    }

    swap_key(srt_lst, ind_lst,high , j+1);
    return (j+1);
}

void paru_qsort (Int *srt_lst, Int *ind_lst, Int low, Int high)
{ //recursive
    DEBUGLEVEL(0);
    PRLEVEL (1, ("%% low=%ld high=%ld  \n",low, high)); 
    if (low < high-15)
    {
        Int piv = partition (srt_lst, ind_lst, low, high);

        paru_qsort(srt_lst, ind_lst, low, piv-1);
        paru_qsort(srt_lst, ind_lst,  piv+1, high);
    }
    else
    {
        for (Int i=low; i <= high; i++)
            for (Int j=i; j <= high; j++)
            {
                if (srt_lst[i] > srt_lst [j])
                    swap_key(srt_lst, ind_lst, i, j);
            }


    }
}

void paru_sort (Int *srt_lst, Int *ind_lst, Int len)
{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    PRLEVEL (1, ("%% Before sort\n")); 
    for (Int i=0; i < len; i++)
    {  
        PRLEVEL (1, ("%% srt_lst[%ld]= %ld ind_lst[%ld]=%ld \n",
                    i, srt_lst[i], i, ind_lst[i]));
    }
#endif
    paru_qsort (srt_lst, ind_lst, 0, len-1);
    DEBUGLEVEL(0);
#ifndef NDEBUG
    PRLEVEL (1, ("%% After sort\n")); 
    for (Int i=0; i < len; i++)
    { 
        PRLEVEL (1, ("%% srt_lst[%ld]= %ld ind_lst[%ld]=%ld \n",
                    i, srt_lst[i], i, ind_lst[i]));
    }
#endif
}


void paru_update_rel_ind_col ( Int f, Int e, paru_matrix *paruMatInfo) 
{
    // updating relative column index 
    // it might be for curent element or for the Upart therefore we might even
    // dont have the curEl
    DEBUGLEVEL(1);
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
    paru_fac *Us =  paruMatInfo->partial_Us;
    Int colCount = Us[f].n;

    //TODO be sure not to need in assemble_row
        for (Int i = 0 ; i < el->lac; i++)
                colRelIndex [i] = -1;

    for (Int i = el->lac; i < nEl; i++)
    {
        Int global_ind = el_Index[i];
        if (global_ind < 0)
        {
            colRelIndex [i] = -1;
            continue;
        }
        PRLEVEL (1, ("%% searching for: cb_index[%ld]=%ld\n",
                    i,  global_ind));
        Int found = bin_srch (fcolList, 0, colCount, global_ind);
        colRelIndex [i] = found;
        ASSERT (found != -1);
    }

    PRLEVEL (1, ("%%update relative in %ld finished\n", f));
    //TODO: update the rVal of el
}
