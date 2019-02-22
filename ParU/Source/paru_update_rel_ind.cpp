/** =========================================================================  /
 * =======================  paru_update_rel_in ==============================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"
/*! @brief updating element's relative indices in regard to another element 
 * @param 
 */

void paru_update_rel_ind (Element *el, Element *cb_el, char rc) 
{
    DEBUGLEVEL(1);
        Int *el_Index, //row/col global index of destination
         *cb_el_Index,  //row/col global index of source
         *RelIndex,    //relative row/col index of source to be updated
            len_cb, len_el;       // number of rows/colum of el and cb_el

    if (rc == 'r'){ //row relative index update
        el_Index = rowIndex_pointer (el);
        cb_el_Index = rowIndex_pointer (cb_el);
        RelIndex = relRowInd (cb_el);
        len_cb = cb_el->nrows;
        len_el = el->nrows;
    }
    else{
        ASSERT (rc == 'c');
        el_Index = colIndex_pointer (el);
        cb_el_Index = colIndex_pointer (cb_el);
        RelIndex = relColInd (cb_el);
        len_cb = cb_el->ncols;
        len_el = el->ncols;
    }


    for (Int i=0; i < len_cb ; i++){
        Int global_ind = cb_el_Index[i];
        if (global_ind < 0) continue;
        PRLEVEL (1, ("%% searching for: cb_index[%ld]=%ld\n", i,  global_ind));
#ifndef NDEBUG
        Int found = -1;
#endif
        //search for global_ind in contribution block
        for(Int i2 =0; i2 < len_el ; i2++)
            if (global_ind == el_Index[i2]){
                RelIndex[i]=i2;
#ifndef NDEBUG
                PRLEVEL (1, ("%%\t found in %ld\n", i2));
                found = 1;
#endif
                continue;
            }
#ifndef NDEBUG
        ASSERT (found == 1);
        if(found != 1)
                PRLEVEL (1, ("%%\t %ld Not found in \n",global_ind));

#endif
    }
}

