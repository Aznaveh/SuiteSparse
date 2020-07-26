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


void paru_update_rel_ind_row (paru_Element *el, paru_Element *cb_el, 
        cholmod_common *cc) 
{
    DEBUGLEVEL(0);

    //row relative index update
    //Int *el_Index = rowIndex_pointer (el);
    Int nEl = el->ncols;
    Int len_el = el->nrows;
    Int *el_Index = (Int*)(el+1)+nEl; // row global index of destination
    //Int *cb_el_Index = rowIndex_pointer (cb_el);
    Int nCbEl = cb_el->ncols;
    Int len_cb = cb_el->nrows;
    Int *cb_el_Index = (Int*)(cb_el+1)+nCbEl; // row global index of source
    // relative rowcol index of source to be updated
    //Int *RelIndex = relRowInd (cb_el);
    Int *RelIndex = (Int*)(cb_el+1)+ 2*nCbEl + len_cb; 
    Int *Work = el->rWork;


    if (len_cb*len_el < (len_cb+len_el)*log2(len_el) && Work == NULL )
    {
        //do a linear search if not sorted already and it is not  worth it
        PRLEVEL (1, ("%% Linear search\n"));
        for (Int i=0; i < len_cb ; i++)
        {
            Int global_ind = cb_el_Index[i];
            if (global_ind < 0) continue;
            PRLEVEL (1, ("%% searching for: cb_index[%ld]=%ld\n",
                        i,  global_ind));
#ifndef NDEBUG
            Int found = -1;
#endif
            //linear search for global_ind in contribution block
            for(Int i2 =0; i2 < len_el ; i2++)
                if (global_ind == el_Index[i2])
                {
                    RelIndex[i]=i2;
#ifndef NDEBUG
                    PRLEVEL (1, ("%%\t found in %ld\n", i2));
                    found = 1;
#endif
                    continue;
                }
#ifndef NDEBUG
            if(found != 1)
                PRLEVEL (0, ("%%t Not found in %ld \n",global_ind));
            ASSERT (found == 1);

#endif
        }
    } 
    else
    {


        PRLEVEL (1, ("%% Sort and Binary search\n"));
        if (Work == NULL )
        {  //Not sorted so far
            Work= (Int*) paru_alloc (2*len_el, sizeof(Int), cc);
            el->rWork = Work;

            Int *srt_lst = Work;      // list of relative indices; keys of sort   
            Int *ind_lst = Work + len_el; //list of indices 

            for (Int i=0; i < len_el ; i++)
            {  //initialize the lists
                srt_lst[i] = el_Index [i];
                ind_lst[i] = i;
            }
            paru_sort (srt_lst, ind_lst, len_el);
        }

        ASSERT (Work != NULL );
        Int *srt_lst = Work;      // list of relative indices; keys of sort   
        Int *ind_lst = Work + len_el; //list of indices 

        for (Int i=0; i < len_cb ; i++)
        {
            Int global_ind = cb_el_Index[i];
            if (global_ind < 0) continue;
            PRLEVEL (1, ("%% searching for: cb_index[%ld]=%ld\n",
                        i,  global_ind));
            Int found = bin_srch_ind 
                (srt_lst, ind_lst, 0, len_el-1, global_ind);
            RelIndex [i] = found;
            ASSERT (found != -1);
        }
    }
}

void paru_update_rel_ind_col (paru_matrix *paruMatInfo, Int f, 
        paru_Element *el, paru_Element *cb_el, cholmod_common *cc) 
{
    DEBUGLEVEL(0);

    //Int *el_Index = colIndex_pointer (el); //col global index of destination
    Int *el_Index = (Int*)(el+1); //col global index of destination
    //Int *cb_el_Index = colIndex_pointer (cb_el); //col global index of source
    Int *cb_el_Index = (Int*)(cb_el+1); //col global index of source
    Int len_cb = cb_el->ncols;
    Int mCbEl= cb_el->nrows;
    // relative col index of source to be updated
    //Int *RelIndex = relColInd (cb_el); 
    Int *RelIndex = (Int*)(cb_el+1)+ len_cb + mCbEl;
    Int len_el = el->ncols;
    

    Int *fcolList = paruMatInfo->fcolList[f] ;

    for (Int i=0; i < len_cb ; i++)
    {
        Int global_ind = cb_el_Index[i];
        if (global_ind < 0) continue;
        PRLEVEL (1, ("%% searching for: cb_index[%ld]=%ld\n",
                    i,  global_ind));
        Int found = bin_srch (fcolList, 0, len_el-1, global_ind);
        RelIndex [i] = found;
        ASSERT (found != -1);
    }
}
