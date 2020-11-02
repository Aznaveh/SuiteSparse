/** =========================================================================  /
 * =======================  paru_bin_search =================================  /
 * ========================================================================== */
/*!
 * @brief binary search in different contexts
 *
 * @author Aznaveh
 *  */
#include "Parallel_LU.hpp"
#define LEN 32

Int bin_srch_ind  (Int *srt_lst, Int *ind_lst, Int l, Int r, Int num)
    // when the original list is not sorted and we also have to send the index
    // used for the row
{
    if ( r >= l) 
    {
        Int mid = l + (r-l)/2;
        if (srt_lst[mid] == num)
            return ind_lst[mid];

        if (srt_lst[mid] >  num)
            return bin_srch_ind (srt_lst, ind_lst, l, mid-1, num);
        return bin_srch_ind (srt_lst, ind_lst, mid+1, r,  num);
    }
    return (-1);
}

Int bin_srch (Int *srt_lst, Int l, Int r, Int num)
    // a simple binary search for when we know all the indices are available
    // TODO: add linear search for small cases
{
    if ( r >= l + LEN) 
    {
        Int mid = l + (r-l)/2;
        if (srt_lst[mid] == num)
            return mid;

        if (srt_lst[mid] >  num)
            return bin_srch_col (srt_lst, l, mid-1, num);
        return bin_srch_col (srt_lst, mid+1, r,  num);
    }
    if ( r >= l ) 
    {
        for (Int i = l; i < r; r++)
            if ( srt_lst [i] == num )
                return i;
    }
    return (-1);
}

Int bin_srch_col (Int *srt_lst, Int l, Int r, Int num)
    // a simple binary search for when it is possible that some columns were 
    // flipped
{
    if ( r >= l + LEN) 
    {
        Int mid = l + (r-l)/2;
        Int srt_lstMid = (srt_lst[mid]<0)? flip(srt_lst[mid]) : srt_lst[mid];
        if (srt_lstMid == num)
            return mid;

        if (srt_lstMid >  num)
            return bin_srch_col (srt_lst, l, mid-1, num);
        return bin_srch_col (srt_lst, mid+1, r,  num);
    }
    if ( r >= l) 
    {
        for (Int i = l; i < r; r++)
            if ( srt_lst [i] == num )
                return i;
    }
    return (-1);
}
