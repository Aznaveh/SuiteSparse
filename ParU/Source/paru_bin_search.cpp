/** =========================================================================  /
 * =======================  paru_bin_search =================================  /
 * ========================================================================== */
/*!
 * @brief binary search in different contexts
 *
 * @author Aznaveh
 *  */
#include "Parallel_LU.hpp"
Int bin_srch  (Int *srt_lst, Int *ind_lst, Int l, Int r, Int num)
{
    if ( r >= l) 
    {
        Int mid = l + (r-l)/2;
        if (srt_lst[mid] == num)
            return ind_lst[mid];

        if (srt_lst[mid] >  num)
            return bin_srch (srt_lst, ind_lst, l, mid-1, num);
        return bin_srch (srt_lst, ind_lst, mid+1, r,  num);
    }
    return (-1);
}

Int bin_srch_col (Int *srt_lst, Int l, Int r, Int num)
{
    if ( r >= l) 
    {
        Int mid = l + (r-l)/2;
        if (srt_lst[mid] == num)
            return mid;

        if (srt_lst[mid] >  num)
            return bin_srch_col (srt_lst, l, mid-1, num);
        return bin_srch_col (srt_lst, mid+1, r,  num);
    }
    return (-1);
}
