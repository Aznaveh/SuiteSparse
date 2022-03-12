////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_bin_search ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief binary search in different contexts.
 *  IMPORTANT: it includes r.
 *
 * @author Aznaveh
 *  */
#include "paru_internal.hpp"
#define LEN 8 //XXX

Int bin_srch(Int *srt_lst, Int l, Int r, Int num)
// a simple binary search for when we know all the indices are available
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% BINSearch %ld,%ld for %ld\n", l, r, num));
    if (r >= l + LEN)
    {
        Int mid = l + (r - l) / 2;
        PRLEVEL(1, ("%% mid is %ld\n", mid));
        if (srt_lst[mid] == num) return mid;

        if (srt_lst[mid] > num)
        {
            PRLEVEL(1, ("%% 1 New %ld,%ld \n", l, mid - 1));
            return bin_srch(srt_lst, l, mid - 1, num);
        }
        PRLEVEL(1, ("%% 2 New %ld,%ld \n", mid + 1, r));
        return bin_srch(srt_lst, mid + 1, r, num);
    }

    if (r >= l)
    {
        for (Int i = l; i <= r; i++)
            if (srt_lst[i] == num) return i;
    }

    return (-1);
}

Int bin_srch_col(Int *srt_lst, Int l, Int r, Int num)
// a simple binary search for when it is possible that some columns were
// flipped
{
    if (r >= l + LEN)
    {
        Int mid = l + (r - l) / 2;
        Int srt_lstMid = (srt_lst[mid] < 0) ? flip(srt_lst[mid]) : srt_lst[mid];
        if (srt_lstMid == num) return mid;

        if (srt_lstMid > num) return bin_srch_col(srt_lst, l, mid - 1, num);
        return bin_srch_col(srt_lst, mid + 1, r, num);
    }

    if (r >= l)
    {
        for (Int i = l; i <= r; i++)
            if (srt_lst[i] == num) return i;
    }
    return (-1);
}
