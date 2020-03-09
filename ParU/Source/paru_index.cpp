#include "Parallel_LU.hpp"
/** =========================================================================  /
 * =======================  paru_tuples =====================================  /
 * ==========================================================================  /
 *
 * @author Aznaveh
 * */
Int paru_add_index(paru_Index Ind, Int ind)
{
    Int size = paru_Index.size;
    Int count = paru_Index.count; 
    Int upperBound = paru_Index.upperBound;
    Int **listp = paru_Index.list;
    if (count == size)
    { // reallocate memeory and update size

    }
    Int *list = *listp;
    list[count++] = ind;
    return 0;
}
Int paru_cut_off_tail(paru_Index Ind)
{
    Int size = paru_Index.size;
    Int count = paru_Index.count; 
    Int upperBound = paru_Index.upperBound;
    Int **listp = paru_Index.list;

    Int *list = *listp;

    Int *relist;
    if (size != count)
    {
        Int sz = sizeof(Int)*size;
        relist = (Int*) paru_realloc (count, sizeof(Int), list, &sz, cc);
        if (relist == NULL)
        {

        }
        *list = relist;
    }

}
