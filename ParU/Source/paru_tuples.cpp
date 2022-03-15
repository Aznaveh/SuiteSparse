////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_tuples ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief adding tuples to row and column list
 *
 * Row and column tuple list add is basically the same in sequential case
 * I have two functions for my parallel algorithm
 *
 * @param  a list of tuples and the the tuple we want to add
 * @return 0 on sucess
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
Int paru_add_rowTuple(ParU_TupleList *RowList, Int row, ParU_Tuple T)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("row =%ld, (%ld,%ld)\n", row, T.e, T.f));

    ParU_TupleList *cur = &RowList[row];

    if (cur->len > cur->numTuple)
        cur->list[cur->numTuple++] = T;

    else
    {
        PRLEVEL(1, ("%%Row ParU_Tuple reallocating space\n"));
        Int newLen = cur->len * 2 + 1;
        ParU_Tuple *newList =
            (ParU_Tuple *)paru_alloc(newLen, sizeof(ParU_Tuple));
        if (newList == NULL)  // Error in allocating memory
            return 1;
        for (Int i = 0; i < cur->numTuple; ++i)  // copy old to new
            newList[i] = cur->list[i];
        paru_free(cur->len, sizeof(ParU_Tuple), cur->list);
        cur->len = newLen;
        cur->list = newList;
        cur->list[cur->numTuple++] = T;
    }
    return 0;
}
