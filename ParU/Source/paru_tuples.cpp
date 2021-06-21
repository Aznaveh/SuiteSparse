#include "Parallel_LU.hpp"
/** =========================================================================  /
 * =======================  paru_tuples =====================================  /
 * ========================================================================== */
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
Int paru_add_rowTuple(tupleList *RowList, Int row, paru_Tuple T)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("row =%ld, (%ld,%ld)\n", row, T.e, T.f));

    tupleList *cur = &RowList[row];

    if (cur->len > cur->numTuple)
        cur->list[cur->numTuple++] = T;

    else
    {
        PRLEVEL(1, ("%%Row paru_Tuple reallocating space\n"));
        Int newLen = cur->len * 2 + 1;
        paru_Tuple *newList =
            (paru_Tuple *)paru_alloc(newLen, sizeof(paru_Tuple));
        if (newList == NULL)  // Error in allocating memory
            return 1;
        for (Int i = 0; i < cur->numTuple; ++i)  // copy old to new
            newList[i] = cur->list[i];
        paru_free(cur->len, sizeof(paru_Tuple), cur->list);
        cur->len = newLen;
        cur->list = newList;
        cur->list[cur->numTuple++] = T;
    }
    return 0;
}
