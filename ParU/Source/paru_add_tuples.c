#include "Parallel_LU.hpp"
/** =========================================================================  /
 * =======================  paru_add_tuples =================================  /
 * ==========================================================================  /
/*! \brief adding tuples to row and column list
 *
 * Row and column tuple list add is basically the same in sequential case
 * I have two functions for my parallel algorithm
 *
 * \param  a list of tuples and the the tuple we want to add
 * \return 0 on sucess 
 */
int paru_add_rowTuple (tupleList *RowList, Int row, 
        Tuple T, cholmod_common *cc)
{
    DEBUGLEVEL(0);
    PRLEVEL (1, ("row =%ld, (%ld,%ld)\n", row,T.e, T.f));

    tupleList *cur = &RowList[row];

    if (cur->len > cur->numTuple)
        cur->list[cur->numTuple++] = T;

    else{
        PRLEVEL (0, ("Row Tuple reallocating space\n"));
        Int newLen = cur->len*2 + 1;
        Tuple *newList = 
            (Tuple*) paru_alloc (newLen, sizeof(Tuple), cc);
        if (newList == NULL)    // Error in allocating memory
            return 1;
        for (int i = 0; i < cur->len; ++i) //copy old to new
            newList [i] = cur->list [i];
        paru_free (cur->len, sizeof(Tuple), cur->list, cc); 
        cur->len = newLen;
        cur->list = newList;
        cur->list [cur->numTuple++] = T;
    }
    return 0;
}

int paru_add_colTuple (tupleList *ColList, Int col, 
        Tuple T, cholmod_common *cc)  /*! TODO: Sort for parallel case	 */
{
    DEBUGLEVEL(1);
    tupleList *cur = &ColList [col];
    if (cur->len > cur->numTuple)
        cur->list [cur->numTuple++] = T;
    else{
        PRLEVEL (1, ("Col Tuple reallocating space\n"));
        Int newLen = cur->len*2 + 1;
        Tuple *newList = 
            (Tuple*) paru_alloc (newLen, sizeof(Tuple), cc);
        if (newList == NULL)    // Error in allocating memory
            return 1;
        for (int i = 0; i < cur->len; ++i) //copy old to new
            newList [i] = cur->list [i];
        paru_free (cur->len, sizeof(Tuple), cur->list, cc); 
        cur->len = newLen;
        cur->list = newList;
        cur->list [cur->numTuple++] = T;
    }
    return 0;
}
