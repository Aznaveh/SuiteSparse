#include "Parallel_LU.hpp"
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
    tupleList cur= RowList[row];
    if (cur.len > cur.numTuple)
        cur.list [++cur.numTuple] = T;
    else{
            Int newLen = (cur.len*2)+1;
            Tuple *newList = 
                (Tuple*) paru_alloc (1, newLen*sizeof(Tuple), cc);
            if( newList == NULL)
                return 1;
            paru_free (cur.len, sizeof(Tuple), cur.list, cc); 
            for (int i = 0; i < cur.len; ++i) //copy old to new
                newList [i] = cur.list [i];
            cur.len = newLen;
            cur.list = newList;
            cur.list [++cur.numTuple] = T;
    }
    return 0;
}
   
int paru_add_colTuple (tupleList *ColList, Int col, 
        Tuple T, cholmod_common *cc)
{
    tupleList cur= ColList [col];
    if (cur.len > cur.numTuple)
        cur.list [++cur.numTuple] = T;
    else{
            Int newLen = (cur.len*2)+1;
            Tuple *newList = 
                (Tuple*) paru_alloc (1, newLen*sizeof(Tuple), cc);
            if( newList == NULL)
                return 1;
            paru_free (cur.len, sizeof(Tuple), cur.list, cc); 
            for (int i = 0; i < cur.len; ++i) //copy old to new
                newList [i] = cur.list [i];
            cur.len = newLen;
            cur.list = newList;
            cur.list [++cur.numTuple] = T;
    }
    return 0;
}
