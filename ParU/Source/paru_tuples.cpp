#include "Parallel_LU.hpp"
/** =========================================================================  /
 * =======================  paru_tuples =====================================  /
 * ==========================================================================  /
 *! @brief adding tuples to row and column list
 *
 * Row and column tuple list add is basically the same in sequential case
 * I have two functions for my parallel algorithm
 *
 * @param  a list of tuples and the the tuple we want to add
 * @return 0 on sucess 
 *
 * @author Aznaveh
 * */
/*! TODO: Cap number of tuples to maximum number of elements; 
 *      How to do this? 
 *      The maximum number is a constant based on the input
 *      matrix. should I send it as an argument to the function or there is a
 *      better way of doing this*/
Int paru_add_rowTuple (
        tupleList *RowList, 
        Int row, 
        paru_Tuple T, 
        cholmod_common *cc)
{
    DEBUGLEVEL(0);
    PRLEVEL (1, ("row =%ld, (%ld,%ld)\n", row,T.e, T.f));

    tupleList *cur = &RowList[row];

    if (cur->len > cur->numTuple)
        cur->list[cur->numTuple++] = T;

    else{
        PRLEVEL (1, ("%%Row paru_Tuple reallocating space\n"));
        Int newLen = cur->len*2 + 1;
        paru_Tuple *newList = 
            (paru_Tuple*) paru_alloc (newLen, sizeof(paru_Tuple), cc);
        if (newList == NULL)    // Error in allocating memory
            return 1;
        for (Int i = 0; i < cur->numTuple; ++i) //copy old to new
            newList [i] = cur->list [i];
        paru_free (cur->len, sizeof(paru_Tuple), cur->list, cc); 
        cur->len = newLen;
        cur->list = newList;
        cur->list [cur->numTuple++] = T;
    }
    return 0;
}

Int paru_add_colTuple (tupleList *ColList, Int col, 
        paru_Tuple T, cholmod_common *cc)  
{
    DEBUGLEVEL(0);
    PRLEVEL (1, ("%% col=%ld, (%ld,%ld)\n", col, T.e, T.f));
    tupleList *cur = &ColList [col];
    PRLEVEL (1, ("%% cur->numTuple =%ld\n", cur->numTuple));
    if (cur->len > cur->numTuple)
        cur->list [cur->numTuple++] = T;
    else{
        PRLEVEL (1, ("%% Col paru_Tuple reallocating space\n"));
        Int newLen = cur->len*2 + 1;
        paru_Tuple *newList = 
            (paru_Tuple*) paru_alloc (newLen, sizeof(paru_Tuple), cc);
        if (newList == NULL)    // Error in allocating memory
            return 1;
        for (Int i = 0; i < cur->numTuple; ++i) //copy old to new
            newList [i] = cur->list [i];
        paru_free (cur->len, sizeof(paru_Tuple), cur->list, cc); 
        cur->len = newLen;
        cur->list = newList;
        cur->list [cur->numTuple++] = T;
    }
    return 0;
}

/*TODO I am not using remove functions so maybe delete them */
Int paru_remove_colTuple(tupleList *ColList, Int col, Int t){
    // pick the last tuple and insert into the one should be removed
    DEBUGLEVEL(0);
    tupleList *curColTupleList = &ColList[col];
    Int  numTuple = curColTupleList -> numTuple;
    PRLEVEL (1, ("%% cur->numTuple =%ld\n", numTuple));
    paru_Tuple *listColTuples = curColTupleList->list;
    listColTuples [t] = listColTuples [numTuple-1];
    curColTupleList -> numTuple--;
    return 0;
}

Int paru_remove_rowTuple(tupleList *RowList, Int row, Int t){
    // pick the last tuple and insert into the one should be removed
    DEBUGLEVEL(0);
    tupleList *curRowTupleList = &RowList[row];
    Int  numTuple = curRowTupleList -> numTuple;
    PRLEVEL (1, ("%% cur->numTuple =%ld\n", numTuple));
    paru_Tuple *listRowTuples = curRowTupleList->list;
    listRowTuples [t] = listRowTuples [numTuple-1];
    curRowTupleList -> numTuple--;
    return 0;
}
