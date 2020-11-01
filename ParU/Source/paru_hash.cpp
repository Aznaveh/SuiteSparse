/** =========================================================================  /
 * =======================  paru_hash  ======================================  /
 * ========================================================================== */

/*! @brief functions to deal with the hash, insert and find
 * 
 * insert: 
 *   gets key and value and put it into the already initailzed table with -1 
 * find:
 *   find for the key and check if it has the correct value
 *
 *
 *  @author Aznaveh
 */
#include "Parallel_LU.hpp"

void paru_insert_hash(Int key, Int value, std::vector<Int> &colHash)
{
    DEBUGLEVEL(0);

#ifndef NDEBUG  
    Int p = 0;
    Int loop_cnt = 0;
    PRLEVEL (p, ("%% Insert hash key=%ld value=%ld ", key, value ));
    PRLEVEL (p, ("size=%ld \n", colHash.size() ));
#endif

    Int size = colHash.size();
    Int  index = key % size;  //hash function
    PRLEVEL (p, ("index =%ld \n", index ));
    while ( colHash[index] != -1 )
    {
#ifndef NDEBUG  
    loop_cnt++;
#endif
        index = (index+1) % size;
    }
    ASSERT (loop_cnt <= size);
    colHash[index] = value;

#ifndef NDEBUG  
    p = 2;
    PRLEVEL (p, ("%%"));
    for (auto i:colHash)
        PRLEVEL (p, (" %ld ", i));
    PRLEVEL (p, ("\n"));
#endif
}

Int paru_find_hash (Int key, std::vector<Int> &colHash, Int *fcolList)
{
    DEBUGLEVEL(1);
#ifndef NDEBUG  
    Int p = 0;
#endif
    Int size = colHash.size();
    Int index = key % size;
    Int value = colHash [index];
    Int loop_cnt = 0;
    while (  value != -1 && fcolList [value] != key && loop_cnt <= size )
    {
        loop_cnt++;
        index = (index+1) % size;
        value = colHash [index];
    }
    return value;
}

