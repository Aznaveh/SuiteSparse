/** =========================================================================  /
 * =======================  paru_hash  ======================================  /
 * ========================================================================== */

/*! @brief functions to deal with the hash, insert and find
 * 
 * insert: 
 *   gets key and value and put it into the already initailzed table with -1 
 *   simple linear probing 
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
    Int p = 1;
    PRLEVEL (p, ("%% Insert hash key=%ld value=%ld ", key, value ));
    PRLEVEL (p, ("size=%ld \n", colHash.size() ));
#endif
    Int loop_cnt = 0;

    Int size = colHash.size();
    Int  index = key % size;  //hash function
    PRLEVEL (p, ("index =%ld \n", index ));
    while ( colHash[index] != -1  )
    { //finding an empty spot
        loop_cnt++;
        if( loop_cnt > log2(size) )
            return; //without inserting inside the hash
        index = (index+1) % size;
    }

    colHash[index] = value;

#ifndef NDEBUG  
    p = 1;
    PRLEVEL (p, ("%%"));
    for (auto i:colHash)
        PRLEVEL (p, (" %ld ", i));
    PRLEVEL (p, ("\n"));
#endif
}

Int paru_find_hash (Int key, std::vector<Int> &colHash, Int *fcolList)
{
    DEBUGLEVEL(0);
#ifndef NDEBUG  
    Int p = 1;
    PRLEVEL (p, ("%% find for hash key=%ld \n", key));
#endif
    Int size = colHash.size();

    Int index = key % size;
    Int value = colHash [index];
    Int loop_cnt = 0;
    while (  fcolList [value] != key  )
    {
        PRLEVEL (p, ("%% index =%ld \n", index ));
        if( loop_cnt++ > log2(size) )
        { // take a long time in the hash; 
            //  guarantees that find takes at most log time
            PRLEVEL (p, ("%% binary search for hash\n"));
            value = bin_srch (fcolList, 0, size-1, key);
            break;
        }

        //++index %=  size;
        index = (index+1) % size;
        value = colHash [index];
    }
#ifndef NDEBUG  
    p = 1;
    PRLEVEL (p, ("%%"));
    for (auto i:colHash)
        PRLEVEL (p, (" %ld ", i));
    PRLEVEL (p, ("\n"));
    PRLEVEL (p, ("%% value is =%ld \n", value));
    Int bsRes = bin_srch (fcolList, 0, size-1, key);
    PRLEVEL (p, ("%% binSearch=%ld \n", bsRes));
    //ASSERT (bin_srch (fcolList, 0, size-1, key) == value );
#endif
    return value;
}
