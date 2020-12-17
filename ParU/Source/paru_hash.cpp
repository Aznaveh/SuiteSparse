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
#define HASH_FACTOR 107 

void paru_insert_hash(Int key, Int value, std::vector<Int> &colHash)
{
    DEBUGLEVEL(1);

#ifndef NDEBUG  
    Int p = 1;
    PRLEVEL (p, ("%% Insert hash key=%ld value=%ld ", key, value ));
    PRLEVEL (p, ("size=%ld \n", colHash.size() ));
    PRLEVEL (p, ("%% before insertion"));
    for (auto i:colHash)
        PRLEVEL (p, (" %ld ", i));
    PRLEVEL (p, ("\n"));

#endif

//     old version sometimes doesn't insert in hash
//    Int loop_cnt = 0;
//    Int size = colHash.size()-1;
//    Int  index = key % size;  //hash function
//    PRLEVEL (p, ("index =%ld \n", index ));
//    while ( colHash[index] != -1  )
//    { //finding an empty spot
//        loop_cnt++;
//        if( loop_cnt > log2(size) )
//            return; //without even inserting inside the hash
//        index = (index+1) % size;
//    }
//    colHash[index] = value;

    // newer version
    Int hash_bits = colHash.size() - 2;
    Int index = (key * HASH_FACTOR) & hash_bits;
    Int loop_cnt = 0;
    while ( colHash [index] != -1 )
    { //finding an empty spot
        index = (index+1) & hash_bits;
        PRLEVEL (p, ("index =%ld colHash=%ld\n", index, colHash [index] ));
        //if( loop_cnt > log2(size) )
        //if( loop_cnt++ > hash_bits )
        //    return; //without even inserting inside the hash
        loop_cnt++;
        ASSERT (loop_cnt < hash_bits);
    }
    colHash[index] = value;

#ifndef NDEBUG  
    p = 1;
    PRLEVEL (p, ("%% hash_bits == %lx ", hash_bits));
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
    Int p = 1;
    PRLEVEL (p, ("%% find for hash key=%ld \n", key));
#endif
    Int hash_bits = colHash.size() - 2;
    //lookup table
    if (colHash.back() == -1 )
    {
        PRLEVEL (p, ("%% LOOKUP key =%ld colHash=%ld \n", key, colHash [key]));
        return colHash [key];
    }

//    old version
//    Int size = colHash.size()-1;
//    Int index = key % size;
//    Int value = colHash [index];
//    Int loop_cnt = 0;
//
//    while (  value == -1 || fcolList [value] != key  )
//    {
//        index = (index+1) % size;
//        PRLEVEL (p, ("%% index =%ld \n", index ));
//        value = colHash [index];
//        if( loop_cnt++ > log2(size) )
//        { // take a long time in the hash; 
//            //  guarantees that find takes at most log time
//            PRLEVEL (p, ("%% binary search for hash\n"));
//            value = bin_srch (fcolList, 0, size-1, key);
//            break;
//        }
//    }


    Int index = (key * HASH_FACTOR) & hash_bits;
    Int value = colHash [index];
    Int loop_cnt = 0;
    Int size = colHash.back();
    while ( value != -1 && fcolList [value] != key )
    {
        index = (index+1) & hash_bits;
        PRLEVEL (p, ("%% index =%ld \n", index ));
        value = colHash [index];
        if( loop_cnt++ > log2(hash_bits) )
        { // take a long time in the hash; 
            //  guarantees that find takes at most log time
            PRLEVEL (p, ("%% binary search for hash\n"));
            value = bin_srch (fcolList, 0, size-1, key);
            break;
        }
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
