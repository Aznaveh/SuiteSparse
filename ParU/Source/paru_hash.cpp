////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_hash  /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
#include "paru_internal.hpp"
// key*257 & mask
#define HASH_FUNCTION(key) (((key << 8) + (key)) & (hash_bits))

void paru_insert_hash(Int key, Int value, std::vector<Int> &colHash)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;

#ifndef NDEBUG
    PRLEVEL(PR, ("%% Insert hash key=%ld value=%ld ", key, value));
    PRLEVEL(PR, ("size=%ld \n", colHash.size()));
    PRLEVEL(PR, ("%% before insertion"));
    for (auto i : colHash) PRLEVEL(PR, (" %ld ", i));
    PRLEVEL(PR, ("\n"));

#endif

    Int hash_bits = colHash.size() - 2;
    Int index = HASH_FUNCTION(key);

    Int loop_cnt = 0;
    while (colHash[index] != -1)
    {  // finding an empty spot
        index = (index + 1) & hash_bits;
        PRLEVEL(PR, ("index =%ld colHash=%ld\n", index, colHash[index]));
        loop_cnt++;
        ASSERT(loop_cnt < hash_bits);
    }
    colHash[index] = value;

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%% hash_bits == %lx ", hash_bits));
    PRLEVEL(PR, ("%%"));
    for (auto i : colHash) PRLEVEL(PR, (" %ld ", i));
    PRLEVEL(PR, ("\n"));
#endif
}

Int paru_find_hash(Int key, std::vector<Int> &colHash, Int *fcolList)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;
#ifndef NDEBUG
    PRLEVEL(PR, ("%% find for hash key=%ld \n", key));
#endif
    // lookup table
    if (colHash.back() == -1)
    {
        PRLEVEL(PR, ("%% LOOKUP key =%ld colHash=%ld \n", key, colHash[key]));
        return colHash[key];
    }

    Int hash_bits = colHash.size() - 2;
    Int index = HASH_FUNCTION(key);
    Int value = colHash[index];
    Int loop_cnt = 0;
    Int size = colHash.back();
    while (value != -1 && fcolList[value] != key)
    {
        index = (index + 1) & hash_bits;
        PRLEVEL(PR, ("%% index =%ld \n", index));
        value = colHash[index];
        if (loop_cnt++ > log2(hash_bits))
        {  // take a long time in the hash;
            //  guarantees that find takes at most log time
            PRLEVEL(PR, ("%% binary search for hash\n"));
            value = bin_srch(fcolList, 0, size - 1, key);
            break;
        }
    }

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%%"));
    for (auto i : colHash) PRLEVEL(PR, (" %ld ", i));
    PRLEVEL(PR, ("\n"));
    PRLEVEL(PR, ("%% value is =%ld \n", value));
    Int bsRes = bin_srch(fcolList, 0, size - 1, key);
    PRLEVEL(PR, ("%% binSearch=%ld \n", bsRes));
#endif
    return value;
}
