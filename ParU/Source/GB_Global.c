//------------------------------------------------------------------------------
// GB_Global: global values in GraphBLAS
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

#include "GB_stuff.h"

//------------------------------------------------------------------------------
// Global storage: for all threads in a user application that uses GraphBLAS
//------------------------------------------------------------------------------

typedef struct
{

    //--------------------------------------------------------------------------
    // for malloc debugging only
    //--------------------------------------------------------------------------

    #ifndef NDEBUG
    #define GB_MEMTABLE_SIZE 10000
    void    *memtable_p [GB_MEMTABLE_SIZE] ;
    size_t   memtable_s [GB_MEMTABLE_SIZE] ;
    #endif
    int nmemtable ;

    //--------------------------------------------------------------------------
    // internal memory pool
    //--------------------------------------------------------------------------

    // free_pool [k] is a pointer to a link list of freed blocks, all of size
    // exactly equal to 2^k.  The total number of blocks in the kth pool is
    // given by free_pool_nblocks [k], and the upper bound on this is given by
    // free_pool_limit [k].  If any additional blocks of size 2^k above that
    // limit are freed by GB_dealloc_memory, they are not placed in the pool,
    // but actually freed instead.

    void *free_pool [64] ;
    int64_t free_pool_nblocks [64] ;
    int64_t free_pool_limit [64] ;

}
GB_Global_struct ;

GB_PUBLIC GB_Global_struct GB_Global ;

GB_Global_struct GB_Global =
{

    // for malloc debugging only
    .nmemtable = 0,     // memtable is empty

    // all free_pool lists start out empty
    .free_pool = {
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },

    .free_pool_nblocks = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },

    // default limits on the number of free blocks in each list:
    .free_pool_limit = {
        0,      // size 2^0 = 1 byte   none
        0,      // size 2^1 = 2        none
        0,      // size 2^2 = 4        none

        16483,  // size 2^3 = 8        (2^14 blocks * 2^3  = 128 KB total)
        16483,  // size 2^4 = 16 bytes (2^14 blocks * 2^4  = 256 KB total)
        16483,  // size 2^5 = 32       (2^14 blocks * 2^5  = 512 KB total)
        16483,  // size 2^6 = 64       (2^14 blocks * 2^6  = 1 MB total)
        16483,  // size 2^7 = 128      (2^14 blocks * 2^7  = 2 MB total)

        16483,  // size 2^8 = 256      (2^14 blocks * 2^8  = 4 MB total)
        8192,   // size 2^9 = 512      (2^13 blocks * 2^9  = 4 MB total)
        4096,   // size 2^10 = 1 KB    (2^12 blocks * 2^10 = 4 MB total)
        2048,   // size 2^11 = 2 KB    (2^11 blocks * 2^11 = 4 MB total)

        1024,   // size 2^12 = 4 KB    (2^10 blocks * 2^12 = 4 MB total)
        512,    // size 2^13 = 8 KB    (2^9  blocks * 2^13 = 4 MB total)
        256,    // size 2^14 = 16 KB   (2^8  blocks * 2^14 = 4 MB total)
        128,    // size 2^15 = 32 KB   (2^7  blocks * 2^15 = 4 MB total)

        64,     // size 2^16 = 64 KB   (2^6  blocks * 2^16 = 4 MB total)
        32,     // size 2^17 = 128 KB  (2^5  blocks * 2^17 = 4 MB total)
        16,     // size 2^18 = 256 KB  (2^4  blocks * 2^18 = 4 MB total)
        8,      // size 2^19 = 512 KB  (2^3  blocks * 2^19 = 4 MB total)

        // maximum total size = about 52 MB
        // by default, no blocks larger than 512 KB are kept in the free_pool

        //Aznaveh: increasing the default
        0,      // size 2^20 = 1 MB
        0,      // size 2^21
        0,      // size 2^22
        0,      // size 2^23
        0,      // size 2^24
        0,      // size 2^25
        0,      // size 2^26
        0,      // size 2^27
        0,      // size 2^28
        0,      // size 2^29

        0,      // size 2^30 (1 GB)
        0,      // size 2^31
        0,      // size 2^32
        0,      // size 2^33
        0,      // size 2^34
        0,      // size 2^35
        0,      // size 2^36
        0,      // size 2^37
        0,      // size 2^38
        0,      // size 2^39

        // These larger sizes are of course unlikely to appear, but adding all
        // 64 possibilities means that the free_pool does not need to check an
        // upper bound.

        0,      // size 2^40 (1 TB)
        0,      // size 2^41
        0,      // size 2^42
        0,      // size 2^43
        0,      // size 2^44
        0,      // size 2^45
        0,      // size 2^46
        0,      // size 2^47
        0,      // size 2^48
        0,      // size 2^49

        0,      // size 2^50 (1 PB)
        0,      // size 2^51
        0,      // size 2^52
        0,      // size 2^53
        0,      // size 2^54
        0,      // size 2^55
        0,      // size 2^56
        0,      // size 2^57
        0,      // size 2^58
        0,      // size 2^59

        0,      // size 2^60 (1 exabyte)
        0,      // size 2^61
        0,      // size 2^62
        0 },    // size 2^63 (4 exabytes!)

} ;

//------------------------------------------------------------------------------
// malloc debuging
//------------------------------------------------------------------------------

// These functions keep a separate record of the pointers to all allocated
// blocks of memory and their sizes, just for sanity checks.

GB_PUBLIC
void GB_Global_memtable_dump (void)
{
    #ifndef NDEBUG
    printf ("\nmemtable dump: %d\n", GB_Global.nmemtable) ;
    for (int k = 0 ; k < GB_Global.nmemtable ; k++)
    {
        printf ("  %4d: %12p : %ld\n", k,
            GB_Global.memtable_p [k],
            GB_Global.memtable_s [k]) ;
    }
    #endif
}

GB_PUBLIC
int GB_Global_memtable_n (void)
{
    return (GB_Global.nmemtable) ;
}

GB_PUBLIC
void GB_Global_memtable_clear (void)
{
    GB_Global.nmemtable = 0 ;
}

// add a pointer to the table of malloc'd blocks
GB_PUBLIC
void GB_Global_memtable_add (void *p, size_t size)
{
    #ifndef NDEBUG
    ASSERT ((p == NULL) == (size == 0)) ;
    if (p == NULL) return ;
    bool fail = false ;
    // printf ("memtable add %p size %ld\n", p, size) ;
    #pragma omp critical(GB_memtable)
    {
        int n = GB_Global.nmemtable  ;
        fail = (n > GB_MEMTABLE_SIZE) ;
        if (!fail)
        {
            for (int i = 0 ; i < n ; i++)
            {
                if (p == GB_Global.memtable_p [i])
                {
                    printf ("\nadd duplicate %p size %ld\n", p, size) ;
                    GB_Global_memtable_dump ( ) ;
                    printf ("Hey %d %p\n", i,p) ;
                    fail = true ;
                    break ;
                }
            }
        }
        if (!fail && p != NULL)
        {
            GB_Global.memtable_p [n] = p ;
            GB_Global.memtable_s [n] = size ;
            GB_Global.nmemtable++ ;
        }
    }
    ASSERT (!fail) ;
    // GB_Global_memtable_dump ( ) ;
    #endif
}

// get the size of a malloc'd block
GB_PUBLIC
size_t GB_Global_memtable_size (void *p)
{
    size_t size = 0 ;
    #ifndef NDEBUG
    if (p == NULL) return (0) ;
    bool found = false ;
    #pragma omp critical(GB_memtable)
    {
        int n = GB_Global.nmemtable  ;
        for (int i = 0 ; i < n ; i++)
        {
            if (p == GB_Global.memtable_p [i])
            {
                size = GB_Global.memtable_s [i] ;
                found = true ;
                break ;
            }
        }
    }
    if (!found)
    {
        printf ("\nFAIL: %p not found\n", p) ;
        GB_Global_memtable_dump ( ) ;
        ASSERT (0) ;
    }
    #endif
    return (size) ;
}

// test if a malloc'd block is in the table
GB_PUBLIC
bool GB_Global_memtable_find (void *p)
{
    bool found = false ;
    #ifndef NDEBUG
    if (p == NULL) return (false) ;
    #pragma omp critical(GB_memtable)
    {
        int n = GB_Global.nmemtable  ;
        for (int i = 0 ; i < n ; i++)
        {
            if (p == GB_Global.memtable_p [i])
            {
                found = true ;
                break ;
            }
        }
    }
    #endif
    return (found) ;
}

// remove a pointer from the table of malloc'd blocks
GB_PUBLIC
void GB_Global_memtable_remove (void *p)
{
    #ifndef NDEBUG
    if (p == NULL) return ;
    bool found = false ;
    // printf ("memtable remove %p ", p) ;
    #pragma omp critical(GB_memtable)
    {
        int n = GB_Global.nmemtable  ;
        for (int i = 0 ; i < n ; i++)
        {
            if (p == GB_Global.memtable_p [i])
            {
                // found p in the table; remove it
                // printf ("size %ld\n", GB_Global.memtable_s [i]) ;
                GB_Global.memtable_p [i] = GB_Global.memtable_p [n-1] ;
                GB_Global.memtable_s [i] = GB_Global.memtable_s [n-1] ;
                GB_Global.nmemtable -- ;
                found = true ;
                break ;
            }
        }
    }
    if (!found)
    {
        printf ("remove %p NOT FOUND\n", p) ;
        GB_Global_memtable_dump ( ) ;
    }
    ASSERT (found) ;
    // GB_Global_memtable_dump ( ) ;
    #endif
}


//------------------------------------------------------------------------------
// free_pool: fast access to free memory blocks
//------------------------------------------------------------------------------

// each free block contains a pointer to the next free block.  This requires
// the free block to be at least 8 bytes in size.
#define GB_NEXT(p) ((void **) p) [0]

// free_pool_init: initialize the free_pool
GB_PUBLIC
void GB_Global_free_pool_init (bool clear)
{ 
    #pragma omp critical(GB_free_pool)
    {
        if (clear)
        {
            // clear the free pool
            for (int k = 0 ; k < 64 ; k++)
            {
                GB_Global.free_pool [k] = NULL ;
                GB_Global.free_pool_nblocks [k] = 0 ;
            }
        }

        int64_t n = 16384 ;
        for (int k = 0 ; k < 64 ; k++)
        {
            GB_Global.free_pool_limit [k] = n ;
        }
        #if 0
        // set the default free_pool_limit
        for (int k = 0 ; k < 64 ; k++)
        {
            GB_Global.free_pool_limit [k] = 0 ;
        }
        int64_t n = 16384 ;
        for (int k = 3 ; k <= 8 ; k++)
        {
            GB_Global.free_pool_limit [k] = n ;
        }
        for (int k = 9 ; k <= 19 ; k++)
        {
            n = n/2 ;
            GB_Global.free_pool_limit [k] = n ;
        }
        //Aznaveh
        for (int k = 20 ; k < 64 ; k++)
        {
            n = n/2 ;
            GB_Global.free_pool_limit [k] = n ;
        }
        #endif
    }
}

#ifndef NDEBUG
// check if a block is valid
static inline void GB_Global_free_pool_check (void *p, int k, char *where)
{
    // check the size of the block
    // printf ("check %p\n", p) ;
    ASSERT (k >= 3 && k < 64) ;
    ASSERT (p != NULL) ;
    size_t size = GB_Global_memtable_size (p) ;
    ASSERT (size == ((size_t) 1) << k) ;
}
#endif

// free_pool_get: get a block from the free_pool, or return NULL if none
    GB_PUBLIC
void *GB_Global_free_pool_get (int k)
{
    void *p = NULL ;
    ASSERT (k >= 3 && k < 64) ;
    #pragma omp critical(GB_free_pool)
    {
        p = GB_Global.free_pool [k] ;
        if (p != NULL)
        {
            // remove the block from the kth free_pool
            GB_Global.free_pool_nblocks [k]-- ;
            GB_Global.free_pool [k] = GB_NEXT (p) ;
        }
    }
    if (p != NULL)
    { 
        // clear the next pointer inside the block, since the block needs
        // to be all zero
        // printf ("got %p k %d\n", p, k) ;
        #ifndef NDEBUG
        GB_Global_free_pool_check (p, k, "get") ;
        #endif
        GB_Global_free_pool_dump (2) ; // printf ("\ndid get\n\n") ;
    }
    return (p) ;
}

// free_pool_put: put a block in the free_pool, unless it is full
    GB_PUBLIC
bool GB_Global_free_pool_put (void *p, int k)
{ 
    #ifndef NDEBUG
    GB_Global_free_pool_check (p, k, "put") ;
    #endif
    bool returned_to_pool = false ;
    #pragma omp critical(GB_free_pool)
    {
        returned_to_pool =
            (GB_Global.free_pool_nblocks [k] < GB_Global.free_pool_limit [k]) ;
        if (returned_to_pool)
        {
            // add the block to the head of the free_pool list
            // printf ("put %p k %d\n", p, k) ;
            GB_Global.free_pool_nblocks [k]++ ;
            GB_NEXT (p) = GB_Global.free_pool [k] ;
            GB_Global.free_pool [k] = p ;
        }
    }
    GB_Global_free_pool_dump (2) ; // printf ("\ndid put\n\n") ;
    return (returned_to_pool) ;
}

// free_pool_dump: check the validity of the free_pool
GB_PUBLIC
void GB_Global_free_pool_dump (int pr)
{
    #ifndef NDEBUG
    bool fail = false ;
    #pragma omp critical(GB_free_pool)
    {
        for (int k = 0 ; k < 64 && !fail ; k++)
        {
            int64_t nblocks = GB_Global.free_pool_nblocks [k] ;
            int64_t limit   = GB_Global.free_pool_limit [k] ;
            if (nblocks != 0 && pr > 0)
            {
                printf ("pool %2d: %8ld blocks, %8ld limit\n",
                    k, nblocks, limit) ;
            }
            int64_t nblocks_actual = 0 ;
            void *p = GB_Global.free_pool [k] ;
            for ( ; p != NULL && !fail ; p = GB_NEXT (p))
            {
                if (pr > 1) printf ("  %16p ", p) ;
                size_t size = GB_Global_memtable_size (p) ;
                if (pr > 1) printf ("size: %ld\n", size) ;
                nblocks_actual++ ;
                fail = fail || (size != ((size_t) 1) << k) ;
                if (fail && pr > 0) printf ("    fail\n") ;
                fail = fail || (nblocks_actual > nblocks) ;
            }
            if (nblocks_actual != nblocks)
            {
                if (pr > 0) printf ("fail: # blocks %ld %ld\n",
                    nblocks_actual, nblocks) ;
                fail = true ;
            }
        }
    }
    ASSERT (!fail) ;
    #endif
}

// free_pool_limit_get: get the limit on the # of blocks in the kth pool
GB_PUBLIC
int64_t GB_Global_free_pool_limit_get (int k)
{
    int64_t nblocks = 0 ;
    if (k >= 3 && k < 64)
    { 
        #pragma omp critical(GB_free_pool)
        {
            nblocks = GB_Global.free_pool_limit [k] ;
        }
    }
    return (nblocks) ;
}

// free_pool_limit_set: set the limit on the # of blocks in the kth pool
GB_PUBLIC
void GB_Global_free_pool_limit_set (int k, int64_t nblocks)
{
    if (k >= 3 && k < 64)
    { 
        #pragma omp critical(GB_free_pool)
        {
            GB_Global.free_pool_limit [k] = nblocks ;
        }
    }
}

// free_pool_nblocks_total:  total # of blocks in free_pool (for debug only)
GB_PUBLIC
int64_t GB_Global_free_pool_nblocks_total (void)
{
    int64_t nblocks = 0 ;
    #pragma omp critical(GB_free_pool)
    {
        for (int k = 0 ; k < 64 ; k++)
        {
            nblocks += GB_Global.free_pool_nblocks [k] ;
        }
    }
    return (nblocks) ;
}

