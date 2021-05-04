//------------------------------------------------------------------------------
// GB_malloc_memory: wrapper for malloc
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

// A wrapper for malloc.  Space is not initialized.

#include "GB_stuff.h"

//------------------------------------------------------------------------------
// GB_malloc_helper:  use malloc to allocate an uninitialized memory block
//------------------------------------------------------------------------------

static inline void *GB_malloc_helper
(
    // input/output:
    size_t *size            // on input: # of bytes requested
                            // on output: # of bytes actually allocated
)
{
    void *p = NULL ;

    // determine the next higher power of 2
    size_t size_requested = (*size) ;
    (*size) = GB_IMAX (*size, 8) ;
    int k = GB_CEIL_LOG2 (*size) ;

    // if available, get the block from the pool
    if (GB_Global_free_pool_limit_get (k) > 0)
    {
        // round up the size to the nearest power of two
        (*size) = ((size_t) 1) << k ;
        p = GB_Global_free_pool_get (k) ;
        if (p != NULL) printf ("malloc from pool: %p %ld\n", p, *size) ;
    }

    if (p == NULL)
    {
        // no block in the free_pool, so allocate it
        p = SuiteSparse_malloc (*size, 1) ;
        #ifndef NDEBUG
        GB_Global_memtable_add (p, *size) ;
        #endif
        printf ("hard malloc %p %ld %ld\n", p, size_requested, *size) ;
    }
    GB_Global_free_pool_dump (2) ; GB_Global_memtable_dump ( ) ;

    return (p) ;
}

//------------------------------------------------------------------------------
// GB_malloc_memory
//------------------------------------------------------------------------------

GB_PUBLIC   // accessed by the MATLAB tests in GraphBLAS/Test only
void *GB_malloc_memory      // pointer to allocated block of memory
(
    size_t nitems,          // number of items to allocate
    size_t size_of_item,    // sizeof each item
    // output
    size_t *size_allocated  // # of bytes actually allocated
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (size_allocated != NULL) ;

    void *p ;
    size_t size ;

    // make sure at least one item is allocated
    nitems = GB_IMAX (1, nitems) ;

    // make sure at least one byte is allocated
    size_of_item = GB_IMAX (1, size_of_item) ;

    bool ok = GB_size_t_multiply (&size, nitems, size_of_item) ;
    if (!ok || nitems > GxB_INDEX_MAX || size_of_item > GxB_INDEX_MAX)
    { 
        // overflow
        (*size_allocated) = 0 ;
        return (NULL) ;
    }

    //--------------------------------------------------------------------------
    // allocate the memory block
    //--------------------------------------------------------------------------

    p = GB_malloc_helper (&size) ;

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    (*size_allocated) = (p == NULL) ? 0 : size ;
    ASSERT (GB_IMPLIES (p != NULL, size == GB_Global_memtable_size (p))) ;
    return (p) ;
}

