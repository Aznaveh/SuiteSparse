//------------------------------------------------------------------------------
// GB_free_memory: wrapper for free
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

// A wrapper for free.  If p is NULL on input, it is not freed.

// The memory is freed using the free() function pointer passed in to GrB_init,
// which is typically the ANSI C free function.  The free_pool is bypassed.

#include "GB_stuff.h"

GB_PUBLIC   // accessed by the MATLAB tests in GraphBLAS/Test only
void GB_free_memory         // free memory, bypassing the free_pool
(
    // input/output
    void **p,               // pointer to allocated block of memory to free
    // input
    size_t size_allocated   // # of bytes actually allocated
)
{

    if (p != NULL && (*p) != NULL)
    { 
        printf ("\nhard free %p %ld %ld\n", *p, size_allocated,
            GB_Global_memtable_size (*p)) ;
        ASSERT (size_allocated == GB_Global_memtable_size (*p)) ;

        SuiteSparse_free (*p) ;

        #ifndef NDEBUG
        GB_Global_memtable_remove (*p) ;
        #endif

        GB_Global_free_pool_dump (2) ; GB_Global_memtable_dump ( ) ;
        (*p) = NULL ;
    }
}

