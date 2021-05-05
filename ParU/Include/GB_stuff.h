
#ifndef GB_STUFF_H
#define GB_STUFF_H

// force debugging off
#ifndef NDEBUG
#define NDEBUG
#endif

// turn on debugging
//#undef NDEBUG

// defined somewhere else
#ifdef ASSERT
#undef ASSERT
#endif
#ifndef NDEBUG
#include <assert.h>
#define ASSERT(e) assert(e)
#else
#define ASSERT(e)
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "SuiteSparse_config.h"

#define GB_PUBLIC extern

// The largest valid dimension permitted in this implementation is 2^60.
#define GxB_INDEX_MAX ((uint64_t) (1ULL << 60))

#define GB_IMAX(x,y) (((x) > (y)) ? (x) : (y))
#define GB_IMIN(x,y) (((x) < (y)) ? (x) : (y))
#define GB_IMPLIES(p,q) (!(p) || (q))

GB_PUBLIC
bool GB_size_t_multiply     // true if ok, false if overflow
(
    size_t *c,              // c = a*b, or zero if overflow occurs
    const size_t a,
    const size_t b
) ;


// # of bits in an unsigned long long, normally 64
#define GB_64 (8 * sizeof (unsigned long long))

// floor and ceiling of the log2 of an integer.
#ifdef __GNUC__
// see https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html
#define GB_CLZLL(k)   __builtin_clzll ((unsigned long long) (k))
#define GB_CEIL_LOG2(k)  ((uint64_t) ((k) < 2) ? 0 : (GB_64 - GB_CLZLL ((k)-1)))
#define GB_FLOOR_LOG2(k) ((uint64_t) ((k) < 2) ? 0 : (GB_64 - GB_CLZLL (k) - 1))
#else
#define GB_CLZLL(k)   not defined, using log2 instead
#define GB_CEIL_LOG2(k)  ((uint64_t) (ceil  (log2 ((double) k))))
#define GB_FLOOR_LOG2(k) ((uint64_t) (floor (log2 ((double) k))))
#endif

// GB_IS_POWER_OF_TWO(k) is true if the unsigned integer k is an exact power of
// two, or if k is zero.  This expression should not be used if k is negative.
#define GB_IS_POWER_OF_TWO(k) (((k) & ((k) - 1)) == 0)


GB_PUBLIC int      GB_Global_memtable_n (void) ;
GB_PUBLIC void     GB_Global_memtable_dump (void) ;
GB_PUBLIC void     GB_Global_memtable_clear (void) ;
GB_PUBLIC void     GB_Global_memtable_add (void *p, size_t size) ;
GB_PUBLIC size_t   GB_Global_memtable_size (void *p) ;
GB_PUBLIC void     GB_Global_memtable_remove (void *p) ;
GB_PUBLIC bool     GB_Global_memtable_find (void *p) ;

GB_PUBLIC void     GB_Global_free_pool_init (bool clear) ;
GB_PUBLIC void    *GB_Global_free_pool_get (int k) ;
GB_PUBLIC bool     GB_Global_free_pool_put (void *p, int k) ;
GB_PUBLIC void     GB_Global_free_pool_dump (int pr) ;
GB_PUBLIC int64_t  GB_Global_free_pool_limit_get (int k) ;
GB_PUBLIC void     GB_Global_free_pool_limit_set (int k, int64_t nblocks) ;
GB_PUBLIC int64_t  GB_Global_free_pool_nblocks_total (void) ;


//------------------------------------------------------------------------------
// memory management
//------------------------------------------------------------------------------

GB_PUBLIC   // accessed by the MATLAB tests in GraphBLAS/Test only
void *GB_calloc_memory      // pointer to allocated block of memory
(
    size_t nitems,          // number of items to allocate
    size_t size_of_item,    // sizeof each item
    // output
    size_t *size_allocated  // # of bytes actually allocated
) ;

GB_PUBLIC   // accessed by the MATLAB tests in GraphBLAS/Test only
void *GB_malloc_memory      // pointer to allocated block of memory
(
    size_t nitems,          // number of items to allocate
    size_t size_of_item,    // sizeof each item
    // output
    size_t *size_allocated  // # of bytes actually allocated
) ;

GB_PUBLIC   // accessed by the MATLAB tests in GraphBLAS/Test only
void GB_free_memory         // free memory, bypassing the free_pool
(
    // input/output
    void **p,               // pointer to allocated block of memory to free
    // input
    size_t size_allocated   // # of bytes actually allocated
) ;

GB_PUBLIC   // accessed by the MATLAB tests in GraphBLAS/Test only
void GB_dealloc_memory      // free memory, return to free_pool or free it
(
    // input/output
    void **p,               // pointer to allocated block of memory to free
    // input
    size_t size_allocated   // # of bytes actually allocated
) ;

GB_PUBLIC   // accessed by the MATLAB tests in GraphBLAS/Test only
void GB_free_pool_finalize (void) ;

#endif
