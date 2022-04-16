////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_cov.hpp ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

#ifndef PARU_COVERAGE_H
#define PARU_COVERAGE_H

#include "paru_internal.hpp"
//!
//  internal libraries that are not visible to the user
//  @author Aznaveh
//

/* sample usage:

    To test:

        int info = paru_stuff ( ...) ;

    do this instead:

        paru_set_malloc_tracking (true) ;
        for (nmalloc = 0 ; ; nmalloc++)
        {
            paru_set_nmalloc (nmalloc) ;
            // do stuff
            int info = paru_stuff (...) ;
            if (info != PARU_OUT_OF_MEMORY) break ;
        }
        paru_set_malloc_tracking (false) ;

    or:

        int info ;
        BRUTAL_ALLOC_TEST (info, paru_sym (...)) ;
        if (info == PARU_SUCCESS)
        {
            BRUTAL_ALLOC_TEST (info, paru_num (...)) ;
        }
        if (info == PARU_SUCCESS)
        {
            BRUTAL_ALLOC_TEST (info, paru_solve (...)) ;
        }
        TEST_ASSERT (info == PARU_SUCCESS)
        free sym, num, soln

    // put in test coverage *.h file:
    #ifdef PARU_ALLOC_TESTING
    #define BRUTAL_ALLOC_TEST(info,method)              \
    {                                                   \
        paru_set_malloc_tracking (true) ;               \
        for (Int nmalloc = 0 ; ; Int nmalloc++)         \
        {                                               \
            paru_set_nmalloc (nmalloc)                  \
            // do stuff                                 \
            info = method ;                             \
            if (info != PARU_OUT_OF_MEMORY) break ;     \
            if (nmalloc > 1000000) { test failure }
        }                                               \
        paru_set_malloc_tracking (false) ;              \
    }
    #else
    #define BRUTAL_ALLOC_TEST(info,method)              \
    {                                                   \
        info = method ;                                 \
    }
    #endif
*/

#endif
