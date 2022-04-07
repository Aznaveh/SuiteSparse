////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_memcpy ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*!  @brief  wrapper around memcpy
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

void paru_memcpy(void *destination, const void *source, 
        size_t num, ParU_Control *Control)
{
    size_t mem_chunk = Control->mem_chunk;
    if (num < mem_chunk)
    {  // single task memcpy
        memcpy(destination, source, num);
    }
    else
    {  // multiple task memcpy
        size_t nchunks = 1 + (num / mem_chunk);

        int64_t k;
        #pragma omp taskloop 
        for (k = 0; k < (int64_t)nchunks; k++)
        {
            size_t start = k * mem_chunk;
            if (start < num)
            {
                size_t chunk = MIN(num - start, mem_chunk);
                // void* arithmetic is illegal it is why I am using this
                unsigned char *pdest = (unsigned char *)destination + start;
                const unsigned char *psrc = (unsigned char *)source + start;
                memcpy(pdest, psrc, chunk);
            }
        }
    }
}
