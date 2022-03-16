////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_memset  ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!  @brief  wrapper around memset
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define MEM_CHUNK (1024 * 1024 * 128)  // hard coded chunk-size //XXX

void paru_memset(void *ptr, Int value, size_t num)
{
    if (num < MEM_CHUNK)
    {  // single task memse
        memset(ptr, value, num);
    }
    else
    {  // multiple task memset
        size_t nchunks = 1 + (num / MEM_CHUNK);

        int64_t k;
        #pragma omp taskloop
        for (k = 0; k < (int64_t)nchunks; k++)
        {
            size_t start = k * MEM_CHUNK;
            if (start < num)
            {
                size_t chunk = MIN(num - start, MEM_CHUNK);
                // void* arithmetic is illegal it is why I am using this
                unsigned char *ptr_chunk = (unsigned char *)ptr + start;
                memset(ptr_chunk, value, chunk);
            }
        }
    }
}

