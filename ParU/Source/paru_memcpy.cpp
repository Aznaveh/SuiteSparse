////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_memcpy ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!  @brief  wrapper around memcpy
 *
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
#define MEM_CHUNK (1024 * 1024)  // hard coded chunk-size

void paru_memcpy(void *destination, const void *source, size_t num)
{
    if (num < MEM_CHUNK)
    {  // single task memcpy
        memcpy(destination, source, num);
    }
    else
    {  // multiple task memcpy
        size_t nchunks = 1 + (num / MEM_CHUNK);

        int64_t k;
        #pragma omp taskloop private(k) num_tasks(nchunks)
        for (k = 0; k < (int64_t)nchunks; k++)
        {
            size_t start = k * MEM_CHUNK;
            if (start < num)
            {
                size_t chunk = MIN(num - start, MEM_CHUNK);
                // void* arithmetic is illegal it is why I am using this
                unsigned char *pdest = (unsigned char *)destination + start;
                const unsigned char *psrc = (unsigned char *)source + start;
                memcpy(pdest, psrc, chunk);
            }
        }
    }
}
