////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_memset  ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!  @brief  wrapper around memset
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

void paru_memset(void *ptr, Int value, size_t num, ParU_Control *Control)
{
    size_t mem_chunk = Control->mem_chunk;
    if (num < mem_chunk)
    {  // single task memse
        memset(ptr, value, num);
    }
    else
    {  // multiple task memset
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
                unsigned char *ptr_chunk = (unsigned char *)ptr + start;
                memset(ptr_chunk, value, chunk);
            }
        }
    }
}

