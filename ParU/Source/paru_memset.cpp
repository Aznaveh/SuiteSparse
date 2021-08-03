/** =========================================================================  /
 * =======================  paru_memset  ====================================  /
 * ========================================================================== */
/*!  @brief  wrapper around memset 
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

void paru_memset (void * ptr, Int value, size_t num)
{
    memset(ptr, value, num);
}
