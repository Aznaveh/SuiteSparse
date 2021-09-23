/** =========================================================================  /
 * =======================  paru_memset  ====================================  /
 * ========================================================================== */
/*!  @brief  wrapper around memcpy
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"

void paru_memcpy(void *destination, const void *source, size_t num)
{
    // TODO: consider make it parallel
    memcpy(destination, source, num);
}
