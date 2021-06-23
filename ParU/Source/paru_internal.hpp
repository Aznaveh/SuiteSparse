/** =========================================================================  /
 * =======================  paru_internal.hpp ===============================  /
 * ========================================================================== */
/*!
 * internal libraries that are not visible to the user
 * @author Aznaveh
 *  */
//#include "ParU.hpp"
#include "Parallel_LU.hpp"

// more here ...

// force debugging off
#ifndef NDEBUG
#define NDEBUG
#endif

// heap related
void paru_make_heap(Int f, Int start_fac, std::vector<Int> &pivotal_elements,
                    heaps_info &hi, std::vector<Int> &colHash,
                    paru_matrix *paruMatInfo);

void paru_make_heap_empty_el(Int f, std::vector<Int> &pivotal_elements,
                             heaps_info &hi, paru_matrix *paruMatInfo);

