
#include "ParU.hpp"

// more here ...

// force debugging off
#ifndef NDEBUG
#define NDEBUG
#endif

// heap related
void paru_make_heap(Int f, Int start_fac, std::vector<Int> &pivotal_elements,
                    heaps_info &hi, std::vector<Int> &colHash,
                    paru_matrix *paruMatInfo);
