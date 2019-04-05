/** =========================================================================  /
 * =======================  paru_assemble   =================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

#ifndef NDEBUG  // using STL for debugging
#include <iostream>
#include <algorithm>
#include <set>
#endif

/*! @brief Computing factorization of current front and doing the numerical
 * assembly that ancestors will assemble. Degree update will be used in this
 * version. Just like ./paru_assemble.cpp
 * 
 *  @author Aznaveh
 *
 *  
 *
 *
 * @param  a list of tuples and the the tuple we want to add
 * @return 0 on sucess 
 */
