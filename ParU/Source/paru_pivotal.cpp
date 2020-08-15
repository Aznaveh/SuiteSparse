/** =========================================================================  /
 * =======================  paru_pivotal ====================================  /
 * ========================================================================== */

/*! @brief 
 *  adding the list of pivotal elements from the heap, computing the list of
 *  rows and assembling pivotal columns
 *
 * @param pivotal_elements list
 *          f and paruMatInfo
 *
 *  @author Aznaveh
 */
#include "Parallel_LU.hpp"

void paru_pivotal (paru_matrix *paruMatInfo, std::vector<Int> &pivotal_elements,
        Int f)
{
    DEBUGLEVEL(1);
#ifndef NDEBUG
    Int p = 2;
#endif 
 
    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    std::vector<Int>** heapList = paruMatInfo->heapList;
    Int eli = snM [f]; 

    Int *Super = LUsym->Super;
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];


    std::vector<Int>* elHeap = heapList[eli] ;
    paru_Element **elementList = paruMatInfo->elementList;
    auto greater = [&elementList](Int a, Int b)
    { return lnc_el(elementList,a) > lnc_el(elementList,b); };


    PRLEVEL (p, ("%% col2= %ld\n", col2));
    while ( lnc_el(elementList, elHeap->front()) < col2 && elHeap->size() > 0)
    {
        PRLEVEL (p, ("%% elHeap->front = %ld ", elHeap->front()));
        PRLEVEL (p, (" lnc_el = %ld ", 
                    lnc_el(elementList, elHeap->front())));
        PRLEVEL (p, ("%% elHeap->size= %ld \n", elHeap->size()));
        pivotal_elements.push_back(elHeap->front());
        std::pop_heap(elHeap->begin(), elHeap->end(), greater);
        elHeap->pop_back();
    }

#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%% pivotal elements of eli(%ld): ", eli));
    for(Int i=0 ; i < pivotal_elements.size(); i++)
        PRLEVEL (p, ("%ld ", pivotal_elements[i]));
    PRLEVEL (p, ("\n"));
#endif 
    return;
}
