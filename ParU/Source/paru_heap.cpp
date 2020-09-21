/** =========================================================================  /
 * =======================  paru_heap.cpp ===================================  /
 * ========================================================================== */
/*! @brief  Wrappers for handling heap
 *
 * @author Aznaveh
 * 
 */
#include "Parallel_LU.hpp"
#define HEAP_ToL 10  //tolerance on how to 

void paru_make_heap(paru_matrix *paruMatInfo, Int f )
{
    DEBUGLEVEL(0);
#ifndef NDEBUG  
    Int p = 1;
#endif

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int * aChild = LUsym->aChild;
    Int * aChildp = LUsym->aChildp;
    Int *snM = LUsym->super2atree;
    paru_Element **elementList = paruMatInfo->elementList;
    Int m = paruMatInfo-> m;

    std::vector<Int>** heapList = paruMatInfo->heapList;

    Int eli = snM [f]; 

    // element id of the biggest child
    Int biggest_Child_id = -1;
    Int biggest_Child_size = -1;
    Int tot_size = 0; 
    work_struct *Work =  paruMatInfo->Work;
    Int *rowMarkp = Work->rowMark;
    Int rowMark = 0;
    for (Int i = aChildp[eli]; i <= aChildp[eli+1]-1; i++) 
    {
        Int chelid = aChild[i];  // element id of the child
        // max(rowMark , child->rowMark)
        Int f_rmark = rowMarkp[chelid];
        rowMark = rowMark >  f_rmark ?  rowMark : f_rmark;

        PRLEVEL (p, ("%% chelid = %ld\n", chelid));
        std::vector<Int>* curHeap = heapList[chelid];
        //ASSERT(curHeap != nullptr);
        PRLEVEL (p, ("%% curHeap= %p\n", curHeap));
        if (curHeap == nullptr) continue;
        Int cur_size = curHeap->size();

        PRLEVEL (p, ("%% cur_size =  %ld\n",cur_size));
        tot_size += cur_size; 
        if (cur_size > biggest_Child_size)
        {
            biggest_Child_id = chelid;
            biggest_Child_size = cur_size;
        }
    }
    rowMarkp[eli] = rowMark;
    auto greater = [&elementList](Int a, Int b)
    { return lnc_el(elementList,a) > lnc_el(elementList,b); };

    PRLEVEL (p, ("%% tot_size =  %ld\n", tot_size ));
    PRLEVEL (p+1, ("%% biggest_Child_id = %ld\n", biggest_Child_id));
    PRLEVEL (p, ("%% biggest_Child_size = %ld\n", biggest_Child_size));

    //shallow copy of the biggest child
    std::vector<Int>* elHeap = heapList[eli] = heapList[biggest_Child_id];
    heapList[biggest_Child_id] = nullptr;

    std::make_heap(elHeap->begin(), elHeap->end(), greater ); 
    //O(n) heapify of all children or O(klgn) add to the biggest child
    if ( biggest_Child_size > HEAP_ToL*(tot_size - biggest_Child_size) )
    { //klogn
        PRLEVEL (p, ("%% klogn algorhtm\n"));
        for (Int i = aChildp[eli]; i <= aChildp[eli+1]-1; i++) 
        {   
            Int chelid = aChild[i];  // element id of the child
            std::vector<Int>*  chHeap = heapList[chelid];
            if (chHeap == nullptr) continue;
            //concatening the child and freeing the memory
            for (Int k = 0; k < chHeap->size(); k++)
            {
                Int elid = (*chHeap)[k];
                if (elementList[elid] != NULL)
                {
                    elHeap->push_back(elid);
                    std::push_heap(elHeap->begin(), elHeap->end(), greater);
                }
            }
            delete heapList[chelid];
            heapList[chelid] = nullptr;
        }
    }
    else
    {  //heapify
        PRLEVEL (p, ("%%heapify with the size %ld\n", tot_size));
        for (Int i = aChildp[eli]; i <= aChildp[eli+1]-1; i++) 
        {
            Int chelid = aChild[i];  // element id of the child
            std::vector<Int>* chHeap = heapList[chelid];
            if (chHeap == nullptr) continue;
            //concatening the child and freeing the memory
            elHeap->insert(elHeap->end(), chHeap->begin(), chHeap->end()); 
            PRLEVEL (1, ("%%Heap free %p id=%ld\n", heapList[chelid], chelid));
            delete heapList[chelid];
            heapList[chelid] = nullptr;
            //heapifying
            std::make_heap(elHeap->begin(), elHeap->end(), greater ); 

        }
    }

#ifndef NDEBUG  
    PRLEVEL (p, ("%% element ids:\n %%"));
    for(Int i = 0; i < elHeap->size(); i++)
        PRLEVEL (p, (" %ld", (*elHeap)[i]));
    PRLEVEL (p, ("\n"));
    PRLEVEL (p, ("%% first cols:\n %%"));
    for(Int i = 0; i < elHeap->size(); i++)
    {
        Int elid = (*elHeap)[i];
        if (elid != NULL)
            PRLEVEL (p, (" %ld", lnc_el(elementList, elid) ));
    }
    PRLEVEL (p, ("\n"));
    //chekcing the heap
    for(Int i = elHeap->size()-1 ; i > 0; i--)
    {
        Int elid = (*elHeap)[i];
        Int pelid = (*elHeap)[(i-1)/2]; //parent id
        ASSERT ( lnc_el(elementList,pelid) <= lnc_el(elementList,elid));
    }

#endif
}
