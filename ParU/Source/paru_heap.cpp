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
    DEBUGLEVEL(1);
#ifndef NDEBUG  
    Int p = 1;
#endif

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int * aChild = LUsym->aChild;
    Int * aChildp = LUsym->aChildp;
    Int *snM = LUsym->super2atree;

    std::vector<Int>** heapList = paruMatInfo->heapList;

    Int eli = snM [f]; 

    // element id of two biggest child
    Int biggest_Child_id = -1;
    Int biggest_Child_size = -1;
    Int next_Biggest_id = -1;
    Int next_Biggest_size = -1;
    Int tot_size = 0; // 
    for (Int i = aChildp[eli]; i <= aChildp[eli+1]-1; i++) 
    {
        Int chelid = aChild[i];  // element id of the child
        PRLEVEL (p+1, ("%% chelid = %ld\n", chelid));
        std::vector<Int>* curHeap = heapList[chelid];
        //ASSERT(curHeap != nullptr);
        if (curHeap == nullptr) continue;
        Int cur_size = curHeap->size();
        PRLEVEL (p+1, ("%% cur_size =  %ld\n",cur_size));
        tot_size += cur_size; 
        if (cur_size > biggest_Child_size)
        {
            biggest_Child_id = chelid;
            biggest_Child_size = cur_size;
        }
        else if (cur_size > next_Biggest_size)
        {
            next_Biggest_id = chelid;
            next_Biggest_size = cur_size;
        }
    }

    PRLEVEL (p, ("%% tot_size =  %ld\n", tot_size ));
    PRLEVEL (p, ("%% biggest_Child_id = %ld\n", biggest_Child_id));
    PRLEVEL (p, ("%% biggest_Child_size = %ld\n", biggest_Child_size));
    PRLEVEL (p, ("%% next_Biggest_id = %ld\n", next_Biggest_id));
    PRLEVEL (p, ("%% next_Biggest_size = %ld\n", next_Biggest_size));

    //shallow copy of the biggest child
    std::vector<Int>* elHeap = heapList[eli] = heapList[biggest_Child_id];
    heapList[biggest_Child_id] = nullptr;

    //O(n) heapify of all children or O(klgn) add to the biggest child

    if ( biggest_Child_size > HEAP_ToL*(tot_size - biggest_Child_size) )
    { //klogn
        for (Int i = aChildp[eli]; i <= aChildp[eli+1]-1; i++) 
        {   
            Int chelid = aChild[i];  // element id of the child
            std::vector<Int>*  chHeap = heapList[chelid];
            //ASSERT(curHeap != nullptr);
            if (chHeap == nullptr) continue;
            //concatening the child and freeing the memory
            for (Int k = 0; k < chHeap->size(); k++)
            {
                elHeap->push_back((*chHeap)[k]);
                //std::push_heap(elHeap->begin(), elHeap->end(), [](){return true;});
            }
            delete heapList[chelid];
            heapList[chelid] = nullptr;
 
        }
    }
    else
    {  //heapify

        for (Int i = aChildp[eli]; i <= aChildp[eli+1]-1; i++) 
        {
            Int chelid = aChild[i];  // element id of the child
            std::vector<Int>* chHeap = heapList[chelid];
            //ASSERT(curHeap != nullptr);
            if (chHeap == nullptr) continue;
            //concatening the child and freeing the memory
            elHeap->insert(elHeap->end(), chHeap->begin(), chHeap->end()); 
            delete heapList[chelid];
            heapList[chelid] = nullptr;

        }


    }



}
