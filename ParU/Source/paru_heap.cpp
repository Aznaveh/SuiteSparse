/** =========================================================================  /
 * =======================  paru_heap.cpp ===================================  /
 * ========================================================================== */
/*! @brief  Wrappers for handling heap
 *
 * @author Aznaveh
 * 
 */
#include "Parallel_LU.hpp"
#define HEAP_ToL 8  //tolerance on when to decide heapify or just add one by one

void perc_down (Int i, Int *lacList, std::vector<Int> &heap)
    // ith position go down into the tree to find its proper position
{
    DEBUGLEVEL(0);
    PRLEVEL (1, ("%%haepsize = %ld\n", heap.size() ));
    while (2*i+1 < heap.size() )
    {
        Int child1 = 2*i+1; // left
        Int child2 = 2*i+2; // right

        PRLEVEL (1, (" lchild=%ld, rchild=%ld\n", child1, child2));
        if (child2 <= heap.size()-1 )
        { // 2 chidren 
            PRLEVEL (1, ("%%Two children\n"));
            if (lacList[heap [i]] < lacList [heap [child1]]  && 
                    lacList[heap [i]] < lacList [heap [child2]] )
                //correct
                break;

            if ( lacList [heap [child1]] < lacList [heap [child2]] )
            {
                PRLEVEL (1, ("%%swapping left child\n"));
                Int tmp = heap [i];
                heap [i] = heap [child1];
                heap [child1] = tmp;
                i = child1;
            }
            else
            {
                PRLEVEL (1, ("%%swapping right child\n"));
                Int tmp = heap [i];
                heap [i] = heap [child2];
                heap [child2] = tmp;
                i = child2;
            }
        }
        else
        { // just left child
            PRLEVEL (1, ("%%Only one child\n"));
            if (lacList[heap [i]] < lacList [heap [child1]] )
                //correct
                break;
            else
            {
                Int tmp = heap [i];
                heap [i] = heap [child1];
                heap [child1] = tmp;
                i = child1;
            }
        }
    }
}

void remove_heap (Int i, Int *lacList, std::vector<Int> &heap)
    // remove ith element from the middle of the heap
    // lacList contatine keys of the min heap
    // i index /heap[i] element/ lacList[heap[i]] key/
{
    DEBUGLEVEL(0);
    PRLEVEL (1, ("%%Removing %ld\n", i));
    Int e = heap [i] = heap.back();
    PRLEVEL (1, (" %ld-%ld(%ld)\n", i, e, lacList [e]));
    Int par = (i-1)/2;
    Int parEl = heap [par];  //parent element
    PRLEVEL (1, ("%ld-", par));
    PRLEVEL (1, ("%ld(%ld)\n", parEl, lacList [parEl]));
    if ( i == 0 || lacList [parEl] < lacList [e] )
    {
        // move down
        PRLEVEL (1, ("%%move down\n"));
        perc_down (i, lacList, heap);
    }
    else
    {
        PRLEVEL (1, ("%%move up\n"));
        while ( lacList [parEl] > lacList [e] )
        {
            //swap
            heap [i] = parEl;
            heap [par] = e;
            i = par;
            par = (i-1)/2;
            parEl = heap [par];  //parent element
            e = heap [i]; 
        }
    }
    heap.pop_back();
#ifndef NDEBUG  
    Int p = 0;
    //chekcing the heap
    for(Int i = heap.size()-1 ; i > 0; i--)
    {
        Int elid = heap [i];
        Int pelid = heap [(i-1)/2]; //parent id
        if ( lacList[pelid] > lacList[elid] )
        {
            PRLEVEL (p, ("%ld-%ld(%ld) <", (i-1)/2, pelid, lacList[pelid] ));
            PRLEVEL (p, ("%ld-%ld(%ld) \n", i, elid, lacList[elid] ));
        }
        ASSERT ( lacList[pelid] <= lacList[elid]);
    }
#endif
}


void paru_make_heap (Int f, std::vector<Int> &pivotal_elements, 
        paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(1);
#ifndef NDEBUG  
    Int p = 1;
#endif

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *aChild = LUsym->aChild;
    Int *aChildp = LUsym->aChildp;
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
    //TODO updating rowMark should go to pivotal update out of here
    // Put a part of this into pivotal and a part into heap for after prior
    // block assembly
    for (Int i = aChildp[eli]; i <= aChildp[eli+1]-1; i++) 
    { //finding the largest child
        Int chelid = aChild[i];  // element id of the child

        PRLEVEL (p, ("%% chelid = %ld\n", chelid));
        std::vector<Int>* curHeap = heapList[chelid];
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
#ifndef NDEBUG  
    PRLEVEL (p, ("%% element ids:\n %%"));
    for(Int i = 0; i < curHeap->size(); i++)
        PRLEVEL (p, (" %ld", (*curHeap)[i]));
    PRLEVEL (p, ("\n"));
#endif
 
    }
    Int *lacList = paruMatInfo -> lacList;
    auto greater = [&lacList](Int a, Int b){ return lacList[a] > lacList[b]; };

    PRLEVEL (p, ("%% tot_size =  %ld\n", tot_size ));
    PRLEVEL (p+1, ("%% biggest_Child_id = %ld\n", biggest_Child_id));
    PRLEVEL (p, ("%% biggest_Child_size = %ld\n", biggest_Child_size));

    if (biggest_Child_id != -1) //everything already freed
    {

        //shallow copy of the biggest child
        std::vector<Int>* elHeap = heapList[eli] = heapList[biggest_Child_id];
        heapList[biggest_Child_id] = nullptr;

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

    }
    else
    {
        //TODO to see if we can put currentheap null
        std::vector<Int>* elHeap = heapList[eli]  = new std::vector<Int>;
    }


    //fixing the current heap
    std::vector<Int>* curHeap = heapList[eli];
    ASSERT (curHeap != nullptr);

    curHeap->push_back(eli);
    std::push_heap(curHeap->begin(), curHeap->end(), greater);
    PRLEVEL (p, ("%% %ld pushed ",eli));

    for(Int i = 0 ; i < pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        paru_Element *el = elementList[e];
        if (el == NULL) continue;
        PRLEVEL (p, ("%ld  ",e));
        curHeap->push_back(e);
        std::push_heap(curHeap->begin(), curHeap->end(), greater);
    }

    PRLEVEL (p, ("After everything eli %ld has %ld elements\n", 
             eli, curHeap->size() ));

#ifndef NDEBUG
    PRLEVEL (p, ("%%Heap after making it(size = %ld) \n", curHeap->size() ));
    for(Int i = 0; i < curHeap->size(); i++)
    {
        Int elid = (*curHeap)[i];
        PRLEVEL (p, (" %ld(%ld) ", elid, lacList[elid] ));
    }
    PRLEVEL (p, ("\n"));
    for(Int i = curHeap->size()-1 ; i > 0; i--)
    {
        Int elid = (*curHeap)[i];
        Int pelid = (*curHeap)[(i-1)/2]; //parent id
        if( lacList[pelid] > lacList[elid])
            PRLEVEL (p, ("ATT %ld(%ld)\n\n ", elid, lacList[elid]));
        ASSERT ( lacList[pelid] <= lacList[elid]);
    }
#endif
}
