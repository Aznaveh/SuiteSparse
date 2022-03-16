////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_heap.cpp //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  Wrappers for handling heap
 *
 * @author Aznaveh
 *
 */
#include "paru_internal.hpp"

void paru_check_prior_element(Int e, Int f, Int start_fac,
                              std::vector<Int> &colHash,
                              paru_matrix *paruMatInfo)
// check if e can be assembeld into f
{
    Paru_Work *Work = paruMatInfo->Work;
    Int *elRow = Work->elRow;

    ParU_Element **elementList = paruMatInfo->elementList;

    ParU_Element *el = elementList[e];
    if (elRow[e] == 0 && el->rValid > start_fac)
    {  // all the rows are inside he current front; maybe assemble some cols
        paru_assemble_cols(e, f, colHash, paruMatInfo);
        return;
    }

    //    if ( (elCol [e] == 0 && el->cValid > start_fac) ||
    //            el->cValid == paruMatInfo->time_stamp[f])
    if (el->rValid == start_fac || el->cValid == paruMatInfo->time_stamp[f])
    {  // all the cols are inside he current front; maybe assemble some rows
        paru_assemble_rows(e, f, colHash, paruMatInfo);
    }
}

ParU_Ret paru_make_heap(Int f, Int start_fac,
                               std::vector<Int> &pivotal_elements,
                               heaps_info &hi, std::vector<Int> &colHash,
                               paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;

    ParU_Symbolic *Sym = paruMatInfo->Sym;
    Int *aChild = Sym->aChild;
    Int *aChildp = Sym->aChildp;
    Int *snM = Sym->super2atree;
    ParU_Element **elementList = paruMatInfo->elementList;
    // Int m = paruMatInfo-> m;

    std::vector<Int> **heapList = paruMatInfo->heapList;

    Int eli = snM[f];

    // element id of the biggest child
    Int biggest_Child_id = -1;
    Int biggest_Child_size = -1;
    Int tot_size = 0;

    biggest_Child_id = hi.biggest_Child_id;
    biggest_Child_size = hi.biggest_Child_size;
    tot_size = hi.sum_size;

    Int *lacList = paruMatInfo->lacList;
    auto greater = [&lacList](Int a, Int b) { return lacList[a] > lacList[b]; };

    PRLEVEL(PR, ("%% tot_size =  %ld\n", tot_size));
    PRLEVEL(PR, ("%% biggest_Child_id = %ld ", biggest_Child_id));
    PRLEVEL(PR, ("%% biggest_Child_size = %ld\n", biggest_Child_size));
    Int size_of_rest = tot_size - biggest_Child_size + pivotal_elements.size();
    PRLEVEL(PR, ("%% the rest size = %ld\n", size_of_rest));

    if (biggest_Child_id != -1)
    // There are still elements remained in the heaps
    {
        // shallow copy of the biggest child
        std::vector<Int> *curHeap = heapList[eli] = heapList[biggest_Child_id];
        heapList[biggest_Child_id] = nullptr;

        // O(n) heapify of all children or O(klgn) add to the biggest child
        if (log2(biggest_Child_size) >
            (biggest_Child_size / (size_of_rest + 1)) + 1)
        {  // klogn
            PRLEVEL(PR, ("%% klogn algorhtm\n"));
            for (Int i = aChildp[eli]; i <= aChildp[eli + 1] - 1; i++)
            {
                Int chelid = aChild[i];  // element id of the child
                std::vector<Int> *chHeap = heapList[chelid];
                if (chHeap == nullptr) continue;
                // concatening the child and freeing the memory
                for (Int k = 0; k < (Int)chHeap->size(); k++)
                {
                    Int e = (*chHeap)[k];
                    if (elementList[e] != NULL)
                    {
                        paru_check_prior_element(e, f, start_fac, colHash,
                                                 paruMatInfo);
                        if (elementList[e] != NULL)
                        {
                            curHeap->push_back(e);
                            std::push_heap(curHeap->begin(), curHeap->end(),
                                           greater);
                        }
                    }
                }
                delete heapList[chelid];
                heapList[chelid] = nullptr;
            }

            for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
            {
                Int e = pivotal_elements[i];
                ParU_Element *el = elementList[e];
#ifndef NDEBUG
                ASSERT(el != NULL);
#endif
                if (el == NULL) continue;
                PRLEVEL(PR, ("%ld  ", e));
                curHeap->push_back(e);
                std::push_heap(curHeap->begin(), curHeap->end(), greater);
            }
            curHeap->push_back(eli);
            std::push_heap(curHeap->begin(), curHeap->end(), greater);
            PRLEVEL(PR, ("%% %ld pushed ", eli));
        }
        else
        {  // heapify
            PRLEVEL(PR, ("%%heapify with the size %ld\n", tot_size));
            for (Int i = aChildp[eli]; i <= aChildp[eli + 1] - 1; i++)
            {
                Int chelid = aChild[i];  // element id of the child
                std::vector<Int> *chHeap = heapList[chelid];
                if (chHeap == nullptr) continue;
                // concatening the child and freeing the memory

                // curHeap->insert(curHeap->end(),
                //      chHeap->begin(), chHeap->end());
                for (Int k = 0; k < (Int)chHeap->size(); k++)
                {
                    Int e = (*chHeap)[k];
                    if (elementList[e] != NULL)
                    {
                        paru_check_prior_element(e, f, start_fac, colHash,
                                                 paruMatInfo);
                        if (elementList[e] != NULL) curHeap->push_back(e);
                    }
                }
                PRLEVEL(1,
                        ("%%Heap free %p id=%ld\n", heapList[chelid], chelid));
                delete heapList[chelid];
                heapList[chelid] = nullptr;
            }
            // adding pivotal elements and the current element
            curHeap->insert(curHeap->end(), pivotal_elements.begin(),
                            pivotal_elements.end());
            curHeap->push_back(eli);
            // heapifying
            std::make_heap(curHeap->begin(), curHeap->end(), greater);
        }
    }
    else
    {
        PRLEVEL(PR, ("Nothing in the heap. size of pivotal %ld \n",
                    pivotal_elements.size()));
        std::vector<Int> *curHeap; 
        try
        {
            curHeap = heapList[eli] = new std::vector<Int>;
        }
        catch (std::bad_alloc const&)
        {  // out of memory
            return PARU_OUT_OF_MEMORY;
        }
        // deep copy
        //*curHeap = pivotal_elements;
        // swap provides a shallow copy
        std::swap(*curHeap, pivotal_elements);
        curHeap->push_back(eli);
        std::make_heap(curHeap->begin(), curHeap->end(), greater);
    }

#ifndef NDEBUG
    std::vector<Int> *curHeap = heapList[eli];
    PRLEVEL(PR, ("After everything eli %ld has %ld elements\n", eli,
                curHeap->size()));
    PRLEVEL(PR, ("%%Heap after making it(size = %ld) \n", curHeap->size()));
    for (Int i = 0; i < (Int)curHeap->size(); i++)
    {
        Int elid = (*curHeap)[i];
        PRLEVEL(PR, (" %ld(%ld) ", elid, lacList[elid]));
    }
    PRLEVEL(PR, ("\n"));
    for (Int i = curHeap->size() - 1; i > 0; i--)
    {
        Int elid = (*curHeap)[i];
        Int pelid = (*curHeap)[(i - 1) / 2];  // parent id
        if (lacList[pelid] > lacList[elid])
            PRLEVEL(PR, ("ATT %ld(%ld)\n\n ", elid, lacList[elid]));
        ASSERT(lacList[pelid] <= lacList[elid]);
    }
#endif
    return PARU_SUCCESS;
}

ParU_Ret paru_make_heap_empty_el(Int f,
                                        std::vector<Int> &pivotal_elements,
                                        heaps_info &hi,
                                        paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;

    ParU_Symbolic *Sym = paruMatInfo->Sym;
    Int *aChild = Sym->aChild;
    Int *aChildp = Sym->aChildp;
    Int *snM = Sym->super2atree;
    ParU_Element **elementList = paruMatInfo->elementList;
    // Int m = paruMatInfo-> m;

    std::vector<Int> **heapList = paruMatInfo->heapList;

    Int eli = snM[f];

    // element id of the biggest child
    Int biggest_Child_id = -1;
    Int biggest_Child_size = -1;
    Int tot_size = 0;

    biggest_Child_id = hi.biggest_Child_id;
    biggest_Child_size = hi.biggest_Child_size;
    tot_size = hi.sum_size;

    Int *lacList = paruMatInfo->lacList;
    auto greater = [&lacList](Int a, Int b) { return lacList[a] > lacList[b]; };

    PRLEVEL(PR, ("%% tot_size =  %ld\n", tot_size));
    PRLEVEL(PR, ("%% biggest_Child_id = %ld ", biggest_Child_id));
    PRLEVEL(PR, ("%% biggest_Child_size = %ld\n", biggest_Child_size));
    Int size_of_rest = tot_size - biggest_Child_size + pivotal_elements.size();
    PRLEVEL(PR, ("%% the rest size = %ld\n", size_of_rest));

    if (biggest_Child_id != -1)
    // There are still elements remained in the heaps
    {
        // shallow copy of the biggest child
        std::vector<Int> *curHeap = heapList[eli] = heapList[biggest_Child_id];
        heapList[biggest_Child_id] = nullptr;

        // O(n) heapify of all children or O(klgn) add to the biggest child
        if (log2(biggest_Child_size) >
            (biggest_Child_size / (size_of_rest + 1)) + 1)
        {  // klogn
            PRLEVEL(PR, ("%% klogn algorhtm\n"));
            for (Int i = aChildp[eli]; i <= aChildp[eli + 1] - 1; i++)
            {
                Int chelid = aChild[i];  // element id of the child
                std::vector<Int> *chHeap = heapList[chelid];
                if (chHeap == nullptr) continue;
                // concatening the child and freeing the memory
                for (Int k = 0; k < (Int)chHeap->size(); k++)
                {
                    Int e = (*chHeap)[k];
                    if (elementList[e] != NULL)
                    {
                        // paru_check_prior_element(e, f, start_fac, colHash,
                        //                         paruMatInfo);
                        if (elementList[e] != NULL)
                        {
                            curHeap->push_back(e);
                            std::push_heap(curHeap->begin(), curHeap->end(),
                                           greater);
                        }
                    }
                }
                delete heapList[chelid];
                heapList[chelid] = nullptr;
            }

            for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
            {
                Int e = pivotal_elements[i];
                ParU_Element *el = elementList[e];
                if (el == NULL) continue;
                PRLEVEL(PR, ("%ld  ", e));
                curHeap->push_back(e);
                std::push_heap(curHeap->begin(), curHeap->end(), greater);
            }
            std::push_heap(curHeap->begin(), curHeap->end(), greater);
            PRLEVEL(PR, ("%% %ld pushed ", eli));
        }
        else
        {  // heapify
            PRLEVEL(PR, ("%%heapify with the size %ld\n", tot_size));
            for (Int i = aChildp[eli]; i <= aChildp[eli + 1] - 1; i++)
            {
                Int chelid = aChild[i];  // element id of the child
                std::vector<Int> *chHeap = heapList[chelid];
                if (chHeap == nullptr) continue;
                // concatening the child and freeing the memory

                // curHeap->insert(curHeap->end(),
                //      chHeap->begin(), chHeap->end());
                for (Int k = 0; k < (Int)chHeap->size(); k++)
                {
                    Int e = (*chHeap)[k];
                    if (elementList[e] != NULL)
                    {
                        if (elementList[e] != NULL) curHeap->push_back(e);
                    }
                }
                PRLEVEL(1,
                        ("%%Heap free %p id=%ld\n", heapList[chelid], chelid));
                delete heapList[chelid];
                heapList[chelid] = nullptr;
            }
            // adding pivotal elements and the current element
            curHeap->insert(curHeap->end(), pivotal_elements.begin(),
                            pivotal_elements.end());
            // heapifying
            std::make_heap(curHeap->begin(), curHeap->end(), greater);
        }
    }
    else
    {
        PRLEVEL(PR, ("Nothing in the heap. size of pivotal %ld \n",
                    pivotal_elements.size()));
        std::vector<Int> *curHeap; 
        try
        {
            curHeap = heapList[eli] = new std::vector<Int>;
        }
        catch (std::bad_alloc const&)
        {  // out of memory
            return PARU_OUT_OF_MEMORY;
        }
        // deep copy
        //*curHeap = pivotal_elements;
        // swap provides a shallow copy
        std::swap(*curHeap, pivotal_elements);
        std::make_heap(curHeap->begin(), curHeap->end(), greater);
    }

#ifndef NDEBUG
    std::vector<Int> *curHeap = heapList[eli];
    PRLEVEL(PR, ("After everything eli %ld has %ld elements\n", eli,
                curHeap->size()));
    PRLEVEL(PR, ("%%Heap after making it(size = %ld) \n", curHeap->size()));
    for (Int i = 0; i < (Int)curHeap->size(); i++)
    {
        Int elid = (*curHeap)[i];
        PRLEVEL(PR, (" %ld(%ld) ", elid, lacList[elid]));
    }
    PRLEVEL(PR, ("\n"));
    for (Int i = curHeap->size() - 1; i > 0; i--)
    {
        Int elid = (*curHeap)[i];
        Int pelid = (*curHeap)[(i - 1) / 2];  // parent id
        if (lacList[pelid] > lacList[elid])
            PRLEVEL(PR, ("ATT %ld(%ld)\n\n ", elid, lacList[elid]));
        ASSERT(lacList[pelid] <= lacList[elid]);
    }
#endif
    return PARU_SUCCESS;
}
