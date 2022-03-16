////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_prior_assemble ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief numerical assemble of prior fronts
 *
 * @author Aznaveh
 */

#include "paru_internal.hpp"

ParU_Ret paru_prior_assemble(Int f, Int start_fac,
                                    std::vector<Int> &pivotal_elements,
                                    std::vector<Int> &colHash, heaps_info &hi,
                                    paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;

    Paru_Work *Work = paruMatInfo->Work;
    Int *elCol = Work->elCol;

    ParU_Element **elementList = paruMatInfo->elementList;
    ParU_Symbolic *Sym = paruMatInfo->Sym;
    Int *snM = Sym->super2atree;

    Int pMark = start_fac;

#ifndef NDEBUG
    Int *elRow = Work->elRow;
    Int el_ind = snM[f];
    PRLEVEL(PR, ("%%Inside prior\n"));
    PRLEVEL(PR, ("%% pivotal size is %ld ", pivotal_elements.size()));

#endif
    Int ii = 0;

    for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        ParU_Element *el = elementList[e];
        PRLEVEL(PR, ("%% element= %ld  \n", e));
        if (el == NULL)
        {
            PRLEVEL(PR, ("%% element= %ld is NULL ii=%ld \n", e, ii));
            continue;
        }
#ifndef NDEBUG
        PRLEVEL(PR, ("%%elRow[%ld]=%ld \n", e, elRow[e]));
        // if (elRow[e] != 0) PRLEVEL(-1, ("%%elRow[%ld]=%ld \n", e, elRow[e]));
        // ASSERT (elRow[e] == 0);
#endif

        if (el->nzr_pc == 0)  // if all the rows are available in current front
        {
            if (el->rValid == pMark || elCol[e] == 0)
            // it can be fully assembled
            // both a pivotal column and pivotal row
            {
                #ifndef NDEBUG
                PRLEVEL(PR, ("%%assembling %ld in %ld\n", e, el_ind));
                PRLEVEL(PR, ("%% size %ld x %ld\n", el->nrows, el->ncols));
                #endif
                paru_assemble_all(e, f, colHash, paruMatInfo);
                #ifndef NDEBUG
                PRLEVEL(PR, ("%%assembling %ld in %ld done\n", e, el_ind));
                #endif
                continue;
            }

            #ifndef NDEBUG
            PRLEVEL(PR, ("%%assembling %ld in %ld\n", e, el_ind));
            #endif
            paru_assemble_cols(e, f, colHash, paruMatInfo);
            #ifndef NDEBUG
            PRLEVEL(PR, ("%%partial col assembly%ld in %ld done\n", e, el_ind));
            #endif
            if (elementList[e] == NULL) continue;
        }
        else
        {
            if (el->rValid == pMark || elCol[e] == 0)
            // This element contributes to both pivotal rows and pivotal columns
            //  However it has zero rows in current pivotal columns therefore
            //  not all rows are there
            // it can be assembled partially
            //       ________________________________
            //       |      |                         |
            //       |      |                         |
            //       ___xxxxxxxxxxx____________________
            //       |  xxxxxxxxxxx                   |
            //       |  oxxo|oxoxox                   | <- assemble rows
            //       |  ooxx|oxoxox                   |
            //       |  oooo|oxoxox                   |
            //       ---------------------------------
            //          ooooooxxxxx  --> outsidie the front
            //          ooooooxxxxx
            //
            {
                paru_assemble_el_with0rows(e, f, colHash, paruMatInfo);
                if (elementList[e] == NULL) continue;
                #ifndef NDEBUG
                PRLEVEL(PR, ("%%assembling %ld in %ld done\n", e, el_ind));
                #endif
            }
            // keeping current element
        }

        pivotal_elements[ii++] = pivotal_elements[i];
    }

    if (ii < (Int)pivotal_elements.size())
    {
        PRLEVEL(PR, ("%% Prior: size was %ld ", pivotal_elements.size()));
        PRLEVEL(PR, (" and now is %ld\n ", ii));
        pivotal_elements.resize(ii);
    }

    /************ Making the heap from list of the immediate children
     * ******/
    PRLEVEL(1, ("%% Next: work on the heap \n"));
    ParU_Ret res_make_heap;
    res_make_heap = paru_make_heap(f, start_fac, pivotal_elements, hi, colHash,
                                   paruMatInfo);
    if (res_make_heap != PARU_SUCCESS) return res_make_heap;
    PRLEVEL(1, ("%% Done: work on the heap \n"));

    Int eli = snM[f];
    std::vector<Int> **heapList = paruMatInfo->heapList;
    std::vector<Int> *curHeap = heapList[eli];

    if (curHeap->empty()) return PARU_SUCCESS;

#ifndef NDEBUG
    PR = 1;
#endif

#ifndef NDEBUG
    Int *lacList = paruMatInfo->lacList;
    PRLEVEL(PR, ("%% current heap:\n %%"));
    for (Int k = 0; k < (Int)curHeap->size(); k++)
    {
        Int ee = (*curHeap)[k];
        ParU_Element *ell = elementList[ee];
        PRLEVEL(PR, ("%ld-%ld", k, ee));
        if (ell != NULL)
        {
            PRLEVEL(PR, ("(%ld) ", lacList[ee]));
        }
        else
        {
            PRLEVEL(PR, ("(*%ld) ", lacList[ee]));
        }
    }
    PRLEVEL(PR, ("\n"));
#endif

#ifndef NDEBUG
    // chekcing the heap
    for (Int i = curHeap->size() - 1; i > 0; i--)
    {
        Int elid = (*curHeap)[i];
        Int pelid = (*curHeap)[(i - 1) / 2];  // parent id
        if (lacList[pelid] > lacList[elid])
        {
            PRLEVEL(PR, ("%ld-%ld(%ld) <", (i - 1) / 2, pelid, lacList[pelid]));
            PRLEVEL(PR, ("%ld-%ld(%ld) \n", i, elid, lacList[elid]));
        }
        ASSERT(lacList[pelid] <= lacList[elid]);
    }

#endif
    return PARU_SUCCESS;
}
