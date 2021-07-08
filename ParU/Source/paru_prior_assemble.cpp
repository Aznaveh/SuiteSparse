/** =========================================================================  /
 * =======================  paru_prior_assemble =============================  /
 * ========================================================================== */
/*! @brief numerical assemble of prior fronts
 *
 * @author Aznaveh
 */

#include "paru_internal.hpp"

ParU_ResultCode paru_prior_assemble(Int f, Int start_fac,
                         std::vector<Int> &pivotal_elements,
                         std::vector<Int> &colHash, heaps_info &hi,
                         paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);

    work_struct *Work = paruMatInfo->Work;
    Int *elCol = Work->elCol;

    paru_Element **elementList = paruMatInfo->elementList;
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;

    Int pMark = start_fac;

#ifndef NDEBUG
    Int p = 1;
    Int *elRow = Work->elRow;
    Int el_ind = snM[f];
    PRLEVEL(p, ("%%Inside prior\n"));
    PRLEVEL(p, ("%% pivotal size is %ld ", pivotal_elements.size()));

#endif
    Int ii = 0;

    for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        paru_Element *el = elementList[e];
        PRLEVEL(p, ("%% element= %ld  \n", e));
        if (el == NULL)
        {
            PRLEVEL(p, ("%% element= %ld is NULL ii=%ld \n", e, ii));
            continue;
        }
#ifndef NDEBUG
        PRLEVEL(p, ("%%elRow[%ld]=%ld \n", e, elRow[e]));
        // if (elRow[e] != 0) PRLEVEL(-1, ("%%elRow[%ld]=%ld \n", e, elRow[e]));
        // ASSERT (elRow[e] == 0);
#endif

        if (el->nzr_pc == 0)  // if all the rows are available in current front
        {
            if (el->rValid == pMark || elCol[e] == 0)
            // it can be fully assembled
            // both a pivotal column and pivotal row
            {
                PRLEVEL(p, ("%%assembling %ld in %ld\n", e, el_ind));
                PRLEVEL(p, ("%% size %ld x %ld\n", el->nrows, el->ncols));
                paru_assemble_all(e, f, colHash, paruMatInfo);
                PRLEVEL(p, ("%%assembling %ld in %ld done\n", e, el_ind));
                continue;
            }

            PRLEVEL(p, ("%%assembling %ld in %ld\n", e, el_ind));
            paru_assemble_cols(e, f, colHash, paruMatInfo);
            PRLEVEL(p, ("%%partial col assembly%ld in %ld done\n", e, el_ind));
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
                PRLEVEL(p, ("%%assembling %ld in %ld done\n", e, el_ind));
            }
            // keeping current element
        }

        pivotal_elements[ii++] = pivotal_elements[i];
    }

    if (ii < (Int)pivotal_elements.size())
    {
        PRLEVEL(p, ("%% Prior: size was %ld ", pivotal_elements.size()));
        PRLEVEL(p, (" and now is %ld\n ", ii));
        pivotal_elements.resize(ii);
    }

    /************ Making the heap from list of the immediate children
     * ******/
    PRLEVEL(1, ("%% Next: work on the heap \n"));
    
    ParU_ResultCode m_h =
    paru_make_heap(f, start_fac, pivotal_elements, hi, colHash, paruMatInfo);
    if (m_h != PARU_SUCCESS) //some problem happend in make_heap
        return m_h;
    PRLEVEL(1, ("%% Done: work on the heap \n"));

    Int eli = snM[f];
    std::vector<Int> **heapList = paruMatInfo->heapList;
    std::vector<Int> *curHeap = heapList[eli];

    if (curHeap->empty()) return PARU_SUCCESS; 

#ifndef NDEBUG
    p = 1;
#endif

#ifndef NDEBUG
    Int *lacList = paruMatInfo->lacList;
    PRLEVEL(p, ("%% current heap:\n %%"));
    for (Int k = 0; k < (Int)curHeap->size(); k++)
    {
        Int ee = (*curHeap)[k];
        paru_Element *ell = elementList[ee];
        PRLEVEL(p, ("%ld-%ld", k, ee));
        if (ell != NULL)
        {
            PRLEVEL(p, ("(%ld) ", lacList[ee]));
        }
        else
        {
            PRLEVEL(p, ("(*%ld) ", lacList[ee]));
        }
    }
    PRLEVEL(p, ("\n"));
#endif

#ifndef NDEBUG
    // chekcing the heap
    for (Int i = curHeap->size() - 1; i > 0; i--)
    {
        Int elid = (*curHeap)[i];
        Int pelid = (*curHeap)[(i - 1) / 2];  // parent id
        if (lacList[pelid] > lacList[elid])
        {
            PRLEVEL(p, ("%ld-%ld(%ld) <", (i - 1) / 2, pelid, lacList[pelid]));
            PRLEVEL(p, ("%ld-%ld(%ld) \n", i, elid, lacList[elid]));
        }
        ASSERT(lacList[pelid] <= lacList[elid]);
    }

#endif
    return PARU_SUCCESS;
}
