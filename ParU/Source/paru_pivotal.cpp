////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_pivotal ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief
 *  adding the list of pivotal elements from the heap, computing the list of
 *  rows and assembling pivotal columns
 *
 * @param pivotal_elements list
 *          f and paruMatInfo
 *
 *  @author Aznaveh
 */
#include "paru_internal.hpp"

ParU_Res paru_pivotal(std::vector<Int> &pivotal_elements,
                             std::vector<Int> &panel_row, Int &zero_piv_rows,
                             Int f, heaps_info &hi, paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;
    ParU_symbolic *Sym = paruMatInfo->Sym;
    Int *snM = Sym->super2atree;
    std::vector<Int> **heapList = paruMatInfo->heapList;
    Int eli = snM[f];

    Int *Super = Sym->Super;
    Int col1 = Super[f]; /* fornt F has columns col1:col2-1 */
    Int col2 = Super[f + 1];
    Int *aChild = Sym->aChild;
    Int *aChildp = Sym->aChildp;

#ifndef NDEBUG
    Int m = paruMatInfo->m;
    PRLEVEL(PR, ("%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"));
    PRLEVEL(PR, ("%% Pivotal assembly of front %ld (eli %ld) cols %ld-%ld\n", f,
                eli, col1, col2));
#endif

    ParU_Element **elementList = paruMatInfo->elementList;

    Int *lacList = paruMatInfo->lacList;

    work_struct *Work = paruMatInfo->Work;
    Int *rowMarkp = Work->rowMark;
    Int rowMark = 0;

    /*****  making the list of elements that contribute to pivotal columns ****/
    Int biggest_Child_id = -1;
    Int biggest_Child_size = -1;
    Int tot_size = 0;

    for (Int i = aChildp[eli]; i <= aChildp[eli + 1] - 1; i++)
    {
        Int chelid = aChild[i];  // element id of the child
        // max(rowMark , child->rowMark)
        Int f_rmark = rowMarkp[chelid];
        rowMark = rowMark > f_rmark ? rowMark : f_rmark;

        PRLEVEL(PR, ("%% chelid = %ld\n", chelid));
        std::vector<Int> *curHeap = heapList[chelid];

        if (curHeap == nullptr) continue;

        while (curHeap->size() > 0)
        // pop from the heap and put it in pivotal_elements
        {
            Int frontEl = curHeap->front();
            Int lacFel = lacList[frontEl];
            PRLEVEL(PR, ("%% element = %ld col1=%ld", frontEl, col1));
            PRLEVEL(PR, (" lac_el = %ld \n", lacFel));
            // ASSERT(lacFel >= col1);
            PRLEVEL(PR, ("%% curHeap->size= %ld \n", curHeap->size()));

            if (lacFel >= col2) break;

            if (elementList[frontEl] != NULL)
                pivotal_elements.push_back(frontEl);
            std::pop_heap(
                curHeap->begin(), curHeap->end(),
                [&lacList](Int a, Int b) { return lacList[a] > lacList[b]; });
            curHeap->pop_back();
        }

        if (curHeap->empty())
        {
            delete heapList[chelid];
            heapList[chelid] = nullptr;
        }
        else
        {
            // IMPORTANT: type conversion is necessary
            Int cur_size = curHeap->size();
            PRLEVEL(PR, ("%% curHeap->size= *%ld \n", curHeap->size()));
            PRLEVEL(PR, ("%% biggest_Child_size = %ld \n", biggest_Child_size));
            tot_size += curHeap->size();
            if (cur_size > biggest_Child_size)
            {
                PRLEVEL(PR, ("%% biggest_Child_id = %ld \n", biggest_Child_id));
                biggest_Child_id = chelid;
                biggest_Child_size = cur_size;
            }
        }
    }

    hi.sum_size = tot_size;
    hi.biggest_Child_id = biggest_Child_id;
    hi.biggest_Child_size = biggest_Child_size;

    PRLEVEL(PR, ("%%Inside pivot tot_size= %ld \n", hi.sum_size));
    PRLEVEL(PR, ("%% biggest_Child_id = %ld \n", hi.biggest_Child_id));
    PRLEVEL(PR, ("%% hi.biggest_Child_size = %ld \n", hi.biggest_Child_size));

    rowMarkp[eli] = rowMark;

    Int *isRowInFront = Work->rowSize;
    if (++rowMark < 0)
    // just look at the children
    {  // in rare case of overflow
        Int *Sleft = Sym->Sleft;
        // first column of current front until first column of next front
        for (Int i = Sleft[col1]; i < Sleft[Super[f + 1]]; i++)
            isRowInFront[i] = -1;
        rowMark = rowMarkp[eli] = 1;
    }
    rowMarkp[eli] = rowMark;
    PRLEVEL(1, ("%% rowMark=%ld;\n", rowMark));

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%% pivotal columns eli(%ld): ", eli));
    for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
        PRLEVEL(PR, ("%ld ", pivotal_elements[i]));
    PRLEVEL(PR, ("\n"));
    std::set<Int> stl_rowSet;
    std::set<Int>::iterator it;
#endif
    Int panel_width = paruMatInfo->panel_width;
    Int fp = col2 - col1; /* first fp columns are pivotal */
    Int num_panels = (Int)ceil((double)fp / panel_width);

    Int *frowList = paruMatInfo->frowList[f];
    Int rowCount = 0;
    // Int zero_piv_rows = 0;

    /*************** finding set of rows in current front *********************/
    for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        ParU_Element *el = elementList[e];
        ASSERT(el != NULL);

#ifndef NDEBUG
        // Int *el_colIndex = colIndex_pointer (curEl);
        Int *el_colIndex = (Int *)(el + 1);
        PRLEVEL(PR, ("current element(%ld) %ld-%ld\n", e, col1, col2));
        PRLEVEL(PR, ("lac = %ld ", el->lac));
        PRLEVEL(PR, ("lac_col = %ld\n ", lacList[e]));
        ASSERT(el_colIndex[el->lac] >= col1);
        if (PR <= 0) paru_print_element(paruMatInfo, e);
#endif

        Int mEl = el->nrows;
        Int nEl = el->ncols;

        // Int *el_rowIndex = rowIndex_pointer (el);
        Int *el_rowIndex = (Int *)(el + 1) + nEl;

        // Int *rowRelIndex = relRowInd (el);
        Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

        PRLEVEL(1, ("%% rowMark=%ld;\n", rowMark));

        el->nzr_pc = 0;  // initializing ; number of zero rows
        Int nrows2bSeen = el->nrowsleft;

        for (Int rEl = 0; rEl < mEl; rEl++)
        {
            Int curRow = el_rowIndex[rEl];
            PRLEVEL(1, ("%% curRow =%ld rEl=%ld\n", curRow, rEl));
            if (nrows2bSeen == 0) break;
            if (curRow < 0) continue;  // that row has already deleted
            nrows2bSeen--;

#ifndef NDEBUG
            //            stl_rowSet.insert(curRow);
            PRLEVEL(1, ("%% %p ---> isRowInFront [%ld]=%ld\n",
                        isRowInFront + curRow, curRow, isRowInFront[curRow]));
#endif

            if (isRowInFront[curRow] < rowMark)
            {  // first time seeing curRow
                // Int *el_colIndex = colIndex_pointer (curEl);
                Int *el_colIndex = (Int *)(el + 1);

                //#if 1
                // checkikng if the numerical values are hard zero
                // look at the pivotal columns and check if there is any
                // nonzeros if there is none I can skip adding this row
                bool nz_found = false;
                for (Int cEl = el->lac; cEl < nEl; cEl++)
                {
                    if (el_colIndex[cEl] < 0)  // already assembled somewhere
                        continue;
                    if (el_colIndex[cEl] >= col2) break;

                    // double *el_Num = numeric_pointer (el);
                    double *el_Num =
                        (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);
                    if (el_Num[cEl * mEl + rEl] != 0.0)
                    {
                        nz_found = true;
                        break;
                    }
                }
                if (!nz_found)
                {
                    el->nzr_pc++;
                    PRLEVEL(1, ("%% Found a row with all zeroes!! "
                                "curRow =%ld el=%ld\n",
                                curRow, e));

                    zero_piv_rows++;
                    rowRelIndex[rEl] = -1;
#ifndef NDEBUG
                    if (PR <= 0) paru_print_element(paruMatInfo, e);
#endif
                    continue;  // Not adding the row
                }
                //#endif
                // Adding curRow to the set
#ifndef NDEBUG
                stl_rowSet.insert(curRow);
#endif
                PRLEVEL(1, ("%%curRow =%ld rowCount=%ld\n", curRow, rowCount));
                frowList[rowCount] = curRow;
                rowRelIndex[rEl] = rowCount;
                PRLEVEL(1, ("%%1st: rowRelIndex[%ld] = %ld\n", rEl, rowCount));
                isRowInFront[curRow] = rowMark + rowCount++;
            }
            else
            {  // already seen curRow
                PRLEVEL(1, ("%%curRow =%ld rowCount=%ld\n", curRow, rowCount));
                PRLEVEL(1, ("%%before updating rowRelIndex[%ld] = %ld\n", rEl,
                            rowRelIndex[rEl]));
                PRLEVEL(1, ("%% rowMark =%ld\n", rowMark));
                rowRelIndex[rEl] = isRowInFront[curRow] - rowMark;
                PRLEVEL(1, ("%%N1st: rowRelIndex[%ld] = %ld\n", rEl,
                            rowRelIndex[rEl]));
            }

            ASSERT(rowCount <= m);
#ifndef NDEBUG
            if (rowCount != (Int)stl_rowSet.size())
            {
                PRLEVEL(1, ("%%curRow =%ld rowCount=%ld\n", curRow, rowCount));
                PRLEVEL(1, ("%%stl_rowSet.size()=%ld \n", stl_rowSet.size()));
            }
#endif
            ASSERT(rowCount == (Int)stl_rowSet.size());
        }
        panel_row[(lacList[e] - col1) / panel_width] = rowCount;
#ifndef NDEBUG
        PRLEVEL(PR, ("%%rowCount=%ld", rowCount));
        PRLEVEL(PR, (" lac=%ld", lacList[e]));
        ASSERT((lacList[e] - col1) / panel_width < num_panels);
        ASSERT((lacList[e] - col1) / panel_width >= 0);
        PRLEVEL(PR, (" ind.=%ld\n", (lacList[e] - col1) / panel_width));
#endif
    }

    if (rowCount < fp)
    {
#ifndef NDEBUG
        // there is a structural problem
        PRLEVEL(PR,
                ("%%STRUCTURAL PROBLEM! rowCount=%ld, fp =%ld", rowCount, fp));
#endif
        if (rowCount + zero_piv_rows > fp)
        {
            PRLEVEL(PR,
                    (" it can be solved by adding %ld zeros", zero_piv_rows));
        }
        else
        {
            PRLEVEL(PR, (" it wil FAIL"));
        }
#ifndef NDEBUG
        PRLEVEL(PR, ("\n"));
#endif
    }

    // make sure that all panel_row is correctly initialized
    PRLEVEL(PR, ("%% num_panels: %ld \n ", num_panels));
    PRLEVEL(PR, ("%% panel_row: \n %%"));
    Int pprow = panel_row[0];
    PRLEVEL(PR, ("%% %ld ", pprow));
    ASSERT(pprow != 0);
    for (Int i = 1; i < num_panels; i++)
    {
        if (pprow > panel_row[i])
        {
            panel_row[i] = pprow;
        }
        else
        {
            pprow = panel_row[i];
        }
        PRLEVEL(1, ("%ld ", panel_row[i]));
        ASSERT(panel_row[i] > 0);
        ASSERT(panel_row[i] <= m);
    }

    paruMatInfo->frowCount[f] = rowCount;

#ifndef NDEBUG /* Checking if pivotal rows are correct */
    PRLEVEL(PR, ("%% panel_row: \n %%"));
    for (Int i = 0; i < num_panels; i++) PRLEVEL(PR, ("%ld ", panel_row[i]));
    PRLEVEL(PR, ("\n"));
    PRLEVEL(PR, ("%%There are %ld rows x %ld columns %ld - %ld "
                "in front %ld with %ld zero rows: \n%%",
                rowCount, fp, col1, col2, f, zero_piv_rows));
    for (Int i = 0; i < rowCount; i++) PRLEVEL(PR, (" %ld", frowList[i]));
    PRLEVEL(PR, ("\n"));
    Int stl_rowSize = stl_rowSet.size();
    if (rowCount != stl_rowSize)
    {
        PRLEVEL(PR, ("%% STL %ld:\n", stl_rowSize));
        for (it = stl_rowSet.begin(); it != stl_rowSet.end(); it++)
            PRLEVEL(PR, ("%% %ld", *it));
        PRLEVEL(PR, ("\n%%My Set %ld:\n", rowCount));
        for (Int i = 0; i < rowCount; i++) PRLEVEL(PR, ("%% %ld", frowList[i]));
        PRLEVEL(PR, ("\n"));
    }
    ASSERT(rowCount == stl_rowSize);
#endif

    Int fm = Sym->Fm[f]; /* Upper bound number of rows of F */
    ASSERT(fm >= rowCount);

    // freeing extra space for rows
    if (rowCount != fm)
    {
        size_t sz = sizeof(Int) * fm;
        frowList = (Int *)paru_realloc(rowCount, sizeof(Int), frowList, &sz);
        paruMatInfo->frowList[f] = frowList;
    }

    double *pivotalFront = (double *)paru_calloc(rowCount * fp, sizeof(double));

    if (pivotalFront == NULL || frowList == NULL )
    {
        printf(
            "Paru: 0ut of memory when tried to allocate for pivotal part %ld\n",
            f);
        return PARU_OUT_OF_MEMORY;
    }

#ifndef NDEBUG
    paruMatInfo->actual_alloc_LUs += rowCount * fp;
    paruMatInfo->actual_alloc_row_int += rowCount;
    if (fm != rowCount) PRLEVEL(PR, ("%% fm=%ld rowCount=%ld ", fm, rowCount));
    PRLEVEL(PR, ("%% LUs=%ld ", paruMatInfo->actual_alloc_LUs));
    PRLEVEL(PR, ("%% pivotalFront = %p size=%ld", pivotalFront, rowCount * fp));
    Int act = paruMatInfo->actual_alloc_LUs + paruMatInfo->actual_alloc_Us +
              paruMatInfo->actual_alloc_row_int;
    Int upp = Sym->Us_bound_size + Sym->LUs_bound_size +
              Sym->row_Int_bound + Sym->col_Int_bound;
    PRLEVEL(PR, ("%% MEM=%ld percent=%lf%%", act, 100.0 * act / upp));
    PRLEVEL(PR, ("%% MEM=%ld percent=%lf%%\n", act, 100.0 * act / upp));
#endif
    paru_fac *LUs = paruMatInfo->partial_LUs;
    paruMatInfo->frowCount[f] = rowCount;

    LUs[f].m = rowCount;
    LUs[f].n = fp;
    ASSERT(LUs[f].p == NULL);
    LUs[f].p = pivotalFront;

    /***************  assembling the pivotal part of the front ****************/
    /*
     *
     *  el           nEl
     *              6, 7, 11, 12
     *             _____________
     *          23 | X  Y  .  .     stored in memory like this:
     *      mEl 17 | X  Y  .  .     ..6, 7,11, 12, 23, 17, 2, X, X, X, Y, Y, Y,
     *           2 | X  Y  .  .
     *
     *    It must be assembled in current pivotal fron like this:
     *                                     fp
     *                                 col1, ... , col
     *
     *                                  6, 7, 8, 9, 10
     *                                  ______________
     *                          0   23 | X  Y  .  .  .
     *               rowCount   1    2 | X  Y  .  .  .
     *                          2    4 | *  *  .  .  .  isRowInFront[4] == 2
     *                          3   17 | X  Y  .  .  .
     * */

    Int ii = 0;  // using for resizing pivotal_elements
    for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        paru_full_summed(e, f, paruMatInfo);
        if (elementList[e] != NULL)
        {  // keeping the element
            pivotal_elements[ii++] = pivotal_elements[i];
        }
    }

   if (ii < (Int)pivotal_elements.size())
    {
        PRLEVEL(PR, ("%% pivotal size was %ld ", pivotal_elements.size()));
        PRLEVEL(PR, ("%% and now is %ld\n ", ii));
        pivotal_elements.resize(ii);
    }

    // second pass through elements with zero rows and
    // growing them to better fit in current front
    // This can possibly help in more assembly

    Int num_children_with0 = 0;
    Int num_children_with0_which_fit = 0;

    for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        ParU_Element *el = elementList[e];
        if (el->nzr_pc > 0)  // an elemen that has at least one zero row
        {
            num_children_with0++;
            Int mEl = el->nrows;
            Int nEl = el->ncols;

            // Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_rowIndex = (Int *)(el + 1) + nEl;

            // Int *rowRelIndex = relRowInd (el);
            Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;
            Int nrows2bSeen = el->nrowsleft;

            for (Int rEl = 0; rEl < mEl; rEl++)
            {
                Int curRow = el_rowIndex[rEl];
                if (nrows2bSeen == 0) break;
                if (curRow < 0) continue;  // that row has already deleted
                nrows2bSeen--;
                if (rowRelIndex[rEl] == -1)  // the zero row
                {
                    if (isRowInFront[curRow] >= rowMark)
                    {
                        el->nzr_pc--;
                        rowRelIndex[rEl] = isRowInFront[curRow] - rowMark;
                    }
                }
            }
#ifndef NDEBUG
            if (el->nzr_pc == 0)
            {  // all the zero rows fit in the front
                PRLEVEL(1, ("%%element %ld totally fit in current front %ld\n",
                            e, f));
                num_children_with0_which_fit++;
            }
            ASSERT(el->nzr_pc >= 0);
#endif
        }
    }

    if (num_children_with0 == num_children_with0_which_fit)
    {  // all the children fit within current front
        zero_piv_rows = 0;
    }

#ifndef NDEBUG
    PRLEVEL(PR, ("%% pivotal columns eli(%ld) after resizing: ", eli));
    for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
        PRLEVEL(PR, ("%ld ", pivotal_elements[i]));
    PRLEVEL(PR, ("\n"));

    PR = 2;
    PRLEVEL(PR, ("%% After all the assemble %ld, z=%ld\n", f, zero_piv_rows));
    PRLEVEL(PR, ("%% x =  \t"));
    for (Int c = col1; c < col2; c++) PRLEVEL(PR, ("%ld\t\t", c));
    PRLEVEL(PR, (" ;\n"));
    for (Int r = 0; r < rowCount; r++)
    {
        PRLEVEL(PR, ("%% %ld\t", frowList[r]));
        for (Int c = col1; c < col2; c++)
            PRLEVEL(PR, (" %2.5lf\t", pivotalFront[(c - col1) * rowCount + r]));
        PRLEVEL(PR, ("\n"));
    }
    PR = 2;
    PRLEVEL(PR, ("x%ld = [ \t", f));
    for (Int r = 0; r < rowCount; r++)
    {
        for (Int c = col1; c < col2; c++)
            PRLEVEL(PR, (" %2.5lf\t", pivotalFront[(c - col1) * rowCount + r]));
        PRLEVEL(PR, ("\n"));
    }
    PRLEVEL(PR, (" ]; %% %ld*%ld\n", rowCount, fp));
    PR = 1;
#endif

    rowMarkp[eli] += rowCount;
    PRLEVEL(1, ("%% rowMarkp[%ld] =%ld\n", eli, rowMarkp[eli]));
    return PARU_SUCCESS;
}
