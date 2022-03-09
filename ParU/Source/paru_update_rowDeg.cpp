////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_update_rowDeg /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  growing current front if necessary and update the row degree of
 *   current front for current panel.
 *
 *  @author Aznaveh
 */
#include "paru_internal.hpp"

void paru_update_rowDeg(Int panel_num, Int row_end, Int f, Int start_fac,
                        std::set<Int> &stl_colSet,
                        std::vector<Int> &pivotal_elements,
                        paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;
#ifndef NDEBUG
    Int n = paruMatInfo->n;
    static Int r1 = 0, r2 = 0, r3 = 0;
#endif
    PRLEVEL(1, ("%%-------ROW degree update of panel %ld of front %ld \n",
                panel_num, f));
    Int panel_width = paruMatInfo->panel_width;
    paru_Element **elementList = paruMatInfo->elementList;
    work_struct *Work = paruMatInfo->Work;

    Int *elRow = Work->elRow;
    Int *elCol = Work->elCol;

    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *Super = LUsym->Super;
    Int col1 = Super[f];  // fornt F has columns col1:col2-1
    Int col2 = Super[f + 1];
    Int fp = col2 - col1;  // first fp columns are pivotal

    Int pMark = start_fac;  // Mark for pivotal rows
    Int npMark =
        ++paruMatInfo->time_stamp[f];  // making all the markings invalid

    Int colCount = stl_colSet.size();

    Int j1 = panel_num * panel_width;  // panel starting column
    Int j2 = (j1 + panel_width < fp) ? j1 + panel_width : fp;

    Int rowCount = paruMatInfo->frowCount[f];
    Int *row_degree_bound = paruMatInfo->row_degree_bound;

    std::set<Int> stl_newColSet;  // the list of new columns

    /*************** finding set of non pivotal cols in current front *********/
    /*
     *    Mark seen elements with pMark
     *        if Marked already added to the list
     *        else all columns are added to the current front
     *
     *                <----------fp--------->
     *                        j1     j2              Update here
     *                         ^     ^                stl_newColSet colCount
     *                         |     | stl_colSet     ^ . . .   ^
     *                         |     |        \       |         |
     *             F           |     |         [QQQQQ|OOOOOOOOO|....
     *              \  ____..._|_  ____...__ _____________________________...
     * ^              |\      |     |       #  ^     |         |
     * |              | \     |     |       #  | old | added   |
     * |              |  \    |     |       #  | list|  columns|
     * |              |___\...|_____|__...__#  |     |         |
     * |   ^    j1--> |       |\* **|       #  fp  oooooooooooo|
     * |   |     H    |       |**\**|       #  |   ooooEloooooo|
     * | panel   E    |       |***\*|       #  |   oooooooooooo|
     * | width   R    |       |***\*|       #  |               |
     * |   |     E    |       |***\*|       #  |    00000000   |
     * |   v          |____...|________..._ #  |    00El0000   |
     * |        j2--> |       |     |       #  v    00000000   |           ...
     * rowCount       |==================================================
     * |              |       |     |       |
     * |              |       |     |       |
     * |              |       |row_end      |
     * |              .       .      .      .
     * |              .       .      .      .
     * |              .       .      .      .
     * v              |___....______________|
     *
     */
    Int *frowList = paruMatInfo->frowList[f];
    std::set<Int>::iterator it;

    tupleList *RowList = paruMatInfo->RowList;
    for (Int i = j1; i < j2; i++)
    {
        Int curFsRow = frowList[i];
#ifndef NDEBUG
        Int curFsRowIndex = i;  // current fully summed row index
        PRLEVEL(1, ("%% 4: curFsRowIndex = %ld\n", curFsRowIndex));
        PRLEVEL(1, ("%% curFsRow =%ld\n", curFsRow));
#endif
        tupleList *curRowTupleList = &RowList[curFsRow];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT(numTuple >= 0);
        paru_Tuple *listRowTuples = curRowTupleList->list;
        PRLEVEL(1, ("%% 4: numTuple = %ld\n", numTuple));

        Int pdst = 0, psrc;
        for (psrc = 0; psrc < numTuple; psrc++)
        {
            paru_Tuple curTpl = listRowTuples[psrc];

            Int e = curTpl.e;
            Int curRowIndex = curTpl.f;
#ifndef NDEBUG
            r1++;
#endif
            if (e < 0 || curRowIndex < 0) continue;

            paru_Element *el = elementList[e];

            if (el == NULL) continue;

            Int nEl = el->ncols;
            // Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_rowIndex = (Int *)(el + 1) + nEl;
            if (el_rowIndex[curRowIndex] < 0) continue;

            Int mEl = el->nrows;
            // Int *rowRelIndex = relRowInd (el);
            Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;
            rowRelIndex[curTpl.f] = curFsRow;

            listRowTuples[pdst++] = curTpl;  // keeping the tuple

            ASSERT(el_rowIndex[curRowIndex] == curFsRow);

            if (el->rValid != pMark)
            {  // an element never seen before

                PRLEVEL(
                        1, ("%%P: first time seen elRow[%ld]=%ld \n", e, elRow[e]));
                PRLEVEL(1, ("%%pMark=%ld  npMark= %ld\n", pMark, npMark));

                // if (el->rValid < pMark)
                if (npMark == pMark + 1)
                    elRow[e] = el->nrowsleft - 1;  // initiaze
                el->rValid = pMark;

                PRLEVEL(1, ("%%changed to elRow[%ld]=%ld \n", e, elRow[e]));
#ifndef NDEBUG
                if (el->rValid > pMark)
                    PRLEVEL(1, ("%%pMark=%ld  rVal= %ld\n", pMark, el->rValid));
#endif
            }
            else  // el->rValid == pMark
            {     // already added to pivotal rows
                if (npMark == pMark + 1) elRow[e]--;
                PRLEVEL(1, ("%%already seen elRow[%ld]=%ld \n", e, elRow[e]));
                continue;
            }

            // Int *el_colIndex = colIndex_pointer (el);
            Int *el_colIndex = (Int *)(el + 1);

            PRLEVEL(1, ("%% element= %ld  nEl =%ld \n", e, nEl));
            for (Int cEl = 0; cEl < nEl; cEl++)
            {
                Int curCol = el_colIndex[cEl];
                PRLEVEL(1, ("%% curCol =%ld\n", curCol));
                ASSERT(curCol < n);

                if (curCol < 0)  // already deleted
                    continue;

                if (curCol < col2 && curCol >= col1) /*is a pivotal col */
                    continue;

                if (stl_colSet.find(curCol) == stl_colSet.end())
                {
                    stl_colSet.insert(curCol);
                    stl_newColSet.insert(curCol);
                    colCount++;
                }
#ifndef NDEBUG
                PR = 1;
                // stl_colSet.insert (curCol);
                for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
                    PRLEVEL(PR, ("%%@  %ld", *it));

#endif

                ASSERT(colCount <= n);
            }
        }

        curRowTupleList->numTuple = pdst;
    }

    Int *snM = LUsym->super2atree;
    Int eli = snM[f];
    if (colCount == 0)
    {  // there is no CB, Nothing to be done
        Work->rowMark[eli] += rowCount;
        return;
    }

#ifndef NDEBUG /* Checking if columns are correct */
    PR = 1;
    PRLEVEL(PR, ("%% There are %ld columns in this contribution block: \n",
                colCount));
    PRLEVEL(PR, ("\n"));
    Int stl_colSize = stl_colSet.size();

    if (colCount != stl_colSize)
    {
        PRLEVEL(PR, ("%% STL %ld:\n", stl_colSize));
        for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
            PRLEVEL(PR, ("%%  %ld", *it));
        PRLEVEL(PR, ("\n%% My Set %ld:\n", colCount));
        PRLEVEL(PR, ("\n"));
    }
    ASSERT(colCount == stl_colSize);
#endif

    // if the front did not grow, there is nothing else to do
    if (stl_newColSet.size() == 0) return;

    paruMatInfo->fcolCount[f] = colCount;

    /**** only travers over elements that contribute to pivotal columns *******/
    /*         to find their intersection
     *
     *         This can be fully skipped. I also can look other children in the
     *         heap
     *
     *            <----------fp--------->
     *                                          stl_newColSet
     *                             stl_colSet     ^        ^
     *                                    \       |        |
     *         F                           [QQQQQ|OOOOOOOOO|....
     *          \  ____..._________...__ _____________________________...
     * ^          |\      |     |       #  ^     |         |
     * |          | \     |     |       #  | old | added   |
     * |          |  \    |     |       #  | list|  columns|
     * |          |___\...|_____|__...__#  |     |         |
     * |   ^      |       |\* **|       #  fp  oooooooooooo|
     * |   |      |       |**\**|       #  |   ooo El ooooo|
     * | panel    |       |***\*|       #  |   oooooooooooo|
     * | width    |       |***\*|       #  |               |
     * |   |      |       |***\*|       #  |    00000000   |
     * |   v      |____...|________..._ #  |    000 El 0   |
     * |          |       |     |       #  v    00000000   |           ...
     * rowCount   |==================================================
     * |          |       |     |     ooooooooooooooooooooo
     * |          |       |     |     ooo HERE oooooooooooo
     * |          |       |row_end    oooooooo EL ooooooooo
     * |          .       .      .    ooooooooooooooooooooo
     * |          .       .      .      .
     * |          .       .      .      .      xxxxxxxxxx
     * v          |___....______________|      xxx EL xxx
     *                                         xxxxxxxxxx
     *                                         xxxxxxxxxx
     *
     *                                     oooooooooooo
     *                                     ooo El ooooo
     *                                     oooooooooooo
     *
     */
#ifndef NDEBUG
    PR = 1;
#endif
    PRLEVEL(1, ("%%Inside pivotal_elements\n"));
    Int ii = 0;
    for (Int i = 0; i < (Int)pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        paru_Element *el = elementList[e];
        if (el == NULL)
        {  // removing the  element from the list
            PRLEVEL(1, ("%% eli = %ld, element= %ld  \n", eli, e));
            continue;
        }
        // keeping other elements inside the list
        pivotal_elements[ii++] = pivotal_elements[i];

#ifndef NDEBUG
        PRLEVEL(PR, ("%% pivotal element= %ld lac=%ld colsleft=%ld \n", e,
                    el->lac, el->ncolsleft));
        if (PR <= 0) paru_print_element(paruMatInfo, e);
#endif
        Int intsct = paru_intersection(e, elementList, stl_newColSet);
        if (el->cValid < pMark)
        {  // first time seeing this element in this front
            elCol[e] = el->ncolsleft - intsct;  // initiaze
        }
        else if (el->cValid != npMark)
        {  // it has been seen
            elCol[e] -= intsct;
        }

        el->cValid = npMark;
    }

#ifndef NDEBUG
    PR = 1;
#endif

    /****** travers over new non pivotal rows of current front *****************
     *          and update the number of rows can be seen in each element.
     *
     *                <----------fp--------->
     *                        j1     j2
     *                         ^     ^
     *                         |     | stl_colSet
     *                         |     |        \       stl_newColSet
     *             F           |     |         [QQQQQ|OOOOOOOOO|....
     *              \  ____..._|_  ____...__ ____________________________...
     * ^              |\      |     |       #  ^     |         |
     * |              | \     |     |       #  | old | added   |
     * |              |  \    |     |       #  | list|  columns|
     * |              |___\...|_____|__...__#  |     |         |
     * |   ^          |       |\* **|       #  fp              |
     * |   |          |       |**\**|       #  |               |
     * | panel        |       |***\*|       #  |               |
     * | width        |       |***\*|       #  |               |
     * |   |          |       |***\*|       #  |               |
     * |   v          |____...|________..._ #  |      vvvvvv   |
     * |        j2--> |       |     |       #  v    00000000   |         ...
     * rowCount       |=============================00 El 00=============
     * |              |       |     |       |       00000000
     * |              |       |     |       |          vvvvvvvvv
     * |          H   |       |     |       |          xxxxxxxxxxxxxxxxxxxx
     * |          E   |       |row_end      |          xxxxxx El xxxxxxxxxx
     * |          R   .       .      .      .          xxxxxxxxxxxxxxxxxxxx
     * |          E   .       .      .      .
     * |              .       .      .      .
     * v  rowCount--->|___....______________|
     *
     */

    if (npMark == pMark + 1)  // just once for each front
    {                         // in the first time calling this function
        PRLEVEL(1, ("UPDATING elRow\n"));
        for (Int k = j2; k < rowCount; k++)
        {
            Int r = frowList[k];
            tupleList *curRowTupleList = &RowList[r];
            Int numTuple = curRowTupleList->numTuple;
            ASSERT(numTuple >= 0);
            paru_Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG
            Int PR = 1;
            PRLEVEL(PR, ("\n %%----r =%ld  numTuple = %ld\n", r, numTuple));
            if (PR <= 0) paru_print_tupleList(RowList, r);
#endif
            Int pdst = 0, psrc;
            for (psrc = 0; psrc < numTuple; psrc++)
            {
                paru_Tuple curTpl = listRowTuples[psrc];
                Int e = curTpl.e;

#ifndef NDEBUG
                if (PR <= 0) paru_print_element(paruMatInfo, e);
#endif
                Int curRowIndex = curTpl.f;

                if (e < 0 || curRowIndex < 0) continue;

                paru_Element *el = elementList[e];
                if (el == NULL) continue;

                // Int *el_rowIndex = rowIndex_pointer (el);
                Int *el_rowIndex = (Int *)(el + 1) + el->ncols;

                if (el_rowIndex[curRowIndex] < 0) continue;

                listRowTuples[pdst++] = curTpl;  // keeping the tuple

                if (el->rValid == pMark)
                {  // already a pivot and wont change the row degree
                    elRow[e]--;
                    PRLEVEL(1, ("%% Pivotal elRow[%ld]=%ld \n", e, elRow[e]));
                }
                else if (el->rValid != npMark)
                {
                    el->rValid = npMark;
                    elRow[e] = el->nrowsleft - 1;  // initiaze
                    PRLEVEL(1, ("%%rValid=%ld \n", el->rValid));
                    PRLEVEL(1, ("%%NP: first time seen elRow[%ld]=%ld \n", e,
                                elRow[e]));
                }
                else
                {  // el->rValid == npMark //it has been seen in this stage
                    elRow[e]--;
                    PRLEVEL(1, ("%%seen before: elRow[e]=%ld \n", elRow[e]));
                }
            }

            curRowTupleList->numTuple = pdst;
        }
    }

    /**** Travers over new non pivotal rows of current panel
     *    and update the row degree. Once at the begining it updates elRow and
     *    then if elRow == 0 it compute the intersection. It might have been
     *    computed or we might need to compute it now. However if it is
     *    computed now it would be mark not to recompute
     *
     *                <----------fp--------->
     *                        j1     j2
     *                         ^     ^
     *                         |     | stl_colSet
     *                         |     |        \       stl_newColSet
     *             F           |     |         [QQQQQ|OOOOOOOOO|....
     *              \  ____..._|_  ____...__ ____________________________...
     * ^              |\      |     |       #  ^     |         |
     * |              | \     |     |       #  | old | added   |
     * |              |  \    |     |       #  | list|  columns|
     * |              |___\...|_____|__...__#  |     |         |
     * |   ^          |       |\* **|       #  fp              |
     * |   |          |       |**\**|       #  |               |
     * | panel        |       |***\*|       #  |               |
     * | width        |       |***\*|       #  |               |
     * |   |          |       |***\*|       #  |               |
     * |   v          |____...|________..._ #  |      vvvvvv   |
     * |        j2--> |       |     |       #  v    00000000   |         ...
     * rowCount   H   |=============================00 El 00=============
     * |          E   |       |     |       |       00000000
     * |          R   |       |     |       |          vvvvvvvvv
     * |          E   |       |     |       |          xxxxxxxxxxxxxxxxxxxx
     * |  row_end --->|       |row_end      |          xxxxxx El xxxxxxxxxx
     * |              .       .      .      .          xxxxxxxxxxxxxxxxxxxx
     * |              .       .      .      .
     * |              .       .      .      .
     * v  rowCount--->|___....______________|
     *
     */

    Int new_row_degree_bound_for_r;

    for (Int k = j2; k < row_end; k++)
    {
        Int r = frowList[k];

        new_row_degree_bound_for_r = colCount;

        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT(numTuple >= 0);
        paru_Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG
        Int PR = 1;
        PRLEVEL(PR,
                ("\n %%--------> 2nd r =%ld  numTuple = %ld\n", r, numTuple));
        if (PR <= 0) paru_print_tupleList(RowList, r);
#endif

        Int pdst = 0, psrc;
        for (psrc = 0; psrc < numTuple; psrc++)
        {
            paru_Tuple curTpl = listRowTuples[psrc];
            Int e = curTpl.e;

#ifndef NDEBUG
            if (PR <= 0) paru_print_element(paruMatInfo, e);
#endif
            Int curRowIndex = curTpl.f;

            paru_Element *el = elementList[e];
            // ASSERT (el != NULL);
            if (el == NULL) continue;

            // Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_rowIndex = (Int *)(el + 1) + el->ncols;

            if (el_rowIndex[curRowIndex] < 0) continue;

            listRowTuples[pdst++] = curTpl;  // keeping the tuple

            if (el->rValid == pMark)
            {  // already a pivot and wont change the row degree
                PRLEVEL(1, ("%% Pivotal elRow[%ld]=%ld \n", e, elRow[e]));
                continue;
            }

            if (elRow[e] != 0)
            {                            // use the upperbound
                if (el->cValid < pMark)  // never seen
                    new_row_degree_bound_for_r += el->ncolsleft;
                else  // tighter upperbound
                    new_row_degree_bound_for_r += elCol[e];
                continue;
            }

            if (el->cValid < pMark)
            {  // first time seeing this element in this front
                el->cValid = npMark;
                Int intsct = paru_intersection(e, elementList, stl_newColSet);
                elCol[e] = el->ncolsleft - intsct;  // initiaze
            }
            else if (el->cValid != npMark)
            {  // it has been seen
                el->cValid = npMark;
                Int intsct;
                if (elCol[e] != 0)
                {
                    intsct = paru_intersection(e, elementList, stl_newColSet);
                    elCol[e] -= intsct;
                }
            }
            new_row_degree_bound_for_r += elCol[e];

            PRLEVEL(1, ("%% pMark=%ld npMark=%ld \n", pMark, npMark));
        }
        curRowTupleList->numTuple = pdst;

        Int old_bound_updated = row_degree_bound[r] + colCount - 1;

#ifndef NDEBUG
        PR = 1;
        PRLEVEL(PR, ("%%old_bound_updated =%ld \n", old_bound_updated));
        PRLEVEL(PR, ("%%new_row_degree_bound_for_r=%ld \n",
                    new_row_degree_bound_for_r));
        PRLEVEL(PR, ("%%row_degroo_bound[%ld]=%ld \n", r, row_degree_bound[r]));
#endif

        row_degree_bound[r] =  // min
            old_bound_updated < new_row_degree_bound_for_r
                ? old_bound_updated
                : new_row_degree_bound_for_r;
    }

    paruMatInfo->time_stamp[f] += 2;  // making all the markings invalid again
#ifndef NDEBUG
    PRLEVEL(1, ("%% Finalized counters r1=%ld r2=%ld r3=%ld sum=%ld\n", r1, r2,
                r3, r1 + r2 + r3));
#endif
}
