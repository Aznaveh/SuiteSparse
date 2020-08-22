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
        Int *panel_row, Int f)
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
    work_struct *Work =  paruMatInfo->Work;
    Int *rowMarkp = Work->rowMark;
    Int rowMark = rowMarkp[eli];


    PRLEVEL (p, ("%% col2= %ld\n", col2));
    while ( lnc_el(elementList, elHeap->front()) < col2 && elHeap->size() > 0)
        // pop from the heap and put it in pivotal_elements
    {
        Int front = elHeap->front(); 
        PRLEVEL (p, ("%% front = %ld ", front));
        PRLEVEL (p, (" lnc_el = %ld ", 
                    lnc_el(elementList, front)));
        PRLEVEL (p, ("%% elHeap->size= %ld \n", elHeap->size()));
        // max(rowMark , fornt->rowMark)
        Int f_rmark = rowMarkp[front];
        rowMark = rowMark >  f_rmark ?  rowMark : f_rmark;

        pivotal_elements.push_back(elHeap->front());
        std::pop_heap
            (elHeap->begin(), elHeap->end(),[&elementList](Int a, Int b)
             { return lnc_el(elementList,a) > lnc_el(elementList,b); }   );
        elHeap->pop_back();
    }

    rowMarkp[eli] = ++ rowMark;

#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%% eli(%ld): ", eli));
    for(Int i=0 ; i < pivotal_elements.size(); i++)
        PRLEVEL (p, ("%ld ", pivotal_elements[i]));
    PRLEVEL (p, ("\n"));
#endif 
    Int *isRowInFront = Work->rowSize; 
    Int panel_width = paruMatInfo->panel_width;



    //TODO: fix it with lnc


    Int *frowList = paruMatInfo->frowList[f];




    for(Int i=0 ; i < pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        paru_Element *el = elementList[e];
#ifndef NDEBUG
        p = -1;
        if (el == NULL) continue;
#endif 
        PRLEVEL (p, ("current element(%ld) ", e ));
        PRLEVEL (p, ("lnc = %ld ",  el->lnc));
        PRLEVEL (p, ("col = %ld\n ", lnc_el(elementList, e) ));


        //Int *el_colIndex = colIndex_pointer (el);
        Int *el_colIndex = (Int*)(el+1);

        Int mEl = el->nrows;
        Int nEl = el->ncols;

        //Int *el_rowIndex = rowIndex_pointer (el); 
        Int *el_rowIndex = (Int*)(el+1)+nEl; 

        //Int *rowRelIndex = relRowInd (el);
        Int *rowRelIndex = (Int*)(el+1)+ 2*nEl + mEl;

        //Int *colRelIndex = relColInd (el);
        Int *colRelIndex =  (Int*)(el+1)+ nEl + mEl;
        //
        //        //TODO: change rowMark completely. Each front should have its own
        //        //private mark
        //        Int rowMark = Work->rowMark;
        //        PRLEVEL (-1, ("rowMark=%ld;\n", rowMark));
        //
        //        if (rowMark < 0) 
        //        {  // in rare case of overflow
        //            memset (isRowInFront, -1, m*sizeof(Int));
        //            rowMark = Work->rowMark = 0;
        //        }
        //
        //        for (Int rEl = 0; rEl < mEl; rEl++)
        //        {
        //            Int curRow = el_rowIndex [rEl]; 
        //            PRLEVEL (1, ("%%@@curRow =%ld rEl=%ld\n", curRow, rEl));
        //            if (curRow < 0 ) continue; // that row has already deleted
        //            ASSERT (curRow < m ) ;
        //#ifndef NDEBUG
        //            Int p = 1;
        //            if (p <= 0) paru_print_element (paruMatInfo, e);
        //            stl_rowSet.insert (curRow);
        //            PRLEVEL (1, ("%% %p ---> isRowInFront [%ld]=%ld\n", 
        //                        isRowInFront+curRow, curRow, isRowInFront[curRow]));
        //#endif
        //
        //            if (isRowInFront[curRow] < rowMark )
        //            {  
        //                // Adding curRow to the set
        //                PRLEVEL (1, ("%%curRow =%ld rowCount=%ld\n", 
        //                            curRow, rowCount));
        //                frowList [rowCount] = curRow;
        //                rowRelIndex [rEl] = rowCount ;
        //                PRLEVEL (1, ("%%1st: rowRelIndex[%ld] = %ld\n",
        //                            rEl, rowCount ));
        //                isRowInFront [curRow] = rowMark + rowCount++; 
        //            }
        //            else
        //            {
        //                rowRelIndex [rEl] = isRowInFront [curRow] - rowMark;
        //                PRLEVEL (1, ("%%N1st: rowRelIndex[%ld] = %ld\n",
        //                            rEl, rowCount ));
        //            }
        //            ASSERT (rowCount <= m); 
        //        }
    }

    return;
}
