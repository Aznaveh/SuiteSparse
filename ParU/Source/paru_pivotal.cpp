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
        Int *panel_row, Int f, cholmod_common *cc)
{
    DEBUGLEVEL(0);
    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    std::vector<Int>** heapList = paruMatInfo->heapList;
    Int eli = snM [f]; 

#ifndef NDEBUG

    Int p = 1;
    PRLEVEL (p, ("%% Pivotal assembly of front %ld (eli %ld)\n",f, eli));
#endif 


    Int *Super = LUsym->Super;
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];


    std::vector<Int>* elHeap = heapList[eli] ;
    paru_Element **elementList = paruMatInfo->elementList;
    
    Int m = paruMatInfo-> m;
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%% Before making the pivotal vector:\n %%"));
    PRLEVEL (p, ("%% element ids:\n %%"));
    for(Int i = 0; i < elHeap->size(); i++)
    {
        Int elid = (*elHeap)[i];
        PRLEVEL (p, (" %ld", elid));
        PRLEVEL (p, (" (%ld) ", lnc_el(elementList, elid) ));
    }
    PRLEVEL (p, ("\n"));
    Int priorEl_lnc = 0;
#endif 

    /*****  making the list of elements that contribute to pivotal columns ****/
    while (elHeap->size() > 0) 
    // pop from the heap and put it in pivotal_elements
    {
        Int frontEl = elHeap->front(); 
        Int lncFel = lnc_el(elementList, frontEl);
        PRLEVEL (p, ("%% element = %ld col1=%ld", frontEl, col1));
        PRLEVEL (p, (" lnc_el = %ld \n", lncFel));
        ASSERT (lncFel >= col1);
        PRLEVEL (p, ("%% elHeap->size= %ld \n", elHeap->size()));

        if ( lncFel >= col2 ) break;

        if (elementList[frontEl] != NULL) 
            pivotal_elements.push_back(frontEl);
        std::pop_heap
            (elHeap->begin(), elHeap->end(),[&elementList](Int a, Int b)
             { return lnc_el(elementList,a) > lnc_el(elementList,b); }   );
        elHeap->pop_back();
#ifndef NDEBUG
        //ensure the vector is sorted based on first column
        if ( lncFel < priorEl_lnc )
            PRLEVEL (p, ("%% cur= %ld prior=%ld\n", lncFel, priorEl_lnc));
        ASSERT ( lncFel >= priorEl_lnc );
        priorEl_lnc = lncFel;
#endif 
    }

#ifndef NDEBUG
    //chekcing the heap
    if (elHeap->size() > 0 )
    {
        for(Int i = elHeap->size()-1 ; i > 0; i--)
        {
            Int elid = (*elHeap)[i];
            Int pelid = (*elHeap)[(i-1)/2]; //parent id
            ASSERT ( lnc_el(elementList,pelid) <= lnc_el(elementList,elid));
        }
    }
#endif


    work_struct *Work =  paruMatInfo->Work;
    Int *rowMarkp = Work->rowMark;
    Int rowMark = rowMarkp[eli];

    Int *isRowInFront = Work->rowSize; 
    if ( ++rowMark < 0) 
        // just look at the children
    {  // in rare case of overflow
        Int *Sleft = LUsym->Sleft;
        // first column of current front until first column of next front
        for (Int i = Sleft[col1]; i < Sleft[Super[f+1]]; i++)
            isRowInFront[i] = -1;
        rowMark = rowMarkp[eli] = 1;
    }
    rowMarkp[eli] = rowMark;
    PRLEVEL (1, ("%% rowMark=%ld;\n", rowMark));


#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%% pivotal columns eli(%ld): ", eli));
    for(Int i=0 ; i < pivotal_elements.size(); i++)
        PRLEVEL (p, ("%ld ", pivotal_elements[i]));
    PRLEVEL (p, ("\n"));
    std::set<Int> stl_rowSet;
    std::set<Int>::iterator it;
#endif 
    Int panel_width = paruMatInfo->panel_width;
    Int fp = col2 - col1;   /* first fp columns are pivotal */ 
    Int num_panels = (Int) ceil( (double)fp/panel_width);


    Int *frowList = paruMatInfo->frowList[f];
    Int rowCount = 0;

    /*************** finding set of rows in current front *********************/
    for(Int i = 0 ; i < pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        paru_Element *el = elementList[e];
        ASSERT (el != NULL);

        PRLEVEL (p, ("current element(%ld) ", e ));
        PRLEVEL (p, ("lnc = %ld ",  el->lnc));
        PRLEVEL (p, ("lnc_col = %ld\n ", lnc_el(elementList, e) ));

        Int mEl = el->nrows;
        Int nEl = el->ncols;

        //Int *el_rowIndex = rowIndex_pointer (el); 
        Int *el_rowIndex = (Int*)(el+1)+nEl; 

        //Int *rowRelIndex = relRowInd (el);
        Int *rowRelIndex = (Int*)(el+1)+ 2*nEl + mEl;


        PRLEVEL (1, ("%% rowMark=%ld;\n", rowMark));

        for (Int rEl = 0; rEl < mEl; rEl++)
        {
            Int curRow = el_rowIndex [rEl]; 
            PRLEVEL (1, ("%% curRow =%ld rEl=%ld\n", curRow, rEl));
            if (curRow < 0 ) continue; // that row has already deleted
            ASSERT (curRow < m ) ;
#ifndef NDEBUG
            Int p = 1;
            if (p <= 0) paru_print_element (paruMatInfo, e);
            stl_rowSet.insert (curRow);
            PRLEVEL (1, ("%% %p ---> isRowInFront [%ld]=%ld\n", 
                        isRowInFront+curRow, curRow, isRowInFront[curRow]));
#endif

            if (isRowInFront[curRow] < rowMark )
            {  //first time seeing curRow 
                // Adding curRow to the set
                PRLEVEL (1, ("%%curRow =%ld rowCount=%ld\n",curRow, rowCount));
                frowList [rowCount] = curRow;
                rowRelIndex [rEl] = rowCount ;
                PRLEVEL (1, ("%%1st: rowRelIndex[%ld] = %ld\n",
                            rEl, rowCount ));
                isRowInFront [curRow] = rowMark + rowCount++; 
            }
            else
            {//already seen curRow
                PRLEVEL (1, ("%%curRow =%ld rowCount=%ld\n",curRow, rowCount));
                PRLEVEL (1, ("%%before updating rowRelIndex[%ld] = %ld\n",
                            rEl, rowRelIndex[rEl]));
                PRLEVEL (1, ("%% rowMark =%ld\n",rowMark));
                rowRelIndex [rEl] = isRowInFront [curRow] - rowMark;
                PRLEVEL (1, ("%%N1st: rowRelIndex[%ld] = %ld\n",
                            rEl, rowRelIndex[rEl]));
            }

            ASSERT (rowCount <= m); 
#ifndef NDEBUG 
            if (rowCount != stl_rowSet.size())
            {
                PRLEVEL (1, ("%%curRow =%ld rowCount=%ld\n",curRow, rowCount));
                PRLEVEL (1, ("%%stl_rowSet.size()=%ld \n", stl_rowSet.size()));
            }
#endif 
            ASSERT (rowCount == stl_rowSet.size());
        }
        panel_row [( lnc_el(elementList,e) - col1) / panel_width] = rowCount;
#ifndef NDEBUG 
        p = 1;
        PRLEVEL (p, ("%%rowCount=%ld", rowCount));
        PRLEVEL (p, (" lnc=%ld", lnc_el(elementList,e)));
        ASSERT (( lnc_el(elementList,e) - col1) / panel_width < num_panels);
        PRLEVEL (p, (" ind.=%ld\n", 
                    (lnc_el(elementList,e) - col1) / panel_width));
#endif 
    }


    // make sure that all panel_row is correctly initialized
    PRLEVEL (p, ("%% num_panels: %ld \n ",num_panels));
    PRLEVEL (p, ("%% panel_row: \n %%"));
    Int pprow = panel_row[0];
    PRLEVEL (p, ("%% %ld ",pprow));
    ASSERT (pprow != 0);
    for (Int i = 1; i < num_panels; i++)
    {
        if (pprow > panel_row[i])
        {
            panel_row[i] = pprow ;
        }
        else
        {
            pprow = panel_row[i];
        }
        PRLEVEL (1, ("%ld ",panel_row[i]));
        ASSERT (panel_row[i] > 0);
        ASSERT (panel_row[i] <= m);

    }

    paruMatInfo->frowCount[f] = rowCount;

#ifndef NDEBUG /* Checking if pivotal rows are correct */
    p = 1;
    PRLEVEL (p, ("%% panel_row: \n %%"));
    for (Int i = 0; i < num_panels; i++)
        PRLEVEL (p, ("%ld ",panel_row[i]));
    PRLEVEL (p, ("\n"));
    PRLEVEL (p, ("%%There are %ld rows x %ld columns %ld - %ld "
                "in this front: \n %%", rowCount, fp, col1, col2));
    for (Int i = 0; i < rowCount; i++)
        PRLEVEL (p, (" %ld", frowList [i]));
    PRLEVEL (p, ("\n"));
    Int stl_rowSize = stl_rowSet.size();
    if (rowCount != stl_rowSize)
    {
        PRLEVEL (p, ("%% STL %ld:\n",stl_rowSize));
        for (it = stl_rowSet.begin(); it != stl_rowSet.end(); it++)
            PRLEVEL (p, ("%% %ld", *it));
        PRLEVEL (p, ("\n%%My Set %ld:\n",rowCount));
        for (Int i = 0; i < rowCount; i++)
            PRLEVEL (p, ("%% %ld", frowList [i]));
        PRLEVEL (p, ("\n"));
    }
    ASSERT (rowCount == stl_rowSize );
#endif 

    double *pivotalFront = 
        (double*) paru_calloc (rowCount*fp, sizeof (double), cc);

    if (pivotalFront == NULL )
    {
        printf ("%% Out of memory when tried to allocate for pivotal part %ld",
                f);
        paru_free ( num_panels, sizeof (Int), panel_row, cc);
        return;
    }

    PRLEVEL (1, ("%% pivotalFront =%p \n", pivotalFront));
    Int fm = LUsym->Fm[f];     /* Upper bound number of rows of F */ 
    PRLEVEL (1, ("%% fm=%ld rowCount=%ld \n", fm, rowCount));
    ASSERT ( fm >= rowCount );
    //freeing extra space for rows
    if (rowCount != fm)
    {
        Int sz = sizeof(Int)*fm; 
        frowList =
            (Int*) paru_realloc (rowCount, sizeof(Int), frowList, &sz, cc);
        paruMatInfo ->frowList[f] = frowList;
    }

    paru_fac *LUs =  paruMatInfo->partial_LUs;
    paruMatInfo->frowCount[f] = rowCount;

    LUs[f].m = rowCount;
    LUs[f].n = fp;
    ASSERT (LUs[f].p == NULL);
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


    Int ii = 0; // using for resizing pivotal_elements
    for(Int i = 0 ; i < pivotal_elements.size(); i++)
    {
#ifndef NDEBUG 
        Int p = 1;
#endif
        Int e = pivotal_elements[i];
        paru_Element *el = elementList[e];
        //ASSERT(el != NULL);
        if (el == NULL) continue;
        PRLEVEL (p, ("current element(%ld) %p", e, el ));
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


#ifndef NDEBUG // print the element which is going to be assembled from
        p = 2;
        PRLEVEL (p, ("%% ASSEMBL element= %ld  mEl =%ld ",e, mEl));
        if (p <= 0)
            paru_print_element (paruMatInfo, e);
#endif

        Int cEl = el->lnc;

        PRLEVEL (p, ("%% cEl =%ld \n", cEl));
        for ( ; cEl < nEl ; cEl++)
        {
            if (el_colIndex[cEl] < 0)  // already assembled somewhere
                continue;

            if (el_colIndex[cEl] >= col2 )  
                break;

            PRLEVEL (p, ("ACCESS to %p\n", el_colIndex));

            Int colIndexF =  el_colIndex[cEl] - col1;
            PRLEVEL (p, ("%% Inside the loop cEl =%ld \n", cEl));
            PRLEVEL (p, ("%% colIndexF =%ld \n", colIndexF));

            //double *el_Num = numeric_pointer (el);
            double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);

            //assemble cEl

            assemble_col (el_Num +cEl*mEl,
                    pivotalFront+colIndexF*rowCount,
                    mEl, rowRelIndex);

            el_colIndex[cEl] = flip (el_colIndex[cEl] );
            el->ncolsleft--;     
            if (el->ncolsleft == 0)
            { //free el
                Int tot_size = sizeof(paru_Element) +
                    sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
                PRLEVEL (-1, ("%% Free %ld  %p size %ld\n",e, el, tot_size));
                paru_free (1, tot_size, el, cc);
                elementList[e] = NULL;
                break;
            }
#ifndef NDEBUG  // Printing the pivotal front
            p = 2;
            PRLEVEL (p, ("%% After Assemble element %ld\n", e));
            PRLEVEL (p, ("%% x =  \t"));
            for (Int c = col1; c < col2; c++) 
                PRLEVEL (p, ("%ld\t\t", c));
            PRLEVEL (p, (" ;\n"));
            for (Int r = 0; r < rowCount; r++)
            {
                PRLEVEL (p, ("%% %ld\t", frowList [r]));
                for (Int c = col1; c < col2; c++)
                    PRLEVEL (p, (" %2.5lf\t", 
                                pivotalFront [(c-col1)*rowCount + r]));
                PRLEVEL (p, ("\n"));
            }
            p = 1;
#endif


        }
        if (elementList[e] != NULL )
        {
            el->lnc = cEl;
            pivotal_elements [ii++] = pivotal_elements [i];
            ASSERT (cEl < nEl);
            PRLEVEL (1, ("%%el->lnc= %ld ",el->lnc));
            PRLEVEL (1, ("el_colIndex[el->lnc]=%ld :\n"
                        , el_colIndex[el->lnc]));
#ifndef NDEBUG // print the element which has been assembled from
            p = 2;
            PRLEVEL (p, ("%% ASSEMBLED element= %ld  mEl =%ld ",e, mEl));
            if (p <= 0)
                paru_print_element (paruMatInfo, e);
#endif
        }

    }

    if ( ii+1 < pivotal_elements.size())
    {
        PRLEVEL (1, ("%% pivotal size was %ld ", pivotal_elements.size()));
        PRLEVEL (1, ("%% and now is %ld\n ", ii+1));
        pivotal_elements.resize(ii+1);
    }

#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%% pivotal columns eli(%ld) after resizing: ", eli));
    for(Int i=0 ; i < pivotal_elements.size(); i++)
        PRLEVEL (p, ("%ld ", pivotal_elements[i]));
    PRLEVEL (p, ("\n"));

    p = 2;
    PRLEVEL (p, ("%% After all the assemble\n"));
    PRLEVEL (p, ("%% x =  \t"));
    for (Int c = col1; c < col2; c++) 
        PRLEVEL (p, ("%ld\t\t", c));
    PRLEVEL (p, (" ;\n"));
    for (Int r = 0; r < rowCount; r++)
    {
        PRLEVEL (p, ("%% %ld\t", frowList [r]));
        for (Int c = col1; c < col2; c++)
            PRLEVEL (p, (" %2.5lf\t", pivotalFront [(c-col1)*rowCount + r]));
        PRLEVEL (p, ("\n"));
    }
#endif

    rowMarkp[eli] += rowCount;
    PRLEVEL (1, ("%% rowMarkp[%ld] =%ld\n", eli, rowMarkp[eli]));
    return;
}
