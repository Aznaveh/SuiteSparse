/** =========================================================================  /
 * =======================  paru_update_rowDeg   ============================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

#ifndef NDEBUG  // using STL for debugging
#include <iostream>
#include <algorithm>
#include <set>
#endif

/*! @brief  growing current front if necessary and update the row degree of
 *   current front for current panel. 
 * 
 *  @author Aznaveh
 *
 *  @param  
 */

void paru_update_rowDeg ( Int panel_num,  Int row_end, 
        Int f, paru_matrix *paruMatInfo){

    DEBUGLEVEL(0);
    PRLEVEL (1, ("%%-------ROW degree update of panel %ld of front %ld \n", 
             panel_num, f ));
    Int panel_width = paruMatInfo->panel_width;
    paru_Element **elementList = paruMatInfo->elementList;
    work_struct *Work =  paruMatInfo->Work;
    Int elCMark = Work -> elCMark;

    Int *elRow = Work -> elRow; 
    Int *elCol = Work -> elCol;

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *Super = LUsym->Super;
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];
    Int fp = col2 - col1;   /* first fp columns are pivotal */ 

    Int time_f = ++paruMatInfo->time_stamp[f]; //making all the markings invalid
    Int npMark = time_f; //Mark for non pivotal rows
    Int pMark = npMark; pMark++;    //Mark for non pivotal rows

    Int colCount = paruMatInfo->fcolCount[f] ;
    Int past_col = colCount;   //saving how many colums are in this front so far

    Int j1 = panel_num*panel_width; // panel starting column
    Int j2 = (j1 + panel_width < fp ) ? 
        j1 + panel_width : fp;


    Int rowCount = paruMatInfo->frowCount[f];
    Int *row_degree_bound = paruMatInfo->row_degree_bound;

    /*************** finding set of non pivotal cols in current front *********/

    /*               
     *    Mark seen elements with pMark or time_f+1
     *        if Marked already added to the list
     *        else all columns are added to the current front
     *
     *                <----------fp--------->
     *                        j1     j2              Update here
     *                         ^     ^             past_col    colCount
     *                         |     | fcolList      ^ . . .   ^
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



    Int *isColInCBcolSet = Work -> colSize;
    Int colMark = Work -> colMark;
    Int m = paruMatInfo-> m;
    Int n = paruMatInfo-> n;

    if (colMark < 0) {  // in rare case of overflow
        //        memset (isRowInFront, -1, n*sizeof(Int)); //Bug detected
        memset (isColInCBcolSet , -1, n*sizeof(Int));
        colMark = Work-> colMark = 0;
    }

    Int *frowList = paruMatInfo->frowList[f];
    Int *fcolList = paruMatInfo->fcolList[f];

#ifndef NDEBUG
    std::set<Int> stl_colSet;
    std::set<Int>::iterator it;
#endif  

    tupleList *RowList = paruMatInfo->RowList;
    for (Int i = j1; i < j2; i++){
        Int curFsRowIndex =(Int) i; //current fully summed row index
        Int curFsRow = frowList [i];
        PRLEVEL (1, ("%% 4: curFsRowIndex = %ld\n", curFsRowIndex));
        PRLEVEL (1, ("%% curFsRow =%ld\n", curFsRow));
        tupleList *curRowTupleList = &RowList [curFsRow];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listRowTuples = curRowTupleList->list;
        PRLEVEL (1, ("%% 4: numTuple = %ld\n", numTuple));
        for (Int i = 0; i < numTuple; i++){
            paru_Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;
            Int curRowIndex = curTpl.f;
            if(e < 0 || curRowIndex < 0) continue;

            paru_Element *el = elementList[e];
            Int mEl = el->nrows;
            Int nEl = el->ncols;
            Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_colIndex = colIndex_pointer (el);
            Int *colRelIndex = relColInd (el);
            Int *rowRelIndex = relRowInd (el);

            if (el_rowIndex [curRowIndex] < 0 ) continue;
            ASSERT (el_rowIndex[curRowIndex] == curFsRow);

            rowRelIndex [curTpl.f] = curFsRow;

            if(el->cValid !=  pMark){// an element never seen before
                el->cValid = pMark;
#ifndef NDEBUG            
                if (el->cValid >  pMark)
                    PRLEVEL (1, ("%%pMark=%ld  cVal= %ld\n", 
                                pMark, el->cValid));
#endif    
                ASSERT(el->cValid <= pMark);
#ifndef NDEBUG
                if ( elCol [e] >= elCMark )
                    PRLEVEL (1, ("%% element %ld can be eaten wholly\n",e));
                //And the rest of e is in U part 
#endif
            }
            else {  //already added to pivotal rows
                //   elRow [e]--;
                continue;
            }

            PRLEVEL (1, ("%% element= %ld  nEl =%ld \n",e, nEl));
            for (Int cEl = 0; cEl < nEl; cEl++){
                Int curCol = el_colIndex [cEl]; 
                PRLEVEL (1, ("%% curCol =%ld\n", curCol));
                ASSERT (curCol < n);

                if (curCol < 0 )  //already deleted
                    continue;
                
                if (curCol < col2 && curCol >= col1 )  /*is a pivotal col */ 
                    continue;
#ifndef NDEBUG
                stl_colSet.insert (curCol);
#endif
                PRLEVEL (1, ("%% %p ---> isColInCBcolSet[%ld]=%ld\n", 
                            isColInCBcolSet+curCol, curCol,
                            isColInCBcolSet[curCol]));

                if (isColInCBcolSet [curCol] < colMark  ){
                    PRLEVEL (1, ("%% curCol = %ld colCount=%ld\n", 
                                curCol, colCount));
                    fcolList [colCount] = curCol;
                    colRelIndex [cEl] = colCount;
                    isColInCBcolSet [curCol] = colMark + colCount++; 
                }
                else{
                    colRelIndex [cEl] = isColInCBcolSet [curCol]- colMark;
                }
                ASSERT (colCount <= n);
            }
        }
    }

    if (colCount == 0){  // there is no CB, Nothing to be done
        Work->rowMark +=  rowCount;
        return;
    }

#ifndef NDEBUG /* Checking if columns are correct */
    Int p = 1;
    PRLEVEL (p, ("%% There are %ld columns in this contribution block: \n",
                colCount));
    for (Int i = 0; i < colCount; i++)
        PRLEVEL (p, ("%%  %ld", fcolList [i]));
    PRLEVEL (p, ("\n"));
    Int stl_colSize = stl_colSet.size();
    if (colCount != stl_colSize){
        PRLEVEL (p, ("%% STL %ld:\n",stl_colSize));
        for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
            PRLEVEL (p, ("%%  %ld", *it));
        PRLEVEL (p, ("\n%% My Set %ld:\n",colCount));
        for (Int i = 0; i < colCount; i++)
            PRLEVEL (p, ("%%  %ld", fcolList [i]));
        PRLEVEL (p, ("\n"));
    }
    //TODO: check this assertion
    //ASSERT (colCount == stl_colSize );
#endif 

    // if the front did not grow, there is nothing else to do
    if (colCount == past_col ){
        return; 
    }

    paruMatInfo->fcolCount[f] = colCount;


    /************* travers over new non pivotal columns ***********************/
    /*               
     *  Marking seen element with pMark or time_f; it would be fine with either
     *  while there is no other column pass
     *   
     *                <----------fp--------->
     *                                                  
     *                                 fcolList       ^         ^         
     *                                        \       |   HERE  |
     *             F                           [QQQQQ|OOOOOOOOO|....
     *              \  ____..._________...__ _____________________________...
     * ^              |\      |     |       #  ^     |         | 
     * |              | \     |     |       #  | old | added   | 
     * |              |  \    |     |       #  | list|  columns|
     * |              |___\...|_____|__...__#  |     |         |
     * |   ^          |       |\* **|       #  fp  oooooooooooo|
     * |   |          |       |**\**|       #  |   ooo El ooooo|
     * | panel        |       |***\*|       #  |   oooooooooooo|
     * | width        |       |***\*|       #  |               |
     * |   |          |       |***\*|       #  |    00000000   |
     * |   v          |____...|________..._ #  |    000 El 0   |
     * |              |       |     |       #  v    00000000   |           ...
     * rowCount       |==================================================
     * |              |       |     |       |     ooooooooooooooooo
     * |              |       |     |       |     ooooooooooooooooo
     * |              |       |row_end      |     oooo EL ooooooooo
     * |              .       .      .      .     ooooooooooooooooo
     * |              .       .      .      . 
     * |              .       .      .      .      xxxxxxxxxx
     * v              |___....______________|      xxx EL xxx        
     *                                             xxxxxxxxxx
     *                                             xxxxxxxxxx
     *                                             
     *                                         oooooooooooo    
     *                                         ooo El ooooo
     *                                         oooooooooooo    
     *                                             
     */
    tupleList *ColList = paruMatInfo->ColList;
    Int *snM = LUsym->super2atree;
    Int el_ind = snM [f]; 
    Int *first = LUsym->first;
    for (Int k = past_col; k < colCount; k++){
        Int c = fcolList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG        
        Int p = 1;
        PRLEVEL (p, ("\n %%--------> 1st: c =%ld  numTuple = %ld\n", 
                    c, numTuple));
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif
        for (Int i = 0; i < numTuple; i++){
            paru_Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            // if (e == el_ind){ //current element}
            if ( e >= el_ind || e < first[el_ind]){ //Not any of descendents
                continue;
            }
#ifndef NDEBUG        
            if (p <= 0)
                paru_print_element (paruMatInfo, e);
#endif
            Int curColIndex = curTpl.f;
            paru_Element *el = elementList[e];
            Int *el_colIndex = colIndex_pointer (el);

            Int *rowRelIndex = relRowInd (el);

            if (el_colIndex [curColIndex] < 0 ){
                continue;  
            }

            if(el->cValid !=  time_f){
                el->cValid =  time_f;
                elCol [e] = el->ncolsleft - 1; //initiaze
                PRLEVEL (1, ("%%cValid=%ld \n",el->cValid));
                PRLEVEL (1, ("%%first time seen elCol[e]=%ld \n", elCol[e]));
            }
            else{ 
                elCol [e]--; 
                PRLEVEL (1, ("%%seen before: elCol[e]=%ld \n", elCol[e]));
                continue;
            }
        }
    }


    /********* travers over new non pivotal rows of current panel *************/
    /*           Marking seen elements by npMark or time_f
     *           if marked with pMark just ignore
     *           if marked with npMark decrease number of rows by 1 
     *           else initialze number of rows seen
     *
     *                <----------fp--------->
     *                        j1     j2
     *                         ^     ^
     *                         |     | fcolList  Update here
     *                         |     |        \
     *             F           |     |         [QQQQQ|OOOOOOOOO|....
     *              \  ____..._|_  ____...__ _____________________________...
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
     * |        j2--> |       |     |       #  v    00000000   |           ...
     * rowCount   H   |=============================00 El 00=============
     * |          E   |       |     |       |       00000000   Update row degree
     * |              |       |     |       |          vvvvvvvvv
     * |          R   |       |     |       |          xxxxxxxxxxxxxxxxxxxx
     * |          E   |       |row_end      |          xxxxxx El xxxxxxxxxx
     * | row_end----->.       .      .      .          xxxxxxxxxxxxxxxxxxxx
     * |              .       .      .      . 
     * |              .       .      .      . 
     * v              |___....______________|              
     *                         
     */
    Int new_row_degree_bound_for_r;
 
    for (Int k = j2; k < row_end; k++){
        Int r = frowList [k];

//        new_row_degree_bound_for_r = rowCount;
        new_row_degree_bound_for_r = colCount;

        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG        
        Int p = 1;
        PRLEVEL (p, ("\n %%--------> 2nd r =%ld  numTuple = %ld\n"
                    , r, numTuple));
        if (p <= 0)
            paru_print_tupleList (RowList, r);
#endif
        for (Int i = 0; i < numTuple; i++){
            paru_Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;

              if (e == el_ind){ //current element}
            //if ( e >= el_ind || e < first[el_ind]){ //Not any of descendents
                continue;
            }


#ifndef NDEBUG        
            if (p <= 0)
                paru_print_element (paruMatInfo, e);
#endif
            Int curRowIndex = curTpl.f;
            paru_Element *el = elementList[e];
            Int *el_colIndex = colIndex_pointer (el);//pointers to row index
            Int *el_rowIndex = rowIndex_pointer (el);

            if (el_rowIndex [curRowIndex] < 0 ){
                continue;  
            }

            if(el->cValid !=  time_f){
                el->cValid =  time_f;
                elCol [e] = el->ncolsleft ; //initiaze
            }
            new_row_degree_bound_for_r += elCol [e] ;

            if(el->rValid == pMark) continue;  //already a pivot

            if(el->rValid != npMark){
                el->rValid =  npMark;
                elRow [e] = el ->nrowsleft - 1; //initiaze
                PRLEVEL (1, ("%%rValid=%ld \n",el->rValid));
                PRLEVEL (1, ("%%first time seen elRow[e]=%ld \n",
                            elRow[e]));
            }
            else{ 
                elRow [e]--;

                PRLEVEL (1, ("%%seen before: elRow[e]=%ld \n", elRow[e]));
            }

        }

//        Int old_bound_updated = row_degree_bound [r] + rowCount - 1 ;
        Int old_bound_updated = row_degree_bound [r] + colCount - 1 ;

#ifndef NDEBUG
        p = 1;
        PRLEVEL (p, ("%%old_bound_updated =%ld \n",old_bound_updated));
        PRLEVEL (p, ("%%new_row_degree_bound_for_r=%ld \n",
                   new_row_degree_bound_for_r));
        PRLEVEL (p, ("%%row_degroo_bound[%ld]=%ld \n",r, row_degree_bound[r]));
#endif
        row_degree_bound [r] =  // min
            old_bound_updated < new_row_degree_bound_for_r ? 
            old_bound_updated : new_row_degree_bound_for_r;
    }

    paruMatInfo->time_stamp[f]+= 2; //making all the markings invalid again
}
