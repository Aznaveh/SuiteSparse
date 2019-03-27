/** =========================================================================  /
 * =======================  paru_update_npColLst ============================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

#ifndef NDEBUG  // using STL for debugging
#include <iostream>
#include <algorithm>
#include <set>
#endif

/*! @brief  growing current front if necessary
 * 
 *  @author Aznaveh
 *
 * @param  
 * @return 0 on sucess 
 */

void paru_update_npColLst (paru_matrix *paruMatInfo, Int panel_num ){

    work_struct *Work =  paruMatInfo->Work;

    /**** 4 ******** finding set of non pivotal cols in current front *********/

    /*               
     *
     *
     *                <----------fp--------->
     *                        j1     j2
     *                         ^     ^
     *                         |     | CBColList  Update here
     *                         |     |        \
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
    m = paruMatInfo-> m;
    n = paruMatInfo-> n;

    if (colMark < 0) {  // in rare case of overflow
        //        memset (isRowInFront, -1, n*sizeof(Int)); //Bug detected
        memset (isColInCBcolSet , -1, n*sizeof(Int));
        colMark = Work-> colMark = 0;
    }
    Int *CBColList = Work -> scratch + 2*rowCount;//scratch=[fsRowList..ipiv..]
    Int *fsRowList = Work->scratch; // fully summed row list
    Int colCount = 0;

#ifndef NDEBUG
    std::set<Int> stl_colSet;
#endif  

    tupleList *RowList = paruMatInfo->RowList;
    for (Int i = 0; i < fp; i++){
        Int curFsRowIndex =(Int) i; //current fully summed row index
        Int curFsRow = fsRowList [i];
        PRLEVEL (1, ("%% 4: curFsRowIndex = %ld\n", curFsRowIndex));
        PRLEVEL (1, ("%% curFsRow =%ld\n", curFsRow));
        tupleList *curRowTupleList = &RowList [curFsRow];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listRowTuples = curRowTupleList->list;
        PRLEVEL (1, ("%% 4: numTuple = %ld\n", numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;
            Int curRowIndex = curTpl.f;
            if(e < 0 || curRowIndex < 0) continue;

            Element *el = elementList[e];
            Int mEl = el->nrows;
            Int nEl = el->ncols;
            Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_colIndex = colIndex_pointer (el);
            Int *colRelIndex = relColInd (el);
            Int *rowRelIndex = relRowInd (el);

            if (el_rowIndex [curRowIndex] < 0 ) continue;
            ASSERT (el_rowIndex[curRowIndex] == curFsRow);

            rowRelIndex [curTpl.f] = curFsRow;

            if(el->cValid !=  time_f){// an element never seen before
                el->cValid = time_f;
#ifndef NDEBUG            
                if (el->cValid >  time_f )
                    PRLEVEL (0, ("%%time_f =%ld  cVal= %ld\n", 
                                time_f , el->cValid));
#endif    
                ASSERT(el->cValid <= time_f);
#ifndef NDEBUG
                if ( elCol [e] >= elCMark )
                    PRLEVEL (1, ("%% element %ld can be eaten wholly\n",e));
                //And the rest of e is in U part 
#endif
            }
            else { // must not happen anyway; it depends on changing strategy
                elRow [e]--;
                continue;
            }

            PRLEVEL (1, ("%% element= %ld  nEl =%ld \n",e, nEl));
            for (Int cEl = 0; cEl < nEl; cEl++){
                Int curCol = el_colIndex [cEl]; 
                PRLEVEL (1, ("%% curCol =%ld\n", curCol));
                ASSERT (curCol < n);
                if (curCol < 0)
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
                    CBColList [colCount] = curCol;
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
        PRLEVEL (1, ("%% pivotalFront =%p\n",pivotalFront));
        return;
    }

#ifndef NDEBUG /* Checking if columns are correct */
    p = 1;
    PRLEVEL (p, ("%% There are %ld columns in this contribution block: \n",
                colCount));
    for (Int i = 0; i < colCount; i++)
        PRLEVEL (p, ("%%  %ld", CBColList [i]));
    PRLEVEL (p, ("\n"));
    Int stl_colSize = stl_colSet.size();
    if (colCount != stl_colSize){
        PRLEVEL (p, ("%% STL %ld:\n",stl_colSize));
        for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
            PRLEVEL (p, ("%%  %ld", *it));
        PRLEVEL (p, ("\n%% My Set %ld:\n",colCount));
        for (Int i = 0; i < colCount; i++)
            PRLEVEL (p, ("%%  %ld", CBColList [i]));
        PRLEVEL (p, ("\n"));
    }
    ASSERT (colCount == stl_colSize );
#endif 
    /************* travers over new non pivotal columns ***********************/

    /*               
     *
     *
     *                <----------fp--------->
     *                                                  
     *                                 CBColList      ^         ^         
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



    /********* travers over new non pivotal rows of current panel *************/


    /*               
     *
     *
     *                <----------fp--------->
     *                        j1     j2
     *                         ^     ^
     *                         |     | CBColList  Update here
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
     * |          E   |       |     |       |       00000000
     * |              |       |     |       |          vvvvvvvvv
     * |          R   |       |     |       |          xxxxxxxxxxxxxxxxxxxx
     * |          E   |       |row_end      |          xxxxxx El xxxxxxxxxx
     * | row_end----->.       .      .      .          xxxxxxxxxxxxxxxxxxxx
     * |              .       .      .      . 
     * |              .       .      .      . 
     * v              |___....______________|              
     *                         
     */


}
