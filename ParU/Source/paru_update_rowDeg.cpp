/** =========================================================================  /
 * =======================  paru_update_rowDeg   ============================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

/*! @brief  growing current front if necessary and update the row degree of
 *   current front for current panel. 
 * 
 *  @author Aznaveh
 *
 *  @param  
 */

void paru_update_rowDeg ( Int panel_num,  Int row_end, Int f, Int *next,
        std::set<Int> &stl_colSet, paru_matrix *paruMatInfo)
{

    DEBUGLEVEL(0);
#ifndef NDEBUG
    Int p = 1;
    static Int r1 = 0, r2 = 0, r3 = 0 ;
#endif
    PRLEVEL (1, ("%%-------ROW degree update of panel %ld of front %ld \n", 
             panel_num, f ));
    Int panel_width = paruMatInfo->panel_width;
    paru_Element **elementList = paruMatInfo->elementList;
    work_struct *Work =  paruMatInfo->Work;

    Int *elRow = Work -> elRow; 
    Int *elCol = Work -> elCol;

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *Super = LUsym->Super;
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];
    Int fp = col2 - col1;   /* first fp columns are pivotal */ 

    Int time_f = ++paruMatInfo->time_stamp[f]; //making all the markings invalid
    Int npMark = time_f; //Mark for non pivotal rows
    Int pMark = npMark; pMark++;    //Mark for pivotal rows

    Int colCount = stl_colSet.size();

    Int j1 = panel_num*panel_width; // panel starting column
    Int j2 = (j1 + panel_width < fp ) ? 
        j1 + panel_width : fp;


    Int rowCount = paruMatInfo->frowCount[f];
    Int *row_degree_bound = paruMatInfo->row_degree_bound;


    std::set<Int> stl_newColSet;  // the list of new columns

    /*************** finding set of non pivotal cols in current front *********/

    /*               
     *    Mark seen elements with pMark or time_f+1
     *        if Marked already added to the list
     *        else all columns are added to the current front
     *
     *                <----------fp--------->           
     *                        j1     j2              Update here
     *                         ^     ^                stl_newColSet colCount
     *                         |     | fcolList       ^ . . .   ^
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



    Int n = paruMatInfo-> n;


    Int *frowList = paruMatInfo->frowList[f];


    std::set<Int>::iterator it;

    tupleList *RowList = paruMatInfo->RowList;
    for (Int i = j1; i < j2; i++)
    {
#ifndef NDEBUG
        Int curFsRowIndex = i; //current fully summed row index
#endif  
        Int curFsRow = frowList [i];
        PRLEVEL (1, ("%% 4: curFsRowIndex = %ld\n", curFsRowIndex));
        PRLEVEL (1, ("%% curFsRow =%ld\n", curFsRow));
        tupleList *curRowTupleList = &RowList [curFsRow];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listRowTuples = curRowTupleList->list;
        PRLEVEL (1, ("%% 4: numTuple = %ld\n", numTuple));


        Int pdst = 0, psrc;
        for (psrc = 0; psrc < numTuple; psrc ++)
        {
        //for (Int i = 0; i < numTuple; i++)
        //{
           // paru_Tuple curTpl = listRowTuples [i];
            paru_Tuple curTpl = listRowTuples [psrc];

            Int e = curTpl.e;
            Int curRowIndex = curTpl.f;
#ifndef NDEBUG
            r1++;
#endif
            if(e < 0 || curRowIndex < 0) continue;

            paru_Element *el = elementList[e];

            if (el == NULL) continue;

            Int nEl = el->ncols;
            //Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_rowIndex = (Int*)(el+1)+nEl;
            if (el_rowIndex [curRowIndex] < 0 ) continue;

            Int mEl = el->nrows;
            //Int *rowRelIndex = relRowInd (el);
            Int *rowRelIndex = (Int*)(el+1) + 2*nEl + mEl;
            rowRelIndex [curTpl.f] = curFsRow;

            listRowTuples [pdst++] = curTpl; //keeping the tuple

            ASSERT (el_rowIndex[curRowIndex] == curFsRow);


            if(el->rValid !=  pMark)
            {// an element never seen before
                el->rValid = pMark;

                PRLEVEL (1, ("%%first time seen elRow[%ld]=%ld \n"
                            , e, elRow[e]));
                if (el->rValid != npMark)
                    elRow [e] = el ->nrowsleft - 1; //initiaze
                else
                    elRow [e]--; //already seen in a non pivotal row

                PRLEVEL (1, ("%%changed to elRow[%ld]=%ld \n", e, elRow[e]));
#ifndef NDEBUG            
                if (el->rValid >  pMark)
                    PRLEVEL (1, ("%%pMark=%ld  rVal= %ld\n", 
                                pMark, el->rValid));
                ASSERT(el->rValid <= pMark);
                if ( elCol [e] >= Work -> elCMark )
                    PRLEVEL (1, ("%% element %ld can be eaten wholly\n",e));
                //And the rest of e is in U part 
#endif
            }
            else 
            {  //already added to pivotal rows
                elRow [e]--;
                PRLEVEL (1, ("%%already seen elRow[%ld]=%ld \n",e, elRow[e]));
                continue;
            }

            //Int *el_colIndex = colIndex_pointer (el);
            Int *el_colIndex = (Int*)(el+1);

            //Int *colRelIndex = relColInd (el);
            Int *colRelIndex = (Int*)(el+1) + mEl + nEl;


            PRLEVEL (1, ("%% element= %ld  nEl =%ld \n",e, nEl));
            for (Int cEl = 0; cEl < nEl; cEl++)
            {
                Int curCol = el_colIndex [cEl]; 
                PRLEVEL (1, ("%% curCol =%ld\n", curCol));
                ASSERT (curCol < n);

                if (curCol < 0 )  //already deleted
                    continue;
                
                if (curCol < col2 && curCol >= col1 )  /*is a pivotal col */ 
                    continue;

                if ( stl_colSet.find(curCol) == stl_colSet.end() )
                {
                    stl_colSet.insert (curCol);
                    stl_newColSet.insert (curCol);
                    colCount++; 
                }
#ifndef NDEBUG
                p=1;
                //stl_colSet.insert (curCol);
                for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
                    PRLEVEL (p, ("%%@  %ld", *it));

#endif

                ASSERT (colCount <= n);
            }
        }

        curRowTupleList->numTuple = pdst;
        }

        if (colCount == 0)
        {  // there is no CB, Nothing to be done
            Work->rowMark +=  rowCount;
            return;
        }

#ifndef NDEBUG /* Checking if columns are correct */
        p = 1;
        PRLEVEL (p, ("%% There are %ld columns in this contribution block: \n",
                    colCount));
        PRLEVEL (p, ("\n"));
        Int stl_colSize = stl_colSet.size();

        if (colCount != stl_colSize)
        {
            PRLEVEL (p, ("%% STL %ld:\n",stl_colSize));
            for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
                PRLEVEL (p, ("%%  %ld", *it));
            PRLEVEL (p, ("\n%% My Set %ld:\n",colCount));
            PRLEVEL (p, ("\n"));
        }
        ASSERT (colCount == stl_colSize );
#endif 

    // if the front did not grow, there is nothing else to do
    if (stl_newColSet.size() == 0) 
        return;

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

    for (it = stl_newColSet.begin(); it != stl_newColSet.end(); it++)
    {
        //TODO: take c with another strategy
        Int c = *it;

        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG        
        p = 1;
        PRLEVEL (p, ("\n %%--------> 1st: c =%ld  numTuple = %ld\n", 
                    c, numTuple));
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif
        for (Int i = 0; i < numTuple; i++)
        {
            paru_Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
#ifndef NDEBUG
            r2++;
#endif

            if ( e >= el_ind || e < first[el_ind])
                continue;
#ifndef NDEBUG        
            if (p <= 0)
                paru_print_element (paruMatInfo, e);
#endif
            paru_Element *el = elementList[e];
            if (el == NULL) continue;

            Int curColIndex = curTpl.f;
            //Int *el_colIndex = colIndex_pointer (el);
            Int *el_colIndex = (Int*)(el+1);

            if (el_colIndex [curColIndex] < 0 )
                continue;  

            if(el->cValid !=  time_f)
            {
                el->cValid =  time_f;
                elCol [e] = el->ncolsleft - 1; //initiaze
                PRLEVEL (1, ("%%cValid=%ld \n",el->cValid));
                PRLEVEL (1, ("%%first time seen elCol[e]=%ld \n", elCol[e]));
            }
            else
            { 
                elCol [e]--; 
                PRLEVEL (1, ("%%seen before: elCol[e]=%ld \n", elCol[e]));
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
 
    for (Int k = j2; k < row_end; k++)
    {
        Int r = frowList [k];

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
        Int pdst = 0, psrc;
        for (psrc = 0; psrc < numTuple; psrc ++)
        {

            paru_Tuple curTpl = listRowTuples [psrc];
            Int e = curTpl.e;
#ifndef NDEBUG
            r3++;
#endif


#ifndef NDEBUG        
            if (p <= 0)
                paru_print_element (paruMatInfo, e);
#endif
            Int curRowIndex = curTpl.f;

            if(e < 0 || curRowIndex < 0) continue;

            paru_Element *el = elementList[e];
            if (el == NULL) continue;

            Int *el_rowIndex = rowIndex_pointer (el);

            if (el_rowIndex [curRowIndex] < 0 )
            {
                continue;  
            }
            listRowTuples [pdst++] = curTpl; //keeping the tuple

            if(el->cValid !=  time_f)
            {
                el->cValid =  time_f;
                elCol [e] = el->ncolsleft ; //initiaze
            }
            
            //TODO: change only if any thing remain
            // if (elRow [e] == 0 && elCol [e] == 0)
            if (elRow [e] == 0)
                new_row_degree_bound_for_r += elCol [e] ;
            else
                new_row_degree_bound_for_r += el->ncolsleft ;


            PRLEVEL (1, ("%% pMark=%ld npMark=%ld \n",pMark, npMark));

            if(el->rValid == pMark)
            {  //already a pivot
                //ASSERT (elRow[e] == 0); //TODO: something is wrong here
                elRow [e]--;
                PRLEVEL (1, ("%% Pivotal elRow[%ld]=%ld \n",e, elRow[e]));
                continue;
            }

            if(el->rValid != npMark)
            {
                el->rValid =  npMark;
                elRow [e] = el ->nrowsleft - 1; //initiaze
                PRLEVEL (1, ("%%rValid=%ld \n",el->rValid));
                PRLEVEL (1, ("%%first time seen elRow[%ld]=%ld \n",
                            e, elRow[e]));
            }
            else
            { 
                elRow [e]--;

                PRLEVEL (1, ("%%seen before: elRow[e]=%ld \n", elRow[e]));
            }
#if 0

            if (elRow [e] == 0) 
            {

                PRLEVEL (1, ("%% elRow[%ld]=%ld \n",e,  elRow[e]));
                if(elCol[e] == 0)
                {
                    // adding the element to the list of children
                    el->next = *next;
                    *next = e;
                }
            }
#endif


        }

        curRowTupleList->numTuple = pdst;
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
        PRLEVEL (1, ("%% Finalized counters r1=%ld r2=%ld r3=%ld sum=%ld\n", 
                    r1, r2, r3, r1+r2+r3));

}
