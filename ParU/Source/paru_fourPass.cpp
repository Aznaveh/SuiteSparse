/** =========================================================================  /
 * =======================  paru_fourPass ===================================  /
 * ========================================================================== */
/*! This function is the same level as paru_assemble 
 *      basically doing the rest of the work for non pivotal rows and columns
 */

#include "Parallel_LU.hpp"

void paru_fourPass (paru_matrix *paruMatInfo,
        Int f,        
        Int fp,             //fp of current front
        cholmod_common *cc)
{

    DEBUGLEVEL(1);
    work_struct *Work =  paruMatInfo->Work;

    Int *isRowInFront = Work->rowSize; 
    Int rowMark = Work->rowMark;

    // Couning how many rows/cols of an element is seen
    Int *elRow = Work -> elRow; 
    Int elRMark = Work -> elRMark;
    Int *elCol = Work -> elCol;
    Int elCMark = Work -> elCMark;

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int el_ind = snM [f]; 


    Element **elementList = paruMatInfo->elementList;
    Element *curCB = elementList[el_ind]; 
    Int rowCount= curCB->nrows + fp;
    Int colCount = curCB->ncols;

    Int *CBColList = Work -> scratch + 2*rowCount; //scratch=[fsRowList..ipiv..]
    Int *isColInCBcolSet = Work -> colSize;
    Int colMark = Work -> colMark;


    Int time_f = ++paruMatInfo->time_stamp[f]; //invalid all the markings




    /*****  1st path: over non pivotal columns to count rows             ******/

    tupleList *ColList = paruMatInfo->ColList;
    for (Int k = 0; k < colCount; k++){
        Int c = CBColList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG        
        Int p = 0;
        PRLEVEL (p, ("\n %%--------> 1st: c =%ld  numTuple = %ld\n", 
                    c, numTuple));
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            ASSERT (e >= 0);
            if (e == el_ind){ //current element
                continue;
            }
#ifndef NDEBUG        
            if (p <= 0)
                paru_print_element (paruMatInfo, e);
#endif
            Int curColIndex = curTpl.f;
            Element *curEl = elementList[e];
            Int *el_colIndex = colIndex_pointer (curEl);

            Int *rowRelIndex = relRowInd (curEl);

            if (el_colIndex [curColIndex] < 0 ){/*!TODO: Dead Delete the tuple*/
                continue;  
            }

            if(curEl->cValid !=  time_f){
                curEl->cValid =  time_f;
                elCol [e] = curEl->ncolsleft - 1; //initiaze
                PRLEVEL (1, ("%%cValid=%ld \n",curEl->cValid));
                PRLEVEL (1, ("%%first time seen elCol[e]=%ld \n", elCol[e]));
            }
            else{ 
                elCol [e]--; 
                PRLEVEL (1, ("%%seen before: elCol[e]=%ld \n", elCol[e]));
                continue;
            }
        }
    }
    /**************************************************************************/

    /*****  2nd path: over rows to count columns                         ******/
    Int *fsRowList = Work->scratch; // fully summed row list
    tupleList *RowList = paruMatInfo->RowList;
    for (Int k = fp; k < rowCount; k++){
        Int r = fsRowList [k];
        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG        
        Int p = 0;
        PRLEVEL (p, ("\n %%--------> 2nd r =%ld  numTuple = %ld\n"
                    , r, numTuple));
        if (p <= 0)
            paru_print_tupleList (RowList, r);
#endif
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;

            if (e == el_ind){ //current element
                continue;
            }
#ifndef NDEBUG        
            if (p <= 0)
                paru_print_element (paruMatInfo, e);
#endif
            Int curRowIndex = curTpl.f;
            Element *curEl = elementList[e];
            Int *el_colIndex = colIndex_pointer (curEl); //pointers to row index
            Int *el_rowIndex = rowIndex_pointer (curEl);

            if (el_rowIndex [curRowIndex] < 0 ){  /*! TODO: Dead Delete it	 */
                continue;  
            }

            if(curEl->rValid != time_f){
                curEl->rValid =  time_f;
                elRow [e] = curEl ->nrowsleft - 1; //initiaze
                PRLEVEL (1, ("%%rValid=%ld \n",curEl->rValid));
                PRLEVEL (1, ("%%first time seen elRow[e]=%ld \n", elRow[e]));
            }
            else{ 
                elRow [e]--;
                PRLEVEL (1, ("%%seen before: elRow[e]=%ld \n", elRow[e]));
                continue;
            }

        }
    }
    /**************************************************************************/

    Int *cb_colIndex = colIndex_pointer (curCB);
    Int *cb_rowIndex = rowIndex_pointer (curCB);
    Int *cb_rowRelIndex = relRowInd (curCB);
    double *cb_numbers = numeric_pointer (curCB);

    time_f = ++paruMatInfo->time_stamp[f]; //invalid all the markings

    /*****                 3rd path: assemble columns                    ******/
    /* it contains adding and deleting tuples*/
    for (Int k = 0; k < colCount; k++){
        Int c = CBColList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG            
        Int p = 0;
        PRLEVEL (p, ("\n %%-------->  3rd: c =%ld  numTuple = %ld\n",
                    c, numTuple));
        if (p <= 0 ){
            paru_print_tupleList (ColList, c);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif
        //  Int pdst = 0,psrc;

        for (Int i = 0; i < numTuple; i++){
            //for (Int psrc = 0; psrc < numTuple; psrc ++){}
            //            Tuple curTpl = listColTuples [psrc];
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            if (e == el_ind){ //current element
                //       listColTuples [pdst++] = curTpl; //keeping the tuple
                continue;
            }
            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curColIndex));

            ASSERT (e >= 0);
            ASSERT (e != el_ind);
            ASSERT (curColIndex >= 0);

            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;

            Int *el_colIndex = colIndex_pointer (curEl);
            Int *el_rowIndex = rowIndex_pointer (curEl);
            Int *rowRelIndex = relRowInd (curEl);
            Int *colRelIndex    = relColInd (curEl);

            if (el_colIndex [curColIndex] < 0 ){ //it will be deleted here
                continue;  
            }

            ASSERT (el_colIndex[curColIndex] == c);
            ASSERT (curColIndex < nEl);

            double *el_Num = numeric_pointer (curEl);
            PRLEVEL (1, ("%% elRow[%ld]=%ld currVal= %ld elRMark=%ld\n", 
                        e, elRow[e], curEl->rValid , elRMark));
            PRLEVEL (1, ("%% time_f =%ld \n", time_f));

            //if (elRow [e] - elRMark == 0){}
            if (elRow [e] == 0 && curEl->rValid == time_f -1 ){
                //all the columns are in CB
                if(curEl->rValid !=  time_f){ /*  Update rowRelIndex	 */
                    PRLEVEL (1, ("%% update row relative element%ld\n",e ));
#ifndef NDEBUG
                    //Printing the contribution block prior index update 
                    p = 0;
                    if (p <= 0){
                        PRLEVEL (p, ("\n%%Before index update %ld:",e));
                        paru_print_element (paruMatInfo, e);
                    }
#endif
                    paru_update_rel_ind (curCB, curEl, 'r') ;
#ifndef NDEBUG            
                    for(Int i=0; i < curEl->nrows; i++){
                        PRLEVEL (1, ("%% rowRelIndex[%ld] =%ld\t", i,
                                    rowRelIndex [i]));
                        ASSERT(rowRelIndex [i] < curCB->nrows);
                        PRLEVEL (1,("\n"));
                    }
#endif

                    curEl->rValid =  time_f;
                }
#ifndef NDEBUG
                //Printing the contribution block after prior blocks assembly
                p = 0;
                if (p <= 0){
                    PRLEVEL (p, ("\n%%Before column assembly of %ld:",e));
                    paru_print_element (paruMatInfo, el_ind);
                    paru_print_element (paruMatInfo, e);
                }
#endif


                //                PRLEVEL (1, ("%% psrc=%ld\n", psrc));
                PRLEVEL (1, ("%%colCount=%ld k=%ld", colCount, k));
                PRLEVEL (1, ("%%curCB->nrows=%ld ", curCB->nrows));
                PRLEVEL (1, ("%% cb_numbers=%2.4lf\n",
                            *(cb_numbers+k*curCB->nrows)));

                assemble_col (el_Num+curColIndex*mEl,cb_numbers+k*curCB->nrows,
                        mEl, rowRelIndex);
                colRelIndex [curColIndex] = -1;
                el_colIndex [curColIndex] = -1;
                curEl->ncolsleft --;
#ifndef NDEBUG
                //Printing the contribution block after prior blocks assembly
                p = 0;
                PRLEVEL (p, ("\n%%After column assembly of %ld:",e));
                if (p <= 0){
                    paru_print_element (paruMatInfo, el_ind);
                    paru_print_element (paruMatInfo, e);
                }
#endif


                //  paru_remove_colTuple (ColList, c, i);
                //  i--; numTuple--;
            } 
            //            else 
            //      listColTuples [pdst++] = curTpl; //keeping the tuple
        }

        //        curColTupleList->numTuple = pdst;
        //        PRLEVEL (1, ("%% pdst=%ld\n", pdst));
        paru_print_tupleList (ColList, c);

    }


    /**********************************************************************/

    /****               4th path: assemble rows                        ****/
    for (Int k = fp; k < rowCount; k++){
        Int r = fsRowList [k];
        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG            
        Int p = 0;
        PRLEVEL (1, ("\n %%-------->  4th: r =%ld  numTuple = %ld\n",
                    r, numTuple));
        if (p <= 0 ){
            paru_print_tupleList (RowList, r);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif

        //       Int pdst = 0,psrc;
        for (Int i = 0; i < numTuple; i++){
            //      for (Int psrc = 0; psrc < numTuple; psrc ++){}
            //     Tuple curTpl = listRowTuples [psrc];
            Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;
            if (e == el_ind){ //current element
                //  listRowTuples [pdst++] = curTpl; //keeping the tuple
                continue;
            }
            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curColIndex));


            Int curRowIndex = curTpl.f;
            ASSERT (e >= 0);
            ASSERT (curRowIndex >= 0);

            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;

            Int *el_colIndex = colIndex_pointer (curEl);
            Int *el_rowIndex = rowIndex_pointer (curEl);
            Int *rowRelIndex = relRowInd (curEl);
            Int *colRelIndex    = relColInd (curEl);


            if (el_rowIndex [curRowIndex] < 0 ){ // it will be deleted here
                continue;  
            }
            if (el_rowIndex[curRowIndex] < 0)     continue; //not valid
            PRLEVEL (1, ("%% el_rowIndex [%ld] =%ld\n", 
                        curRowIndex, el_rowIndex [curRowIndex]));
            ASSERT (el_rowIndex[curRowIndex] == r);
            ASSERT (curRowIndex < mEl);

            double *el_Num = numeric_pointer (curEl);
            PRLEVEL (1, ("%% elCol[%ld]=%ld elCMark=%ld\n", 
                        e, elCol[e], elCMark));

            //if (elCol [e] - elCMark == 0){}
            if (elCol [e] == 0 && curEl->cValid == time_f -1){

#ifndef NDEBUG            
                Int p = 0;
                PRLEVEL (1, ("%% Before row assembly: \n" ));
                if (p <= 0 ){
                    paru_print_element (paruMatInfo, e);
                    paru_print_element (paruMatInfo, el_ind);
                }
#endif
                //all the row is in CB
                //
                if(curEl->cValid !=  time_f){
                    /* Update colRelIndex	 */
                    PRLEVEL (1, ("%% update column relative index %ld\n"
                                ,e ));
                    paru_update_rel_ind (curCB, curEl, 'c') ;
#ifndef NDEBUG            
                    for(Int i=0 ; i <curEl->ncols ; i++){
                        PRLEVEL (1, ("%% colRelIndex[%ld] =%ld\t", i,
                                    colRelIndex [i]));
                        ASSERT(colRelIndex [i] < curCB->ncols);
                        PRLEVEL (1,("\n"));
                    }
#endif
                    curEl->cValid =  time_f;
                }
                assemble_row (el_Num, cb_numbers, mEl, nEl, 
                        curCB->nrows, curRowIndex , k-fp, colRelIndex );
#ifndef NDEBUG            
                PRLEVEL (1, ("%% after row assembly: \n" ));
                if (p <= 0 ){
                    paru_print_element (paruMatInfo, e);
                    paru_print_element (paruMatInfo, el_ind);
                }
#endif
                //                    PRLEVEL (1, ("%% psrc=%ld", psrc));
                rowRelIndex [curRowIndex] = -1;
                el_rowIndex [curRowIndex] = -1;
                curEl->nrowsleft --;
            } 
            //                else 
            //       listRowTuples [pdst++] = curTpl; //keeping the tuple
        }

        //            PRLEVEL (1, ("%%   pdst=%ld", pdst));
        //            curRowTupleList->numTuple = pdst;
    }


    /********************* 5th path: clearing column tuples and unhceck *******/
    for (Int k = 0; k < colCount; k++){
        Int c = CBColList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG            
        Int p = 1;
        PRLEVEL (p, ("\n %%-------->  5th: c =%ld  numTuple = %ld\n",
                    c, numTuple));
        if (p <= 0 ){
            paru_print_tupleList (ColList, c);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif
        Int pdst = 0,psrc;

        for (Int psrc = 0; psrc < numTuple; psrc ++){
            Tuple curTpl = listColTuples [psrc];
            Int e = curTpl.e;
            if (e == el_ind){ //current element
                listColTuples [pdst++] = curTpl; //keeping the tuple
                continue;
            }

            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curColIndex));



            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;

            Int *el_colIndex = colIndex_pointer (curEl);
            Int *el_rowIndex = rowIndex_pointer (curEl);
            Int *rowRelIndex = relRowInd (curEl);
            Int *colRelIndex    = relColInd (curEl);

            elCol[e]= -1;

            if (el_colIndex [curColIndex] < 0 ){ //it will be deleted here
                PRLEVEL (1, ("%% psrc=%ld\n", psrc));
                continue;  
            }
            else 
                listColTuples [pdst++] = curTpl; //keeping the tuple

            PRLEVEL (1, ("%% time_f =%ld \n", time_f));

        }
        curColTupleList->numTuple = pdst;
        PRLEVEL (1, ("%% new number of tuple=%ld\n", pdst));
        paru_print_tupleList (ColList, c);

    }
    /**************************************************************************/


    /********************* 6th path: clearing row tuples and unhceck **********/
    for (Int k = fp; k < rowCount; k++){
        Int r = fsRowList [k];
        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG            
        Int p = 1;
        PRLEVEL (p, ("\n %%-------->  6th: r =%ld  numTuple = %ld\n",
                    r, numTuple));
        if (p <= 0 ){
            paru_print_tupleList (RowList, r);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif
        Int pdst = 0,psrc;

        for (Int psrc = 0; psrc < numTuple; psrc ++){
            Tuple curTpl = listRowTuples [psrc];
            Int e = curTpl.e;
            if (e == el_ind){ //current element
                listRowTuples [pdst++] = curTpl; //keeping the tuple
                continue;
            }

            Int curRowIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curRowIndex));


            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;

            Int *el_rowIndex = rowIndex_pointer (curEl);
            elRow[e]= -1;

            if (el_rowIndex [curRowIndex] < 0 ){ //it will be deleted here
                PRLEVEL (1, ("%% psrc=%ld\n", psrc));
                continue;  
            }
            else 
                listRowTuples [pdst++] = curTpl; //keeping the tuple
            PRLEVEL (1, ("%% time_f =%ld \n", time_f));

        }
        curRowTupleList->numTuple = pdst;
        PRLEVEL (1, ("%% new number of tuples=%ld\n", pdst));
        paru_print_tupleList (RowList, r);

    }


}
