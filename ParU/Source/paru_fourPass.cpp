/** =========================================================================  /
 * =======================  paru_fourPass ===================================  /
 * ========================================================================== */
/*! This function is the same level as paru_assemble 
 *      basically doing the rest of the work for non pivotal rows and columns
 */

#include "Parallel_LU.hpp"

void paru_fourPass (paru_matrix *paruMatInfo,
        Int el_ind,         // Index of current CB
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

    Element **elementList = paruMatInfo->elementList;


    Element *curCB = elementList[el_ind]; 
    Int rowCount= curCB->nrows + fp;
    Int colCount = curCB->ncols;

    Int *CBColList = Work -> scratch + 2*rowCount; //scratch=[fsRowList..ipiv..]
    Int *isColInCBcolSet = Work -> colSize;
    Int colMark = Work -> colMark;



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
        PRLEVEL (p, ("%% 1st: c =%ld  numTuple = %ld\n", c, numTuple));
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            ASSERT (e >= 0);
            Int curColIndex = curTpl.f;
            Element *curEl = elementList[e];
            Int *el_colIndex = colIndex_pointer (curEl);
            Int *el_rowIndex = rowIndex_pointer (curEl); //pointers to row index
            //Int *rowRelIndValid = rowRelIndVal (curEl);
            Int rowRelIndValid = curEl->rValid;

            Int *rowRelIndex = relRowInd (curEl);
 
            if (el_colIndex [curColIndex] < 0 ){/*!TODO: Dead Delete the tuple*/
                continue;  
            }

            if(rowRelIndValid !=  el_ind){//BUGGY, What if I already have seen
                rowRelIndValid =  el_ind;
                elCol [e] = curEl->ncolsleft - 1; //initiaze
            }
            else{ 
                elCol [e]--; 
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
        tupleList *curColTupleList = &RowList[r];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG        
        Int p = 0;
        PRLEVEL (p, ("%% 2nd r =%ld  numTuple = %ld\n", r, numTuple));
        if (p <= 0)
            paru_print_tupleList (RowList, r);
#endif
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Int curRowIndex = curTpl.f;
            Element *curEl = elementList[e];
            Int *el_colIndex = colIndex_pointer (curEl); //pointers to row index
            //Int *colRelIndValid = rowRelIndVal (curEl);
            Int colRelIndValid = curEl->cValid;
            Int *colRelIndex    = relColInd (curEl);
            Int *el_rowIndex = rowIndex_pointer (curEl);

            if (el_rowIndex [curRowIndex] < 0 ){  /*! TODO: Dead Delete it	 */
                continue;  
            }


            if(colRelIndValid != el_ind ){
                colRelIndValid =  el_ind ;
                elRow [e] = curEl ->nrowsleft - 1; //initiaze
            }
            else{ 
                elRow [e]--;
                continue;
            }

       }
    }
    /**************************************************************************/

    Int *cb_colIndex = colIndex_pointer (curCB);
    Int *cb_rowIndex = rowIndex_pointer (curCB);
    Int *cb_rowRelIndex = relRowInd (curCB);
    double *cb_numbers = numeric_pointer (curCB);


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
        PRLEVEL (p, ("%% 3rd: c =%ld  numTuple = %ld\n", c, numTuple));
        if (p <= 0 ){
            paru_print_element (paruMatInfo, el_ind);
            paru_print_tupleList (ColList, c);
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
            PRLEVEL (1, ("%%   numTuple =%ld\n",   numTuple));
            PRLEVEL (1, ("%% element= %ld  mEl =%ld \n",e, mEl));
            PRLEVEL (1, ("%% f =%ld\n", curColIndex));
            PRLEVEL (1, ("%% CB: %ld x %ld\n", rowCount, colCount));
            if (elRow [e] - elRMark == 0){
                //all the column is in CB

/*             /*  Update rowRelIndex	 */
/*            Int mEl = curEl->nrows;
/*           // Updating row relative indices 
/*            PRLEVEL (1, ("%% elRow[%ld]=%ld", e, elRow [e]));  
/*            //            *rowRelIndValid = f ;//current front
/*            for (Int rEl = 0; rEl < mEl; rEl++){   
/*                rowRelIndex [rEl] = isRowInFront [el_rowIndex [rEl]] 
/*                    - rowMark;
/*                ASSERT (rowRelIndex [rEl] < fp );
/*                PRLEVEL (1, ("%% ^^rowRelIndex[%ld] = %ld\n",
/*                            rEl,  rowRelIndex[rEl]));
/*            }
/*            PRLEVEL (1, ("\n%%e=%ld elCol[e]=%ld \n", e, elCol[e]));
/*                
 */
#ifndef NDEBUG
                //Printing the contribution block after prior blocks assembly
                p = 0;
                if (p <= 0){
                    PRLEVEL (p, ("\n%%Before column assembly of %ld:",e));
                    paru_print_element (paruMatInfo, e);
                }
#endif


                PRLEVEL (1, ("%% psrc=%ld\n", psrc));
                PRLEVEL (1, ("%%colCount=%ld k=%ld", colCount, k));

                assemble_col (el_Num+curColIndex*mEl,cb_numbers+k*colCount,
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
            else 
                listColTuples [pdst++] = curTpl; //keeping the tuple
        }

        curColTupleList->numTuple = pdst;
        PRLEVEL (1, ("%% pdst=%ld\n", pdst));
        paru_print_tupleList (ColList, c);

    }


    /**************************************************************************/

    /****                   4th path: assemble rows                        ****/
    for (Int k = fp; k < rowCount; k++){
        Int r = fsRowList [k];
        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG            
        Int p = 0;
        PRLEVEL (1, ("%% 4th: r =%ld  numTuple = %ld\n", r, numTuple));
        if (p <= 0 ){
            paru_print_element (paruMatInfo, el_ind);
            paru_print_tupleList (RowList, r);
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
            PRLEVEL (1, ("%% r =%ld\n", r));
            PRLEVEL (1, ("%% el_rowIndex [%ld] =%ld\n", 
                        curRowIndex, el_rowIndex [curRowIndex]));
            ASSERT (el_rowIndex[curRowIndex] == r);
            ASSERT (curRowIndex < mEl);

            double *el_Num = numeric_pointer (curEl);

            if (elCol [e] - elCMark == 0){
                //all the row is in CB
                //
            /* Update colRelIndex	 */
/*            Int nEl = curEl->ncols;
/*            // Updating row relative indices 
/*            PRLEVEL (1, ("%% elCol[%ld]=%ld", e, elCol [e]));  
/*
/*            for (Int cEl = 0; cEl < nEl; cEl++)   {
/*                colRelIndex [cEl] = isColInCBcolSet [el_colIndex [cEl]] 
/*                    - colMark;
/*                ASSERT (colRelIndex [rEl] < fp );
/*            }
/*            //            *colRelIndValid = f ;//current front
/*
/*            PRLEVEL (1, ("\n%%e=%ld elCol[e]=%ld \n", e, elCol[e]));
*/ 
                
                assemble_row (el_Num, cb_numbers, mEl, nEl, 
                        colCount-fp, curRowIndex , k-fp, colRelIndex );

                PRLEVEL (1, ("%% psrc=%ld", psrc));
                rowRelIndex [curRowIndex] = -1;
                el_rowIndex [curRowIndex] = -1;
                curEl->nrowsleft --;
            } 
            else 
                listRowTuples [pdst++] = curTpl; //keeping the tuple
        }

        PRLEVEL (1, ("%%   pdst=%ld", pdst));
        curRowTupleList->numTuple = pdst;
        paru_print_element (paruMatInfo, el_ind);
   }
}

