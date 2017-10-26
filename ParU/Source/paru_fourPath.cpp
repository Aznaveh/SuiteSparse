/** =========================================================================  /
 * =======================  paru_fourPath ===================================  /
 * ========================================================================== */
/*! This function is the same level as paru_assemble 
 *      basically doing the rest of the work for non pivotal rows and columns
 */

#include "Parallel_LU.hpp"

void paru_fourPath (paru_matrix *paruMatInfo,
        Int el_ind,         // Index of current CB
        Int fp,             //fp of current front
        cholmod_common *cc)
{

    DEBUGLEVEL(0);
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
        PRLEVEL (0, ("1st: c =%ld  numTuple = %ld\n", c, numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;
            if(e < 0 || curColIndex < 0 ){ 
                PRLEVEL (0, ("removing e =%ld\n", e));
                paru_remove_colTuple (ColList, c, i);
                i--; numTuple--;
                continue;  
            }

            if (elCol [e] < elCMark) // an element never seen before
                elCol [e] = elCMark + 1;
            else{ 
                elCol [e]++; 
                continue;
            }

            /*  Update rowRelIndex	 */
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int *el_rowIndex = rowIndex_pointer (curEl); //pointers to row index
            Int *rowRelIndValid = rowRelIndVal (curEl);
            Int *rowRelIndex = relRowInd (curEl);
            Int *colRelIndex    = relColInd (curEl);
            // Updating row relative indices 
            PRLEVEL (1, ("elRow[%ld]=%ld", e, elRow [e]));  
            for (Int rEl = 0; rEl < mEl; rEl++)   
                rowRelIndex [rEl] = isRowInFront [el_rowIndex [rEl]] 
                    - rowMark;
//            *rowRelIndValid = f ;//current front

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
        PRLEVEL (1, ("r =%ld  numTuple = %ld\n", r, numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Int curRowIndex = curTpl.f;
            if(e < 0 || curRowIndex < 0){ 
                paru_remove_rowTuple (RowList, r, i);
                i--; numTuple--;
                continue;  
            }


            if (elRow [e] < elRMark) // an element never seen before
                elRow [e] = elRMark + 1;
            else{ 
                elRow [e]++; 
                continue;
            }
            /* Update colRelIndex	 */
             Element *curEl = elementList[e];
            Int nEl = curEl->ncols;
            Int *el_colIndex = colIndex_pointer (curEl); //pointers to row index
            Int *colRelIndValid = rowRelIndVal (curEl);
            Int *colRelIndex    = relColInd (curEl);
            // Updating row relative indices 
            PRLEVEL (1, ("elCol[%ld]=%ld", e, elCol [e]));  

            for (Int cEl = 0; cEl < nEl; cEl++)   
                colRelIndex [cEl] = isColInCBcolSet [el_colIndex [cEl]] 
                    - colMark;
//            *colRelIndValid = f ;//current front


        }
    }
    /**************************************************************************/

    Int *cb_colIndex = colIndex_pointer (curCB);
    Int *cb_rowIndex = rowIndex_pointer (curCB);
    Int *cb_rowRelIndex = relRowInd (curCB);
    Int *cb_rowRelIndValid = rowRelIndVal (curCB);
    Int *cb_colRelIndex    = relColInd (curCB);
    double *cb_numbers = numeric_pointer (curCB);


    /*****                 3rd path: assemble columns                    ******/
    /* it contains adding and deleting tuples*/
    for (Int k = 0; k < colCount; k++){
        Int c = CBColList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
        PRLEVEL (0, ("3rd: c =%ld  numTuple = %ld\n", c, numTuple));
        Int pdst = 0,psrc;
        for (Int psrc = 0; psrc < numTuple; psrc ++){
            Tuple curTpl = listColTuples [psrc];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;
            ASSERT (e >= 0);
            ASSERT (curColIndex >= 0);

            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;

            Int *el_colIndex = colIndex_pointer (curEl);
            Int *el_rowIndex = rowIndex_pointer (curEl);
            Int *rowRelIndex = relRowInd (curEl);
            Int *rowRelIndValid = rowRelIndVal (curEl);
            Int *colRelIndex    = relColInd (curEl);

            ASSERT (el_colIndex[curColIndex] == c);
            ASSERT (curColIndex < nEl);
            double *el_Num = numeric_pointer (curEl);
            PRLEVEL (1, ("  numTuple =%ld\n",   numTuple));
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));
            PRLEVEL (1, ("f =%ld\n", curColIndex));
            PRLEVEL (1, ("CB: %ld x %ld\n", rowCount, colCount));
            if (elRow [e] - elRMark == curEl->nrowsleft){
                PRLEVEL (0, ("psrc=%ld", psrc));
                assemble_col (el_Num+curColIndex*mEl,cb_numbers+k*colCount,
                        mEl, rowRelIndex);
                colRelIndex [curColIndex] = -1;
                el_colIndex [curColIndex] = -1;
                curEl->ncolsleft --;

                //  paru_remove_colTuple (ColList, c, i);
                //  i--; numTuple--;
            } 
            else 
                listColTuples [pdst++] = curTpl; //keeping the tuple

        }

        curColTupleList->numTuple = pdst;
        PRLEVEL (0, ("pdst=%ld\n", pdst));


        //adding tuple for CB
        Tuple T; T.e = el_ind; T.f= k;
        if( paru_add_colTuple (ColList, c, T, cc)){
            printf("Out of memory: add_colTuple for CB \n");
        }

    }



    /**************************************************************************/

    /*! TODO: 4th path: assemble rows */

}
