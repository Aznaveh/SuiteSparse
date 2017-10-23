/** =========================================================================  /
 * =======================  paru_fourPath ===================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

void paru_fourPath (paru_matrix *paruMatInfo,
        double *dest_el_numbers, //current CB that want to be assembled
        Int fp,             //fp of current front
        Int rowCount,
        Int colCount)
{

    DEBUGLEVEL(0);
    work_struct *Work =  paruMatInfo->Work;

    Int *isRowInFront = Work->rowSize; 
    Int rowMark = Work->rowMark;

    Int *CBColList = Work -> scratch + 2*rowCount; //scratch=[fsRowList..ipiv..]
    Int *isColInCBcolSet = Work -> colSize;
    Int colMark = Work -> colMark;

    // Couning how many rows/cols of an element is seen
    Int *elRow = Work -> elRow; 
    Int elRMark = Work -> elRMark;
    Int *elCol = Work -> elCol;
    Int elCMark = Work -> elCMark;


    Element **elementList = paruMatInfo->elementList;

    /*****  1st path: over non pivotal columns to count rows             ******/

    tupleList *ColList = paruMatInfo->ColList;
    for (Int k = 0; k < colCount; k++){
        Int c = CBColList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
        PRLEVEL (1, ("c =%ld  numTuple = %ld\n", c, numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;
            if(e < 0 || curColIndex < 0 ){ 
                paru_remove_colTuple (ColList, c, i);
                i--; numTuple--;
                continue;  
            }

            if (elCol [e] < elCMark) // an element never seen before
                elCol [e] = elCMark + 1;
            else 
                elCol [e]++; 
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
            else 
                elRow [e]++; 
        }
    }
    /*! TODO: 3st path: assemble columns  
     * it contains adding and deleting tuples*/
    for (Int k = 0; k < colCount; k++){
        Int c = CBColList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
        PRLEVEL (1, ("c =%ld  numTuple = %ld\n", c, numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;
            ASSERT (e > 0);
            ASSERT (curColIndex > 0);

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
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));

            if (elRow [e] - elRMark == curEl->nrowsleft);// if I can take the
                                                       //  whole column
                        assemble_col (dest_el_numbers+c*colCount,
                                el_Num+curColIndex*mEl,mEl, rowRelIndex);
        }
    }



    /**************************************************************************/

    /*! TODO: 4st path: assemble rows */

}
