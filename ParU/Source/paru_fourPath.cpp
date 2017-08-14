/** =========================================================================  /
 * =======================  paru_fourPath ===================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

void paru_fourPath (paru_matrix *paruMatInfo,
        Int rowCount,
        Int colCount){

    DEBUGLEVEL(0);
    work_struct *Work =  paruMatInfo->Work;

    Int *fsRowList = Work->scratch; // fully summed row list
    Int *isRowInFront = Work->rowSize; 
    Int rowMark = Work->rowMark;

    Int *CBColList = Work -> scratch + 2*rowCount; //scratch=[fsRowList..ipiv..]
    Int *isColInCBcolSet = Work -> colSize;
    Int colMark = Work -> colMark;


    Int *elRow = Work -> elRow; 
    Int elRMark = Work -> elRMark;
    Int *elCol = Work -> elCol;
    Int elCMark = Work -> elCMark;


    Element **elementList = paruMatInfo->elementList;

    /*!  1st path: over non pivotal columns to count rows  */

    tupleList *ColList = paruMatInfo->ColList;
    for (Int c = 0; c < colCount; c++){
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
        PRLEVEL (1, ("c =%ld  numTuple = %ld\n", c, numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;
            Int *el_colIndex = (Int*)(curEl+1);  // pointers to element index 
            Int *el_rowIndex = el_colIndex + nEl;// pointers to row indices
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));
            for (Int rEl = 0; rEl < mEl; rEl++){
                Int curRow = el_rowIndex [rEl]; 
                PRLEVEL (1, ("curRow =%ld\n", curRow));
                if (elRow [curRow] < elRMark) // an element never seen before
                    elRow [curRow] = elRMark + 1;
                else 
                    elRow [curRow]++; 
            }
        }
    }

    /*! TODO: 2st path: over rows to count columns */
    tupleList *RowList = paruMatInfo->RowList;
    for (Int r = 0; r < rowCount; r++){
        tupleList *curColTupleList = &RowList[r];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
        PRLEVEL (1, ("r =%ld  numTuple = %ld\n", r, numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;
            Int *el_colIndex = (Int*)(curEl+1);  // pointers to element index 
            Int *el_rowIndex = el_colIndex + nEl;// pointers to row indices
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));
            for (Int cEl = 0; cEl < mEl; cEl++){
                Int curCol = el_colIndex [cEl]; 
                PRLEVEL (1, ("curCol =%ld\n", curCol));
                if (elCol [curCol] < elCMark) // an element never seen before
                    elCol [curCol] = elCMark + 1;
                else 
                    elCol [curCol]++; 
            }
        }
    }
    /*! TODO: 3st path: assemble columns   */
    /*! TODO: 4st path: assemble rows */



    PRLEVEL (0, ("Hi\n"));
}
