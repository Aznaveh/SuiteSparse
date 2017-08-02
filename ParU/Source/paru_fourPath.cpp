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

    /*! TODO: 1st path: over non pivotal columns to count rows  */

    tupleList *ColList = paruMatInfo->ColList;
    for (Int c = 0; c < colCount; c++){
        tupleList *curColTupleList = &ColList[c];
    }
    /*! TODO: 2st path: over rows to count columns */
    /*! TODO: 3st path: assemble columns   */
    /*! TODO: 4st path: assemble rows */



    PRLEVEL (0, ("Hi\n"));
}
