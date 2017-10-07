/** =========================================================================  /
 * =======================  paru_assemble   =================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

#ifndef NDEBUG  // using STL for debugging
#include <iostream>
#include <algorithm>
#include <set>
#endif

/*! \brief assembling a front and updating correspoing elelment
 *
 *
 * \param  a list of tuples and the the tuple we want to add
 * \return 0 on sucess 
 */
void paru_assemble (
        paru_matrix *paruMatInfo,
        /* RowCol list/tuples and LUsym handle */
        Int f,
        /* front need to be assembled */
        cholmod_common *cc)

{
    DEBUGLEVEL(0);
    paru_symbolic *LUsym =  paruMatInfo->LUsym;

    Int m,n,nf;
    m = paruMatInfo-> m;
    n = paruMatInfo-> n;
    nf = paruMatInfo->LUsym->nf;
    Int *Super = LUsym->Super;
    /* ---------------------------------------------------------------------- */
    /* get the front F  */
    /* ---------------------------------------------------------------------- */


    PRLEVEL (0, ("Assemble Front %ld\n", f));
    /* pivotal columns Super [f] ... Super [f+1]-1 */
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];
    Int *Rp = LUsym->Rp;
    Int *Rj = LUsym->Rj;
    Int p1 = Rp [f];        /* Rj [p1:p2-1] = columns in F */
    Int p2 = Rp [f+1];
    Int fp = col2 - col1;   /* first fp columns are pivotal */ 
    Int fn = p2 - p1;          /* Upper bound number of columns of F */ 
    Element **elementList = paruMatInfo->elementList;
    work_struct *Work =  paruMatInfo->Work;

    // Couning how many rows/cols of an element is seen
    /*! TODO: grasping a prior CB sooner     */
    Int *elRow = Work -> elRow; 
    Int elRMark = Work -> elRMark;
    Int *elCol = Work -> elCol;
    Int elCMark = Work -> elCMark;

    PRLEVEL (0, ("fp=%ld pivotal columns:clo1=%ld...col2=%ld\n", 
                fp, col1, col2-1));
    PRLEVEL (1, ("Upper bound number of columns: Rj[%ld]=%ld ... Rj[%ld]=%ld\n", 
                p1, Rj [p1], p2, Rj [p2-1]));
    ASSERT (fp > 0 );
    ASSERT (fp <= fn );

    /* computing number of rows, set union */
    tupleList *ColList = paruMatInfo->ColList;
#if 0
    Int rowsP = 0;
    Int setSize = fn;  // it can be initialized with fm
    Int *rowSet = (Int*) paru_alloc (setSize, sizeof (Int), cc);
#endif

    Int rowCount= 0;
    Int *isRowInFront = Work->rowSize; 
    Int rowMark = Work->rowMark;
    if (rowMark < 0) {  // in rare case of overflow
        memset (isRowInFront, -1, m*sizeof(Int));
        rowMark = Work->rowMark = 0;
    }

    Int *fsRowList = Work->scratch; // fully summed row list
    PRLEVEL (1, ("fsRowList(scratch)=%p isRowInFront(all_initialized)=%p\n", 
                fsRowList, isRowInFront));

#ifndef NDEBUG /* chekcing first part of Work to be zero */
    for (Int i = 0; i < m; i++)  
        ASSERT ( isRowInFront [i] < rowMark);
#endif 


#ifndef NDEBUG
    std::set<Int> stl_rowSet;
    std::set<Int>::iterator it;
#endif 

    /**** 1 ******** finding set of rows in current front *********************/
    for (Int c = col1; c < col2; c++){
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
            Int *el_rowIndex = rowIndex_pointer (curEl); //pointers to row index
            Int *rowRelIndValid = rowRelIndVal (curEl);
            Int *rowRelIndex = relRowInd (curEl);


            /*! TODO: Is it possible to have columns already ate e<0 or f<0 */
            //if (curTpl.f<0) continue;

            //counting prior element's columns
            /*! TODO:Update Mark somewhere    !!! Here is an issue !!! */
            if (elCol [e] < elCMark) // an element never seen before
                elCol [e] = elCMark + 1;
            else { 
                elCol [e]++;    //keep track of number of cols
                continue;       // already have the set of rows
            }

            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));
            for (Int rEl = 0; rEl < mEl; rEl++){
                Int curRow = el_rowIndex [rEl]; 
                PRLEVEL (1, ("curRow =%ld\n", curRow));
                ASSERT (curRow < m ) ;
#ifndef NDEBUG
                stl_rowSet.insert (curRow);
#endif
                PRLEVEL (1, ("%p ---> isRowInFront [%ld]=%ld\n", 
                            isRowInFront+curRow, curRow, isRowInFront[curRow]));

                if (isRowInFront[curRow] < rowMark ){
                    PRLEVEL (1, ("curRow =%ld rowCount=%ld\n", 
                                curRow, rowCount));
                    fsRowList [rowCount] = curRow;
                    isRowInFront [curRow] = rowMark + rowCount++; 
                }
                ASSERT (rowCount <= m); 

                rowRelIndex [rEl] = isRowInFront [el_rowIndex [rEl]] 
                    - rowMark;

            }
            //            // Updating row relative indices 
            //            PRLEVEL (1, ("elRow[%ld]=%ld", e, elRow [e]));  
            //            for (Int rEl = 0; rEl < mEl; rEl++)   
            //                rowRelIndex [rEl] = isRowInFront [el_rowIndex [rEl]] 
            //                    - rowMark;
            *rowRelIndValid = f ;//current front


        }
    }

#ifndef NDEBUG /* Checking if pivotal rows are correct */
    Int p = 1;
    PRLEVEL (p, ("There are %ld rows in this front: \n", rowCount));
    for (Int i = 0; i < rowCount; i++)
        PRLEVEL (p, (" %ld", fsRowList [i]));
    PRLEVEL (p, ("\n"));
    Int stl_rowSize = stl_rowSet.size();
    if (rowCount != stl_rowSize){
        PRLEVEL (p, ("STL %ld:\n",stl_rowSize));
        for (it = stl_rowSet.begin(); it != stl_rowSet.end(); it++)
            PRLEVEL (p, (" %ld", *it));
        PRLEVEL (p, ("\nMy Set %ld:\n",rowCount));
        for (Int i = 0; i < rowCount; i++)
            PRLEVEL (p, (" %ld", fsRowList [i]));
        PRLEVEL (p, ("\n"));
    }
    ASSERT (rowCount == stl_rowSize );
    //   ASSERT (rowCount >= fp ); // otherwise it is a singular matrix

#endif 



    double *pivotalFront = 
        (double*) paru_calloc (rowCount*fp, sizeof (double), cc);

    if (pivotalFront == NULL ){
        printf ("Out of memory when tried to allocate for pivotal part %ld",f);
        return;
    }

    /**** 2 ********  assembling the pivotal part of the front ****************/
    /* 
     *                  
     *  curEl           nEl           
     *              6, 7, 3, 12
     *             ____________
     *          23 | X  Y  .  .     stored in memory like this:
     *      mEl 17 | X  Y  .  .     ...6, 7, 3, 10, 23, 17, 2, X, X, X, Y, Y, Y,
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

    /*  pivot assembly */
    /*! TODO: check if any row/col nulified     */
    for (Int c = col1; c < col2; c++){
        tupleList *curTupleList = &ColList[c];
        Int numTuple = curTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *l = curTupleList->list;
        PRLEVEL (1, ("c =%ld numTuple = %ld\n", c, numTuple));

        Int colIndexF = c - col1;  // relative column index

        for (Int i = 0; i < numTuple; i++){

            Tuple curTpl = l [i];
            Int e = curTpl.e;

            // Assembly of column f of e in colIndexF
            PRLEVEL (1, ("col=%ld, (%ld,%ld)\n", c, e, f));
            //!! WRONG do not delete tuples here
       //     FLIP (curTpl.e); //Nullifying tuple  /*! TODO: Deleting tuple     */
       //     curTupleList->numTuple--;

            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;

            Int *el_colIndex = colIndex_pointer (curEl);
            Int *el_rowIndex = rowIndex_pointer (curEl);
            Int *rowRelIndex = relRowInd (curEl);
            Int *rowRelIndValid = rowRelIndVal (curEl);
            Int *colRelIndex    = relColInd (curEl);

            Int curColIndex = curTpl.f;
            ASSERT (el_colIndex[curColIndex] == c);


            /*! TODO: I can change number of left cols/rows in the end       */
            //            curEl->ncolsleft--;
            PRLEVEL (1, ("curColIndex =%ld\n", curColIndex));

            ASSERT (curColIndex < nEl);
            double *el_Num = numeric_pointer (curEl);
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));

            assemble_col (pivotalFront+colIndexF*rowCount, 
                    el_Num +curColIndex*mEl, mEl, rowRelIndex);

            FLIP(el_colIndex[curColIndex]); //Marking column as assembled
            colRelIndex [curColIndex] = -1;

#if 0
            //This part has been moved to assemble_col function
                  for (Int rEl = 0; rEl < mEl; rEl++){   
                      Int curRow = el_rowIndex [rEl]; 

                      PRLEVEL (1, ("curRow =%ld\n", curRow));
                      ASSERT (curRow < m ) ;
                      ASSERT (isRowInFront [curRow] != -1);

                      PRLEVEL (1, ("rowRelIndValid = %ld, f=%ld,\
                                  rowRelIndex[%ld]= %ld,##%ld\n ", 
                                  *rowRelIndValid, f, rEl, rowRelIndex [rEl], 
                                  isRowInFront [curRow] - rowMark));
                      // relative row index
                      Int rowIndexF = rowRelIndex [rEl];
                      ASSERT (rowIndexF == isRowInFront [curRow] - rowMark);

                      PRLEVEL (1, ("rowIndexF = %ld\n", rowIndexF));
                      PRLEVEL (1, (" colIndexF*rowCount + rowIndexF=%ld\n",
                                  colIndexF*rowCount + rowIndexF));
                      ASSERT ( colIndexF*rowCount + rowIndexF < rowCount * fp);
                      ASSERT ( curColIndex*mEl + rEl < mEl*nEl);
                      pivotalFront [colIndexF*rowCount + rowIndexF] += 
                          el_Num [ curColIndex*mEl + rEl];
                  }
#endif
        }
    }

#ifndef NDEBUG  // Printing the pivotal front
    p = 1;
    PRLEVEL (p, ("x\t"));
    for (Int c = col1; c < col2; c++) {
        PRLEVEL (p, ("%ld\t", c));
    }
    PRLEVEL (p, ("\n"));
    for (Int r = 0; r < rowCount; r++){
        PRLEVEL (p, ("%ld\t", fsRowList [r]));
        for (Int c = col1; c < col2; c++){
            PRLEVEL (p, (" %3.1lf\t", pivotalFront [(c-col1)*rowCount + r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
    /**** 3 ********  factorizing the fully summed part of the matrix    
     *****  a set of pivot is found in this part that is crucial to assemble  **/
    PRLEVEL (1, ("rowCount =%ld\n", rowCount));
    /* using the rest of scratch for permutation; Not sure about 1  */
    int *ipiv =(int *) (Work->scratch+rowCount);
    paru_factorize (pivotalFront, rowCount, fp, ipiv, cc );

    /* To this point fully summed part of the front is computed and L and U    /  
     *  The next part is to find columns of nonfully summed then rows
     *  the rest of the matrix and doing TRSM and GEMM,                       */
    if (ipiv [0] < 0){
        printf ("Singular Matrix\n");
        return;
    }
#ifndef NDEBUG  // Printing the permutation
    p = 0;
    PRLEVEL (p, ("permutation:\n"));
    for (int i = 0; i < fp; i++){
        PRLEVEL (p, ("ipiv[%d] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (p, ("\n"));
#endif


    /*  Searching for columns */
    Int *isColInCBcolSet = Work -> colSize;
    Int colMark = Work -> colMark;
    if (colMark < 0) {  // in rare case of overflow
        memset (isRowInFront, -1, n*sizeof(Int));
        colMark = Work-> colMark = 0;
    }
    Int *CBColList = Work -> scratch + 2*rowCount;//scratch=[fsRowList..ipiv..]
    Int colCount = 0;

#ifndef NDEBUG
    std::set<Int> stl_colSet;
#endif  

    /**** 4 ******** finding set of non pivotal cols in current front *********/
    tupleList *RowList = paruMatInfo->RowList;
    for (Int i = 0; i < fp; i++){
        Int curFsRowIndex =(Int) ipiv [i]; //current fully summed row index
        PRLEVEL (1, ("4: curFsRowIndex = %ld\n", curFsRowIndex));
        ASSERT (curFsRowIndex < m);
        Int curFsRow = fsRowList [curFsRowIndex];
        PRLEVEL (1, ("curFsRow =%ld\n", curFsRow));
        tupleList *curRowTupleList = &RowList [curFsRowIndex];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listRowTuples = curRowTupleList->list;
        PRLEVEL (1, ("numTuple = %ld\n", numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;
            Int *el_rowIndex = rowIndex_pointer (curEl);
            Int *el_colIndex = colIndex_pointer (curEl);
            Int *colRelIndValid = colRelIndVal (curEl);
            Int *colRelIndex = relColInd (curEl);

            //counting prior element's rows
            if (elRow [e] < elRMark) {// an element never seen before
                elRow [e] = elRMark + 1;
#ifndef NDEBUG
                if ( elCol [e] >= elCMark )
                    /*! TODO: Write a phase to assemble them     */
                    PRLEVEL (0, ("element %ld must be eaten wholly\n",e));
                //And the rest of e is in U part 
#endif
            }
            else { // must not happen anyway; it depends on changing strategy
                elRow [e]++; 
                continue;
            }

            PRLEVEL (1, ("element= %ld  nEl =%ld \n",e, nEl));
            for (Int cEl = 0; cEl < nEl; cEl++){
                Int curCol = el_colIndex [cEl]; 
                PRLEVEL (1, ("curCol =%ld\n", curCol));
                ASSERT (curCol < n);
                /*! TODO: implement this part better     */
                if (curCol >= col1 && curCol < col2) //if in pivotal front
                    continue;
#ifndef NDEBUG
                stl_colSet.insert (curCol);
#endif
                PRLEVEL (1, ("%p ---> isColInCBcolSet[%ld]=%ld\n", 
                            isColInCBcolSet+curCol, curCol,
                            isColInCBcolSet[curCol]));
                /*! TODO: Check for elements too in Pass 1 and 3  */
                if (isColInCBcolSet [curCol] < colMark  ){
                    PRLEVEL (1, ("curCol = %ld colCount=%ld\n", 
                                curCol, colCount));
                    CBColList [colCount] = curCol;
                    isColInCBcolSet [curCol] = colMark + colCount++; 
                }
                ASSERT (colCount <= n);

                colRelIndex [cEl] = isColInCBcolSet [el_colIndex [cEl]] 
                    - colMark;
            }
        }
    }

#ifndef NDEBUG /* Checking if columns are correct */

    p = 1;
    PRLEVEL (p, ("There are %ld columns in this contribution block: \n",
                colCount));
    for (Int i = 0; i < colCount; i++)
        PRLEVEL (p, (" %ld", CBColList [i]));
    PRLEVEL (p, ("\n"));
    Int stl_colSize = stl_colSet.size();
    if (colCount != stl_colSize){
        PRLEVEL (p, ("#######################\n"));
        PRLEVEL (p, ("STL %ld:\n",stl_colSize));
        for (it = stl_rowSet.begin(); it != stl_colSet.end(); it++)
            PRLEVEL (p, (" %ld", *it));
        PRLEVEL (p, ("\nMy Set %ld:\n",colCount));
        for (Int i = 0; i < colCount; i++)
            PRLEVEL (p, (" %ld", CBColList [i]));
        PRLEVEL (p, ("\n"));
    }
    ASSERT (colCount == stl_colSize );

#endif 


#ifdef NotUsingMark
    /*Not used now, I am using rowMark to avoid this*/
    /* setting W for next iteration
       for (Int i = 0; i < rowCount; i++){
       Int curRow = fsRowList [i];
       ASSERT (curRow < m );
       ASSERT (isRowInFront [curRow] != -1);
       isRowInFront  [curRow] = -1;
       } */
#endif

    if(colCount !=0 && rowCount != 0)
        paru_fourPath (paruMatInfo, rowCount, colCount);

    double *uPart = 
        (double*) paru_calloc (fp*colCount, sizeof (double), cc);
    if ( uPart == NULL ){
        printf ("Out of memory when tried to allocate for U part %ld",f);
        return;
    }
    /*! TODO: Bug found ipiv changes!! Decument Work data structure*/
#ifndef NDEBUG  // Printing the permutation
    p = 1;
    PRLEVEL (p, ("permutation:\n"));
    for (int i = 0; i < fp; i++){
        PRLEVEL (p, ("ipiv[%d] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (p, ("\n"));
#endif


    /**** 5 ** assemble U part         Row by Row                          ****/ 
    for (Int i = 0; i < fp; i++){
        Int curFsRowIndex =(Int) ipiv [i]; //current fully summed row index
        PRLEVEL (1, ("5(i=%ld): curFsRowIndex = %ld\n",i, curFsRowIndex));
        ASSERT (curFsRowIndex < m);
        Int curFsRow = fsRowList [curFsRowIndex];
        PRLEVEL (1, ("curFsRow =%ld\n", curFsRow));
        tupleList *curRowTupleList = &RowList [curFsRowIndex];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        ASSERT (numTuple <= m);
       Tuple *listRowTuples = curRowTupleList->list;
       PRLEVEL (0, ("numTuple = %ld\n", numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;
            Int *el_rowIndex = rowIndex_pointer (curEl);
            Int *el_colIndex = colIndex_pointer (curEl);
            Int *colRelIndValid = colRelIndVal (curEl);
            Int *colRelIndex = relColInd (curEl);

            Int curRowIndex = curTpl.f;
            ASSERT (el_rowIndex[curRowIndex] == curFsRowIndex);
            ASSERT (curRowIndex < mEl)
            PRLEVEL (1, ("curColIndex =%ld\n", curRowIndex));

            double *el_Num = numeric_pointer (curEl);
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));



            /// Grab the row and add them to the u part row by row
//*            //counting prior element's rows
//*            if (elRow [e] < elRMark) {// an element never seen before
//*                elRow [e] = elRMark + 1;
//*           }
//*            else { // must not happen anyway; it depends on changing strategy
//*                elRow [e]++; 
//*                continue;
//*            }
//*
//*            PRLEVEL (0, ("element= %ld  nEl =%ld \n",e, nEl));
//*            for (Int cEl = 0; cEl < nEl; cEl++){
//*                Int curCol = el_colIndex [cEl]; 
//*                PRLEVEL (1, ("curCol =%ld\n", curCol));
//*                ASSERT (curCol < n);
//*                /*! TODO: implement this part better     */
//*                if (curCol >= col1 && curCol < col2) //if in pivotal front
//*                    continue;
//*               PRLEVEL (1, ("%p ---> isColInCBcolSet[%ld]=%ld\n", 
//*                            isColInCBcolSet+curCol, curCol,
//*                            isColInCBcolSet[curCol]));
//*                /*! TODO: Check for elements too in Pass 1 and 3  */
//*                if (isColInCBcolSet [curCol] < colMark  ){
//*                    PRLEVEL (1, ("curCol = %ld colCount=%ld\n", 
//*                                curCol, colCount));
//*                    CBColList [colCount] = curCol;
//*                    isColInCBcolSet [curCol] = colMark + colCount++; 
//*                }
//*                ASSERT (colCount <= n);
//*
//*                colRelIndex [cEl] = isColInCBcolSet [el_colIndex [cEl]] 
//*                    - colMark;
//*            }

        }
    }



///////////////////////////////////////////////////////////////////////////////
    Int *snM = LUsym->super2atree;
    Int eli = snM [f]; // Element index of the one that is going to be assembled
    if (fp < rowCount ){ // otherwise nothing will remain of this front
        Element *el = elementList[eli] = paru_create_element (rowCount-fp,
                colCount, 0 ,cc);
        if ( el == NULL ){
            printf ("Out of memory when tried to allocate current CB %ld",eli);
            return;
        }
    }


    /*! TODO:TRSM and DGEMM can be here*/
    /**** 6,7 ** Count number of rows and columsn of prior CBs to asslemble ***/ 
    /*
     * I need to either assemble rows or Columns                            ***/





    Work->rowMark += rowCount;
    rowMark = Work -> rowMark;

    Work->colMark += colCount;
    colMark = Work -> colMark;


#ifndef NDEBUG /* chekcing isRowInFront to be zero */
    for (Int i = 0; i < m; i++)
        ASSERT ( isRowInFront [i] < rowMark);
#endif


    /*! TODO: This should be stored somewhere     */
    paru_free (rowCount*fp, sizeof (double), pivotalFront, cc);
    paru_free (fp*colCount,  sizeof (double), uPart, cc);

#if 0
    paru_free (setSize, sizeof (Int), rowSet, cc);
#endif
}
