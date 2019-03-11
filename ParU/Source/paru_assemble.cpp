/** =========================================================================  /
 * =======================  paru_assemble   =================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

#ifndef NDEBUG  // using STL for debugging
#include <iostream>
#include <algorithm>
#include <set>
#endif

/*! @brief assembling a front and updating correspoing elelment
 *  @author Aznaveh
 *
 *
 * @param  a list of tuples and the the tuple we want to add
 * @return 0 on sucess 
 */
void paru_assemble (
        paru_matrix *paruMatInfo,
        /* RowCol list/tuples and LUsym handle */
        Int f,
        /* front need to be assembled */
        cholmod_common *cc)

{
    DEBUGLEVEL(0);
    /* 
     * -2 Print Nothing
     * -1 Just Matlab
     *  0 Detailed
     *  > 0 Everything
     */
    paru_symbolic *LUsym =  paruMatInfo->LUsym;

    Int m,n,nf;
    m = paruMatInfo-> m;
    n = paruMatInfo-> n;
    nf = paruMatInfo->LUsym->nf;
    Int *Super = LUsym->Super;
    /* ---------------------------------------------------------------------- */
    /* get the front F  */
    /* ---------------------------------------------------------------------- */


    PRLEVEL (-1, ("%%~~~~~~~  Assemble Front %ld start ~~~~\n", f));
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

    Int *elRow = Work -> elRow; 
    Int elRMark = Work -> elRMark;
    Int *elCol = Work -> elCol;
    Int elCMark = Work -> elCMark;

    PRLEVEL (0, ("%% fp=%ld pivotal columns:clo1=%ld...col2=%ld\n", 
                fp, col1, col2-1));
    PRLEVEL (1, ("%%Upper bound number of columns: Rj[%ld]=%ld ... Rj[%ld]=%ld\n", 
                p1, Rj [p1], p2, Rj [p2-1]));
    ASSERT (fp > 0 );
    ASSERT (fp <= fn );

    /* computing number of rows, set union */
    tupleList *ColList = paruMatInfo->ColList;

    Int rowCount= 0;
    Int *isRowInFront = Work->rowSize; 
    Int rowMark = Work->rowMark;
    PRLEVEL (-1, ("rowMark=%ld;\n", rowMark));

    if (rowMark < 0) {  // in rare case of overflow
        memset (isRowInFront, -1, m*sizeof(Int));
        rowMark = Work->rowMark = 0;
    }

    Int *fsRowList = Work->scratch; // fully summed row list
    PRLEVEL (1, ("%%fsRowList(scratch)=%p isRowInFront(all_initialized)=%p\n", 
                fsRowList, isRowInFront));

#ifndef NDEBUG /* chekcing first part of Work to be zero */
    for (Int i = 0; i < m; i++){  
        if ( isRowInFront [i] >= rowMark)
            PRLEVEL (1, ("%%rowMark = %ld, isRowInFront[%ld] = %ld\n", 
                        rowMark ,i,
                        isRowInFront [i]));
        ASSERT ( isRowInFront [i] < rowMark);
    }
#endif 


#ifndef NDEBUG
    std::set<Int> stl_rowSet;
    std::set<Int>::iterator it;
#endif 
    // Initializing relative index validation flag of current front
    paru_init_rel (paruMatInfo, f);
    Int time_f = paruMatInfo->time_stamp[f];

    /**** 1 ******** finding set of rows in current front *********************/
    for (Int c = col1; c < col2; c++){
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;

#ifndef NDEBUG            
        Int p = 1;
        PRLEVEL (p, ("%%c =%ld  numTuple = %ld\n", c, numTuple));
        if (p <= 0 )
            paru_print_tupleList (ColList, c);
#endif
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%%e =%ld  curColIndex = %ld\n", e, curColIndex));

            /*! TODO: Never negate e or f for now... keep it?	 */
            if(e < 0 || curColIndex < 0 ) continue;  //already deleted

            Element *el = elementList[e];
            Int mEl = el->nrows;
            Int *el_rowIndex = rowIndex_pointer (el); //pointers to row index
            Int *rowRelIndex = relRowInd (el);
            Int *colRelIndex    = relColInd (el);
            Int *el_colIndex = colIndex_pointer (el);

            PRLEVEL (1, ("%%point to col = %ld\n", el_colIndex[curColIndex]));
            /*! TODO: Keep the tuple or delete it?	 */
            if (el_colIndex [curColIndex]< 0 ) continue; // already assembled

            ASSERT (el_colIndex[curColIndex] == c);


            colRelIndex [curTpl.f] = c - col1; //Initialzing relative index
            // neede for row assembly
            // function

            PRLEVEL (1, ("%%1:element= %ld elCol=%ld elCMark=%ld \n",
                        e, elCol[e], elCMark));


            if(el->rValid !=  time_f){  // an element never seen before
                el->rValid = time_f;
#ifndef NDEBUG            
                if (el->rValid >  time_f )
                    PRLEVEL (0, ("%%time_f =%ld  rVal= %ld\n",
                                time_f , el->rValid));
#endif               
                ASSERT(el->rValid <= time_f);
            }
            else { 
                elCol [e]--;    //keep track of number of cols
                PRLEVEL (1, ("%%  element= %ld is seen before \n",e));
                //counting prior element's columns
                continue;       // already have the set of rows
            }

            PRLEVEL (1, ("%%2:element= %ld  mEl =%ld \n",e, mEl));
            //Set union of all the rows of the current element
            for (Int rEl = 0; rEl < mEl; rEl++){
                Int curRow = el_rowIndex [rEl]; 
                PRLEVEL (1, ("%%@@curRow =%ld rEl=%ld\n", curRow, rEl));
                Int p=1;
                if (p <= 0) paru_print_element (paruMatInfo, e);
                if (curRow < 0 ) continue; // that row has already deleted
                PRLEVEL (1, ("%%**curRow =%ld rEl=%ld\n", curRow, rEl));
                ASSERT (curRow < m ) ;
#ifndef NDEBUG
                stl_rowSet.insert (curRow);
#endif
                PRLEVEL (1, ("%% %p ---> isRowInFront [%ld]=%ld\n", 
                            isRowInFront+curRow, curRow, isRowInFront[curRow]));

                if (isRowInFront[curRow] < rowMark ){  
                    // Adding curRow to the set
                    PRLEVEL (1, ("%%curRow =%ld rowCount=%ld\n", 
                                curRow, rowCount));
                    fsRowList [rowCount] = curRow;
                    rowRelIndex [rEl] = rowCount ;
                    PRLEVEL (1, ("%%1st: rowRelIndex[%ld] = %ld\n",
                                rEl, rowCount ));
                    isRowInFront [curRow] = rowMark + rowCount++; 
                }
                else{
                    rowRelIndex [rEl] = isRowInFront [curRow] - rowMark;
                    PRLEVEL (1, ("%%N1st: rowRelIndex[%ld] = %ld\n",
                                rEl, rowCount ));
                }
                ASSERT (rowCount <= m); 
            }
        }
    }

#ifndef NDEBUG /* Checking if pivotal rows are correct */
    Int p = 1;
    PRLEVEL (p, ("%%There are %ld rows in this front: \n", rowCount));
    for (Int i = 0; i < rowCount; i++)
        PRLEVEL (p, ("%% %ld", fsRowList [i]));
    PRLEVEL (p, ("\n"));
    Int stl_rowSize = stl_rowSet.size();
    if (rowCount != stl_rowSize){
        PRLEVEL (p, ("%% STL %ld:\n",stl_rowSize));
        for (it = stl_rowSet.begin(); it != stl_rowSet.end(); it++)
            PRLEVEL (p, ("%% %ld", *it));
        PRLEVEL (p, ("\n%%My Set %ld:\n",rowCount));
        for (Int i = 0; i < rowCount; i++)
            PRLEVEL (p, ("%% %ld", fsRowList [i]));
        PRLEVEL (p, ("\n"));
    }
    ASSERT (rowCount == stl_rowSize );
    //   ASSERT (rowCount >= fp ); // otherwise it is a singular matrix

#endif 



    double *pivotalFront = 
        (double*) paru_calloc (rowCount*fp, sizeof (double), cc);

    if (pivotalFront == NULL ){
        printf ("%% Out of memory when tried to allocate for pivotal part %ld",
                f);
        return;
    }

    paru_fac *LUs =  paruMatInfo->partial_LUs;
    LUs[f].m =rowCount;
    LUs[f].n=fp;
    ASSERT (LUs[f].p == NULL);
    LUs[f].p = pivotalFront;

    /**** 2 ********  assembling the pivotal part of the front ****************/
    /* 
     *                  
     *  el           nEl           
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
    for (Int c = col1; c < col2; c++){
        tupleList *curTupleList = &ColList[c];
        Int numTuple = curTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *l = curTupleList->list;
        PRLEVEL (1, ("%% c =%ld numTuple = %ld\n", c, numTuple));

        Int colIndexF = c - col1;  // relative column index

        for (Int i = 0; i < numTuple; i++){

            Tuple curTpl = l [i];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;

            // Assembly of column curColIndex of e in colIndexF

            Element *el = elementList[e];
            Int mEl = el->nrows;
            Int nEl = el->ncols;

            Int *el_colIndex = colIndex_pointer (el);
            Int *el_rowIndex = rowIndex_pointer (el);
            Int *rowRelIndex = relRowInd (el);
            Int *colRelIndex    = relColInd (el);
            if (el_colIndex [curColIndex]< 0 ) continue;
            ASSERT (el_colIndex[curColIndex] == c);


            PRLEVEL (1, ("%% curColIndex =%ld\n", curColIndex));

            ASSERT (curColIndex < nEl);
            double *el_Num = numeric_pointer (el);
#ifndef NDEBUG // print the element which is going to be assembled from
            p = 1;
            PRLEVEL (p, ("%% col=%ld, element=%ld,curColIndex=%ld\n", c, e,
                        curColIndex));
            PRLEVEL (p, ("%% ASSEMBL element= %ld  mEl =%ld ",e, mEl));
            PRLEVEL (p, ("%% into column %ld of current front\n",colIndexF ));
            PRLEVEL (p, ("%%assembling from col %ld", curColIndex ));
            if (p <= 0)
                paru_print_element (paruMatInfo, e);
#endif

            assemble_col (el_Num +curColIndex*mEl,
                    pivotalFront+colIndexF*rowCount,
                    mEl, rowRelIndex);

            //FLIP(el_colIndex[curColIndex]); //marking column as assembled
            el_colIndex[curColIndex] = -1;
            colRelIndex [curColIndex] = -1;
            el->ncolsleft--;     // After each assembly Az: Added later


#ifndef NDEBUG  // Printing the pivotal front
            p = 1;
            PRLEVEL (p, ("%%Inside the assembly loop \n"));
            PRLEVEL (p, ("%%x =  \t"));
            for (Int c = col1; c < col2; c++) {
                PRLEVEL (p, ("%% %ld\t\t", c));
            }
            PRLEVEL (p, (" ;\n"));
            for (Int r = 0; r < rowCount; r++){
                PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
                for (Int c = col1; c < col2; c++){
                    PRLEVEL (p, (" %2.5lf\t", 
                                pivotalFront [(c-col1)*rowCount + r]));
                }
                PRLEVEL (p, ("\n"));
            }
#endif

        }
    }

#ifndef NDEBUG  // Printing the pivotal front
    p = 0;
    PRLEVEL (p, ("%% Before pivoting\n"));
    PRLEVEL (p, ("%% x =  \t"));
    for (Int c = col1; c < col2; c++) {
        PRLEVEL (p, ("%ld\t\t", c));
    }
    PRLEVEL (p, (" ;\n"));
    for (Int r = 0; r < rowCount; r++){
        PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
        for (Int c = col1; c < col2; c++){
            PRLEVEL (p, (" %2.5lf\t", pivotalFront [(c-col1)*rowCount + r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
    /**** 3 ********  factorizing the fully summed part of the matrix        ***
     *****  a set of pivot is found in this part that is crucial to assemble **/
    PRLEVEL (1, ("%% rowCount =%ld\n", rowCount));

#ifndef NDEBUG  // Printing the list of rows
    p = 1;
    PRLEVEL (p, ("%% Befor factorization (inside assemble): \n"));
    for (int i = 0; i < rowCount; i++){
        PRLEVEL (p, ("%% fsRowList [%d] =%d\n",i, fsRowList [i]));
    }
    PRLEVEL (p, ("\n"));
#endif



    /* using the rest of scratch for permutation; Not sure about 1  */
    BLAS_INT *ipiv = (BLAS_INT*) (Work->scratch+rowCount);
    //Int fac = paru_dgetrf (pivotalFront, fsRowList, rowCount, fp, ipiv);
    Int fac = paru_factorize(pivotalFront, fsRowList, rowCount, 
            fp, paruMatInfo);

    /* To this point fully summed part of the front is computed and L and U    /  
     *  The next part is to find columns of nonfully summed then rows
     *  the rest of the matrix and doing TRSM and GEMM,                       */
    if (fac < 0){
        printf ("%% m > n\n");
        exit(0);
        return;
    }
    if (ipiv [0] < 0){
        printf ("%% Singular Matrix\n");
        return;
    }
#ifndef NDEBUG  // Printing the list of rows
    p = 1;
    PRLEVEL (p, ("%% After factorization (inside assemble): \n"));
    for (int i = 0; i < rowCount; i++){
        PRLEVEL (p, ("%% fsRowList [%d] =%d\n",i, fsRowList [i]));
    }
    PRLEVEL (p, ("\n"));
#endif


#ifndef NDEBUG  // Printing the permutation
    p = 1;
    PRLEVEL (p, ("%% pivotal rows:\n"));
    for (int i = 0; i < fp; i++){
        PRLEVEL (p, ("%% fsRowList[%d] =%d\n",i, fsRowList[i]));
    }
    PRLEVEL (p, ("%% =======\n"));
    for (int i = fp; i < rowCount; i++){
        PRLEVEL (p, ("%% fsRowList[%d] =%d\n",i, fsRowList[i]));
    }
    PRLEVEL (p, ("\n"));
#endif

#ifndef NDEBUG  // Printing the pivotal front
    p = -1;
    PRLEVEL (p, ("%%L part:\n"));

    //col permutatin
    PRLEVEL (p, ("cols{%ld} = [",f+1));
    for (Int c = col1; c < col2; c++) {
        PRLEVEL (p, ("%ld ", c+1));
    }
    PRLEVEL (p, ("];\n"));

    //row permutatin
    PRLEVEL (p, ("rows{%ld} = [",f+1));
    for (Int r = 0; r < rowCount; r++)
        PRLEVEL (p, ("%ld ", fsRowList [r]+1)); //Matlab is base 1
    PRLEVEL (p, ("];\n"));

    //inv row permutatin

    PRLEVEL (p, ("Luf{%ld}= [",f+1));
    for (Int r = 0; r < rowCount; r++){
        PRLEVEL (p, (" "));
        for (Int c = col1; c < col2; c++){
            PRLEVEL (p, (" %.16lf ", pivotalFront [(c-col1)*rowCount + r]));
        }
        PRLEVEL (p, (";\n   "));
    }
    PRLEVEL (p, ("];\n"));
    //just in cases that there is no U
    PRLEVEL (p, ("Us{%ld} =[];\n", f+1));
    PRLEVEL (p, ("Ucols{%ld}=[];\n",f+1));
    PRLEVEL (p, ("Urows{%ld}=[];\n",f+1));
#endif


    /**** 4 ******** finding set of non pivotal cols in current front *********/
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

    tupleList *RowList = paruMatInfo->RowList;
    for (Int i = 0; i < fp; i++){
        Int curFsRowIndex =(Int) i; //current fully summed row index
        Int curFsRow = fsRowList [i];
        PRLEVEL (1, ("%% 4: curFsRowIndex = %ld\n", curFsRowIndex));
        PRLEVEL (1, ("%% curFsRow =%ld\n", curFsRow));
        //tupleList *curRowTupleList = &RowList [curFsRowIndex]; //BUG DETECTED
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


    /**** 5 ** assemble U part         Row by Row                          ****/ 

    double *uPart = 
        (double*) paru_calloc (fp*colCount, sizeof (double), cc);
    if ( uPart == NULL ){
        printf ("%% Out of memory when tried to allocate for U part %ld",f);
        return;
    }

    paru_fac *Us =  paruMatInfo->partial_Us;
    Us[f].m=fp;
    Us[f].n =colCount;
    ASSERT (Us[f].p == NULL);
    Us[f].p = uPart;


    for (Int i = 0; i < fp; i++){
        Int curFsRowIndex = i; //current fully summed row index
        Int curFsRow = fsRowList [curFsRowIndex];
        PRLEVEL (1, ("%% curFsRow =%ld\n", curFsRow));
        tupleList *curRowTupleList = &RowList [curFsRow];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        ASSERT (numTuple <= m);
        Tuple *listRowTuples = curRowTupleList->list;
        PRLEVEL (1, ("%% numTuple = %ld\n", numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;
            Int curRowIndex = curTpl.f;

            Element *el = elementList[e];
            Int mEl = el->nrows;
            Int nEl = el->ncols;
            Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_colIndex = colIndex_pointer (el);
            Int *rowRelIndex = relRowInd (el);
            Int *colRelIndex = relColInd (el);

            if(el_rowIndex[curRowIndex] < 0 ) continue; 

            PRLEVEL (1, ("%% curFsRowIndex =%ld\n", curFsRowIndex));
            ASSERT (el_rowIndex[curRowIndex] == curFsRow);
            ASSERT (curRowIndex < mEl);
            PRLEVEL (1, ("%% curColIndex =%ld\n", curRowIndex));

            double *el_Num = numeric_pointer (el);
            PRLEVEL (1, ("%% element= %ld  nEl =%ld \n",e, nEl));

            assemble_row (el_Num, uPart, mEl, nEl, fp, 
                    curRowIndex, curFsRowIndex, colRelIndex);

            //FLIP(el_rowIndex[curRowIndex]); //marking column as assembled
            el_rowIndex[curRowIndex] = -1;
            rowRelIndex [curRowIndex] = -1;
            el->nrowsleft--;  // After each assembly Az: Added later

        }
    }

#ifndef NDEBUG  // Printing the  U part
    p = 1;
    PRLEVEL (p, ("%% U part Before TRSM: %ld x %ld\n", fp, colCount));
    PRLEVEL (p, ("%% U\t"));
    for (Int i = 0; i < colCount; i++){
        PRLEVEL (p, ("%ld\t\t", CBColList[i]));
    }
    PRLEVEL (p, ("\n"));
    for (Int i = 0; i < fp; i++){
        PRLEVEL (p, ("%% %ld\t",  fsRowList [i]));
        for (Int j = 0; j < colCount; j++){
            PRLEVEL (p, (" %2.5lf\t", uPart[j*fp+i]));
        }
        PRLEVEL (p, ("\n"));
    }

#endif


    /**** 6 ****                     TRSM and DGEMM                         ***/ 

    paru_trsm(pivotalFront , uPart, fp, rowCount, colCount);

#ifndef NDEBUG  // Printing the  U part
    p = -1;
    PRLEVEL (p, ("%% rowCount=%ld;\n",rowCount));
    PRLEVEL (p, ("%% U part After TRSM: %ld x %ld\n", fp, colCount));

    PRLEVEL (p, ("Ucols{%ld} = [",f+1));
    for (Int i = 0; i < colCount; i++){
        PRLEVEL (p, ("%ld ", CBColList[i]+1));
    }
    PRLEVEL (p, ("];\n"));


    PRLEVEL (p, ("Urows{%ld} = [",f+1));
    for (Int i = 0; i < fp; i++)
        PRLEVEL (p, ("%ld ",  fsRowList [i]+1));
    PRLEVEL (p, ("];\n"));


    PRLEVEL (p, ("Us{%ld} = [",f+1));

    for (Int i = 0; i < fp; i++){
        for (Int j = 0; j < colCount; j++){
            PRLEVEL (p, (" %.16lf ", uPart[j*fp+i]));
        }
        PRLEVEL (p, (";\n    "));
    }
    PRLEVEL (p, ("];\n"));
#endif


    Int *snM = LUsym->super2atree;
    Int eli = snM [f]; // Element index of the one that is going to be assembled
    Element *curEl;
    PRLEVEL (1, ("%% rowCount=%ld, colCount=%ld, fp=%ld\n",
                rowCount, colCount, fp));
    PRLEVEL (1, ("%% curEl is %ld by %ld\n",rowCount-fp,colCount));
    if (fp <= rowCount ){ 
        curEl = elementList[eli] = paru_create_element (rowCount-fp,
                colCount, 0 ,cc); // allocating an un-initialized part of memory
        // While insided the DGEMM BETA == 0
        if ( curEl == NULL ){
            printf ("%% Out of memory when tried to allocate current CB %ld",
                    eli);
            return;
        }
        PRLEVEL (1, ("%% curEl =%p\n", curEl));
    }
    // Initializing curEl global indices
    Int *el_colIndex = colIndex_pointer (curEl);
    for (Int i = 0; i < colCount; ++ i) {
        el_colIndex [i] = CBColList[i];
        Tuple colTuple;
    }
    Int *el_rowIndex = rowIndex_pointer (curEl);
    for (Int i = fp; i < rowCount; ++ i) {
        Int locIndx = i-fp; 
        el_rowIndex [locIndx] = fsRowList[i];
        PRLEVEL (1, ("%% el_rowIndex [%ld] =%ld\n",
                    locIndx, el_rowIndex [locIndx]));
    }
#ifndef NDEBUG
//    //Printing the contribution block before dgemm
//    // Valgrind complains if activated, but can be helpful sometimes
//    p = 0;
//    PRLEVEL (p, ("\n%%Before DGEMM:"));
//    if (p <= 0)
//        paru_print_element (paruMatInfo, eli);
#endif


    double *el_numbers = numeric_pointer (curEl);
    paru_dgemm(pivotalFront, uPart, el_numbers, fp, rowCount, colCount);

#ifndef NDEBUG
    //Printing the contribution block after dgemm
    p = 1;
    PRLEVEL (p, ("\n%%After DGEMM:"));
    if (p <= 0)
        paru_print_element (paruMatInfo, eli);
#endif


    /**** 7 **** Count number of rows and columsn of prior CBs to asslemble ***/ 

    paruMatInfo->time_stamp[f]++; //invalidating all the marks
    PRLEVEL (-1, ("\n%%||||  Start FourPass %ld ||||\n", f));
    paru_fourPass (paruMatInfo, f, fp, cc);
    PRLEVEL (-1, ("\n%%||||  Finish FourPass %ld ||||\n", f));


    ////////////////////////////////////////////////////////////////////////////

    // adding tuples for current front
    for (Int i = 0; i < colCount; ++ i) {
        Tuple colTuple;
        colTuple.e = eli;
        colTuple.f = i;
        if (paru_add_colTuple (ColList, CBColList[i], colTuple, cc) ){
            printf("%% Out of memory: add_colTuple \n");
            return;
        }
    }
    for (Int i = fp; i < rowCount; ++ i) {
        Int locIndx = i-fp; 
       Tuple rowTuple;
        rowTuple.e = eli;
        rowTuple.f = locIndx;
        if (paru_add_rowTuple (RowList, fsRowList[i], rowTuple, cc) ){
            printf("%% Out of memory: add_colTuple \n");
            return; 
        }

    }

#ifndef NDEBUG
    //Printing the contribution block after prior blocks assembly
    p = 1;
    PRLEVEL (p, ("\n%%After prior blocks assembly:"));
    if (p <= 0)
        paru_print_element (paruMatInfo, eli);
#endif



    Work->rowMark +=  rowCount;
    rowMark = Work -> rowMark;

    Work->colMark += colCount;
    colMark = Work -> colMark;
    /* Trying to DEBUG */ 
    Work -> elCMark += colCount;
    Work -> elRMark += colCount;


#ifndef NDEBUG /* chekcing isRowInFront to be zero */
    for (Int i = 0; i < m; i++)
        ASSERT ( isRowInFront [i] < rowMark);
#endif


    PRLEVEL (1, ("%%rowCount =%ld\n", rowCount));
    PRLEVEL (1, ("%%colCount =%ld\n", colCount));
    PRLEVEL (-1, ("fp =%ld;\n", fp));
    PRLEVEL (1, ("%%~~~~~~~Assemble Front %ld finished\n", f));
}
