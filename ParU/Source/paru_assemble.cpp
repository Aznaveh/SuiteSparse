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

    PRLEVEL (1, ("fp=%ld pivotal columns:clo1=%ld...col2=%ld\n", 
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

    Int listP= 0;
    work_struct *Work =  paruMatInfo->Work;
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
    for (Int c = col1; c < col2; c++){
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listColTuples = curColTupleList->list;
        PRLEVEL (1, ("c =%ld  numTuple = %ld\n", c, numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Element **elementList = paruMatInfo->elementList;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;
            Int *el_colrowIndex = (Int*)(curEl+1);  // pointers to element index 
            Int *el_rowIndex = el_colrowIndex + nEl;// pointers to row indices
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));
            /*! TODO: check if any row/col nulified     */
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
                    PRLEVEL (1, ("curRow =%ld listP=%ld\n", curRow, listP));
                    fsRowList [listP] = curRow;
                    PRLEVEL (1, ("listP=%ld FLIP(listP)=%ld\n", 
                                listP, FLIP (listP) ));
                    isRowInFront [curRow] = rowMark + listP++; 
               }
                ASSERT (listP <= m); 

#if 0
            Int rS;
                for (rS = 0; rS < rowsP; rS++){
                    PRLEVEL (1, ("rS =%ld rEl=%ld\n", rS, rEl));
                    if (curRow == rowSet [rS])
                        break; 
                }
                if ( rS == rowsP){ // count the new row
                    if (rowsP >= setSize){
                        PRLEVEL (1, ("rowsP =%ld setSize=%ld\n", 
                                    rowsP, setSize));
                        PRLEVEL (1, ("setSize*2 =%ld\n", setSize*2));
                        Int *newSet= (Int*) paru_realloc (setSize*2, 
                                sizeof (Int), rowSet, &setSize, cc);

                        if(newSet== NULL){
                            printf("Error in allocating memory for rows\n");
                            return;
                        }
                        rowSet = newSet;
                    }
                    rowSet [rowsP++] = curRow;
                }
#endif                
            }
        }
    }

#if 0    
    //shrinking allocated space
    rowSet= (Int*) paru_realloc (rowsP, sizeof (Int), rowSet, &setSize, cc);
#endif    

#ifndef NDEBUG /* Checking if pivotal rows are correct */
    Int p = 1;
    PRLEVEL (p, ("There are %ld rows in this front: \n", listP));
    for (Int i = 0; i < listP; i++)
        PRLEVEL (p, (" %ld", fsRowList [i]));
    PRLEVEL (p, ("\n"));
    Int stl_size = stl_rowSet.size();
    if (listP != stl_size){
        PRLEVEL (1, ("#######################\n"));
        PRLEVEL (1, ("STL %ld:\n",stl_size));
        for (it = stl_rowSet.begin(); it != stl_rowSet.end(); it++)
            PRLEVEL (1, (" %ld", *it));
        PRLEVEL (1, ("\nMy Set %ld:\n",listP));
        for (Int i = 0; i < listP; i++)
            PRLEVEL (1, (" %ld", fsRowList [i]));
        PRLEVEL (1, ("\n"));
    }
    ASSERT (listP == stl_size );
 //   ASSERT (listP >= fp ); // otherwise it is a singular matrix

#endif 



    double *pivotalFront = 
        (double*) paru_calloc (listP*fp, sizeof (double), cc);
    /* assembling the pivotal part of the front */
    /* 
     *                  
     *  curEl           nEl           
     *              6, 7, 3, 10
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
     *               listP      1    2 | X  Y  .  .  .
     *                          2    4 | *  *  .  .  .  isRowInFront[4] == 2
     *                          3   17 | X  Y  .  .  . 
     * */

    for (Int c = col1; c < col2; c++){
        tupleList *curTupleList = &ColList[c];
        Int numTuple = curTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *l = curTupleList->list;
        PRLEVEL (1, ("c =%ld numTuple = %ld\n", c, numTuple));

        for (Int i = 0; i < numTuple; i++){

            Tuple curTpl = l [i];
            Int e = curTpl.e;
            PRLEVEL (1, ("col=%ld, (%ld,%ld)\n", c, e, f));
            FLIP (curTpl.e); //Nullifying tuple
            curTupleList->numTuple--;
            
            Element **elementList = paruMatInfo->elementList;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;
            Int *el_colIndex = (Int*)(curEl+1);    // pointers to element index
            Int *el_rowIndex = el_colIndex + nEl;  // pointers to row indices

            Int curColIndex = curTpl.f;
            ASSERT (el_colIndex[curColIndex] == c);
            FLIP(el_colIndex[curColIndex]); //Nullifying the column
            curEl->ncolsleft--;
            PRLEVEL (1, ("curColIndex =%ld\n", curColIndex));

            ASSERT (curColIndex < nEl);
            double *el_colrowNum = (double*)(el_colIndex + mEl + nEl); 
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));

            for (Int rEl = 0; rEl < mEl; rEl++){   
                Int curRow = el_rowIndex [rEl]; 
                PRLEVEL (1, ("curRow =%ld\n", curRow));
                ASSERT (curRow < m ) ;
                ASSERT (isRowInFront [curRow] != -1);
                Int rowIndexF = isRowInFront [curRow] - rowMark;
                Int colIndexF = c - col1;
                PRLEVEL (1, ("rowIndexF = %ld\n", rowIndexF));
                PRLEVEL (1, (" colIndexF*listP + rowIndexF=%ld\n",
                            colIndexF*listP + rowIndexF));
                ASSERT ( colIndexF*listP + rowIndexF < listP * fp);
                ASSERT ( curColIndex*mEl + rEl < mEl*nEl);
                pivotalFront [colIndexF*listP + rowIndexF] += 
                    el_colrowNum [ curColIndex*mEl + rEl];
            }
        }
    }
    
#ifndef NDEBUG  // Printing the pivotal front
    p = 0;
    PRLEVEL (p, ("x\t"));
    for (Int c = col1; c < col2; c++) {
        PRLEVEL (p, ("%ld\t", c));
    }
    PRLEVEL (p, ("\n"));
    for (Int r = 0; r < listP; r++){
        PRLEVEL (p, ("%ld\t", fsRowList [r]));
        for (Int c = col1; c < col2; c++){
            PRLEVEL (p, (" %3.1lf\t", pivotalFront [(c-col1)*listP + r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif

    /*     factorizing the fully summed part of the matrix                     /
     *     a set of pivot is found in this part that is crucial to assemble   */
    PRLEVEL (1, ("listP =%ld\n", listP));
    int *ipiv =(int *) (Work->scratch+listP+1); /* using the rest of scratch for 
                                               permutation; Not sure about 1  */
    paru_factorize (pivotalFront, listP, fp, ipiv, cc );

    /* To this point fully summed part of the front is computed and L and U    /  
     *  The next part is to find columns of nonfully summed then rows
     *  the rest of the matrix and doing TRSM and GEMM,                       */
    if (ipiv [0] < 0){
        printf ("Singular Matrix\n");
        return;
    }
#ifndef NDEBUG  // Printing the permutation
    p = 1;
    PRLEVEL (p, ("permutation:\n"));
    for (int i = 0; i < fp; i++){
        PRLEVEL (p, ("ipiv[%d] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (p, ("\n"));
#endif
    /* Set union for pivotal columns */
    Int *CBcolSet = Work -> colSize;
    Int colMark = Work -> colMark;
    if (colMark < 0) {  // in rare case of overflow
        memset (isRowInFront, -1, n*sizeof(Int));
        colMark = Work-> colMark = 0;
    }

#ifndef NDEBUG
    std::set<Int> stl_colSet;
#endif 
    
    tupleList *RowList = paruMatInfo->RowList;
    for (Int i = 0; i < fp; i++){
        Int curFsRowIndex =(Int) ipiv [i]; //current fully summed row index
        PRLEVEL (0, ("curFsRowIndex = %ld\n", curFsRowIndex));
        ASSERT (curFsRowIndex < m);
        Int curFsRow = fsRowList [curFsRowIndex];
        PRLEVEL (0, ("curFsRow =%ld\n", curFsRow));
        tupleList *curRowTupleList = &RowList [curFsRowIndex];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *listRowTuples = curRowTupleList->list;
        PRLEVEL (1, ("numTuple = %ld\n", numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Element **elementList = paruMatInfo->elementList;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int nEl = curEl->ncols;
            Int *el_colIndex = (Int*)(curEl+1);  // pointers to element index 
            Int *el_rowIndex = el_colrowIndex + nEl;// pointers to row indices
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));
            for (Int cEl = 0; cEl < nEl; cEl++){
                Int curCol = el_colIndex [cEl]; 
                PRLEVEL (1, ("curRow =%ld\n", curRow));
    /*! TODO: for each CBrow in ipiv 0:fp add all columns to the set     */
#ifndef NDEBUG
                stl_colSet.insert (curCol);
#endif
 
            }
        }


    }




#ifdef NotUsingMark
    /*Not used now, I am using rowMark to avoid this*/
    /* setting W for next iteration
       for (Int i = 0; i < listP; i++){
       Int curRow = fsRowList [i];
       ASSERT (curRow < m );
       ASSERT (isRowInFront [curRow] != -1);
       isRowInFront  [curRow] = -1;
       } */
#endif

    Work->rowMark += listP;
    rowMark = Work->rowMark;


#ifndef NDEBUG /* chekcing isRowInFront to be zero */
    for (Int i = 0; i < m; i++)
        ASSERT ( isRowInFront [i] < rowMark);
#endif


    paru_free (listP*fp, sizeof (Int), pivotalFront, cc);

#if 0
    paru_free (setSize, sizeof (Int), rowSet, cc);
#endif
}
