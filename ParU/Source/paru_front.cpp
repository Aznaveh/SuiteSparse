/** =========================================================================  /
 * =======================  paru_front   ====================================  /
 * ========================================================================== */

/*! @brief Computing factorization of current front and doing the numerical
 * assembly that ancestors will assemble. Degree update will be used in this
 * version. Just like ./paru_assemble.cpp
 *  
 *
 * @param  a list of tuples and the the tuple we want to add
 * @return 0 on sucess 
 *
 *  @author Aznaveh
 */
#include "Parallel_LU.hpp"
int paru_front ( paru_matrix *paruMatInfo, 
        /* RowCol list/tuples and LUsym handle */
        Int f, /* front need to be assembled */
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
    Int m = paruMatInfo-> m;
    Int *Super = LUsym->Super;
    /* ---------------------------------------------------------------------- */
    /* get the front F  */
    /* ---------------------------------------------------------------------- */


    PRLEVEL (-1, ("%%~~~~~~~  Assemble Front %ld start ~~~~\n", f));
    /* pivotal columns Super [f] ... Super [f+1]-1 */
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];
    Int fp = col2 - col1;   /* first fp columns are pivotal */ 


    paru_Element **elementList = paruMatInfo->elementList;
    work_struct *Work =  paruMatInfo->Work;

    // Int elRMark = Work -> elRMark;
    // Int elCMark = Work -> elCMark;

    Int *elCol = Work -> elCol;

    PRLEVEL (1, ("%% fp=%ld pivotal columns:clo1=%ld...col2=%ld\n", 
                fp, col1, col2-1));
    ASSERT (fp > 0 );

    /* computing number of rows, set union */
    tupleList *ColList = paruMatInfo->ColList;

    Int panel_width = paruMatInfo->panel_width;
    Int num_panels = (Int) ceil( (double)fp/panel_width);
    // panel_row shows number of rows in each panel. Needs to be initialized in
    // my new algorithm
    Int *panel_row = (Int*) paru_calloc ( num_panels , sizeof (Int), cc);
    if (panel_row == NULL )
    {
        printf ("%% Out of memory when tried to allocate for panel_row%ld",f);
        return 1;
    }

    Int *snM = LUsym->super2atree;
    Int eli = snM [f]; 

    Int *isRowInFront = Work->rowSize; 
    Int *rowMarkp = Work->rowMark;
    Int rowMark = rowMarkp[eli];


    Int fm = LUsym->Fm[f];     /* Upper bound number of rows of F */ 
    PRLEVEL (1, ("%% the size of fm is %ld\n",fm));
    Int *frowList = (Int*) paru_alloc (fm, sizeof (Int), cc);
    if (frowList == NULL )
    {
        printf ("%% Out of memory when tried to allocate for frowList %ld",f);
        paru_free ( num_panels, sizeof (Int), panel_row, cc);
        return 1;
    }
    paruMatInfo->frowList[f] = frowList;

    std::set<Int>::iterator it;
#ifndef NDEBUG
    std::set<Int> stl_rowSet;
#endif 

    // Initializing relative index validation flag of current front
    paru_init_rel (paruMatInfo, f);
    Int time_f = paruMatInfo->time_stamp[f];

    PRLEVEL (0, ("%% Begin of Front %ld time_f = %ld\n", f, time_f));

    //Int panel_num = 0; 
    paruMatInfo->frowCount[f] = 0;
    //Int rowCount = 0;


    /************ Making the heap from list of the immediate children ******/
 //   PRLEVEL (1, ("%% Next: work on the heap \n"));
 //   paru_make_heap(f, paruMatInfo);
 //   PRLEVEL (1, ("%% Done: work on the heap \n"));

    /********************** pivotal column assembly  **************************/
    /***************  assembling the pivotal part of the front ****************/
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

    std::vector<Int> pivotal_elements;
    PRLEVEL (1, ("%% Next: work on pivotal column assembly\n"));
    paru_pivotal (paruMatInfo, pivotal_elements, panel_row , f, cc);
    PRLEVEL (1, ("%% Done: work on pivotal column assembly\n"));

    Int rowCount = paruMatInfo->frowCount[f];
    frowList = paruMatInfo->frowList[f];


#ifndef NDEBUG /* chekcing first part of Work to be zero */
    rowMark = rowMarkp[eli];
    PRLEVEL (1, ("%% rowMark=%ld;\n", rowMark));
    for (Int i = 0; i < m; i++)
    {  
        if ( isRowInFront [i] >= rowMark)
            PRLEVEL (1, ("%%rowMark = %ld, isRowInFront[%ld] = %ld\n", 
                        rowMark ,i,
                        isRowInFront [i]));
        //ASSERT ( isRowInFront [i] < rowMark);
    }
#endif 



    paru_fac *LUs =  paruMatInfo->partial_LUs;
    double *pivotalFront = LUs[f].p;
    LUs[f].m = rowCount;
    LUs[f].n = fp;

    /***************  factorizing the fully summed part of the matrix        ***
     *****  a set of pivot is found in this part that is crucial to assemble **/
    PRLEVEL (1, ("%% rowCount =%ld\n", rowCount));

#ifndef NDEBUG  // Printing the list of rows
    Int p = 1;
    PRLEVEL (p, ("%% Befor factorization (inside assemble): \n"));
    for (Int i = 0; i < rowCount; i++)
        PRLEVEL (p, ("%% frowList [%ld] =%ld\n",i, frowList [i]));
    PRLEVEL (p, ("\n"));

#endif

    Int fn = LUsym->Cm[f];     /* Upper bound number of cols of F */ 
    std::set<Int> stl_colSet;  /* used in this scope */



    if (rowCount < fp)
    {
        PRLEVEL (1, ("%% %ldx%ld \n",rowCount, fp));
        printf ("structural problem\n");
        paru_free ( num_panels, sizeof (Int), panel_row, cc);
        return 1;
    }

    Int start_fac = paruMatInfo->time_stamp[f]; 
    PRLEVEL (1, ("%% start_fac= %ld\n",start_fac));

    Int fac = paru_factorize(pivotalFront, frowList, rowCount, f, start_fac,
            panel_row, stl_colSet, pivotal_elements, paruMatInfo);
    time_f = paruMatInfo->time_stamp[f]; 
    PRLEVEL (1, ("%%After factorization time_f = %ld\n",time_f));

    /* To this point fully summed part of the front is computed and L and U    /  
     *  The next part is to find columns of nonfully summed then rows
     *  the rest of the matrix and doing TRSM and GEMM,                       */

    PRLEVEL (0, ("%% num_panels = %ld\n", num_panels));
    paru_free (num_panels, sizeof (Int), panel_row, cc);
    PRLEVEL (0, ("%% After free num_panels = %ld\n", num_panels));

    if (fac < 0)
    {
        printf ("%% Some problem in factorization \n");
        return 1;
    }

#ifndef NDEBUG  // Printing the list of rows
    p = 1;
    PRLEVEL (p, ("%% After factorization (inside assemble): \n"));
    for (Int i = 0; i < rowCount; i++)
        PRLEVEL (p, ("%% frowList [%ld] =%ld\n",i, frowList [i]));
    PRLEVEL (p, ("\n"));
#endif


#ifndef NDEBUG  // Printing the permutation
    p = 1;
    PRLEVEL (p, ("%% pivotal rows:\n"));
    for (Int i = 0; i < fp; i++)
        PRLEVEL (p, ("%% frowList[%ld] =%ld\n",i, frowList[i]));
    PRLEVEL (p, ("%% =======\n"));
    for (Int i = fp; i < rowCount; i++)
        PRLEVEL (p, ("%% frowList[%ld] =%ld\n",i, frowList[i]));
    PRLEVEL (p, ("\n"));
#endif

#ifndef NDEBUG  // Printing the pivotal front
    p = -1;
    PRLEVEL (p, ("%%L part:\n"));

    //col permutatin
    PRLEVEL (p, ("cols{%ld} = [",f+1));
    for (Int c = col1; c < col2; c++) 
        PRLEVEL (p, ("%ld ", c+1));
    PRLEVEL (p, ("];\n"));

    //row permutatin
    PRLEVEL (p, ("rows{%ld} = [",f+1));
    for (Int r = 0; r < rowCount; r++)
        PRLEVEL (p, ("%ld ", frowList [r]+1)); //Matlab is base 1
    PRLEVEL (p, ("];\n"));

    //inv row permutatin

    PRLEVEL (p, ("Luf{%ld}= [",f+1));
    for (Int r = 0; r < rowCount; r++)
    {
        PRLEVEL (p, (" "));
        for (Int c = col1; c < col2; c++)
            PRLEVEL (p, (" %.16g ", pivotalFront [(c-col1)*rowCount + r]));
        PRLEVEL (p, (";\n   "));
    }
    PRLEVEL (p, ("];\n"));
    //just in cases that there is no U for MATLAB
    PRLEVEL (p, ("Us{%ld} =[];\n", f+1));
    PRLEVEL (p, ("Ucols{%ld}=[];\n",f+1));
    PRLEVEL (p, ("Urows{%ld}=[];\n",f+1));
#endif


    Int colCount = stl_colSet.size();
    ASSERT ( fn >= colCount );

    Int *fcolList = NULL;

    if (fn != 0) 
    {
        PRLEVEL (1, ("%% fp=%ld fn=%ld \n", fp, fn));
        fcolList = (Int*) paru_alloc (stl_colSet.size(), sizeof (Int), cc);

        if (fcolList == NULL)
        {
            printf ("%% Out of memory when tried to allocate for fcolList=%ld" 
                    "with the size %ld", f, fn);
            return 1;
        }
    }

    paruMatInfo->fcolList[f] = fcolList;


    std::vector<Int>** heapList = paruMatInfo->heapList;
    std::vector<Int>* curHeap = heapList[eli];

    // EXIT point HERE 
    if (colCount == 0 )
    {  // there is no CB, Nothing to be done
        //Work->rowMark +=  rowCount;
        paruMatInfo->fcolCount[f] = 0;
        PRLEVEL (1, ("%%Heap freed inside front %p id=%ld\n",curHeap, eli ));
        delete curHeap;
        paruMatInfo->heapList[eli] = nullptr;
        PRLEVEL (1, ("%% pivotalFront =%p\n",pivotalFront));
        return 0;
    }

    //fcolList copy from the stl_colSet
    //hasing from fcolList indices to column index 
    std::unordered_map <Int, Int> colHash; 
    Int i = 0;
    for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
    {
        colHash.insert({*it , i});
        fcolList[i++] = *it;
    }


    /**** 5 ** assemble U part         Row by Row                          ****/ 

    double *uPart = 
        (double*) paru_calloc (fp*colCount, sizeof (double), cc);
    if ( uPart == NULL )
    {
        printf ("%% Out of memory when tried to allocate for U part %ld",f);
        return 1;
    }

    paru_fac *Us =  paruMatInfo->partial_Us;
    Us[f].m = fp;
    Us[f].n = colCount;
    paruMatInfo->fcolCount[f] = colCount;
    ASSERT (Us[f].p == NULL);
    Us[f].p = uPart;


    tupleList *RowList = paruMatInfo->RowList;
    for (Int i = 0; i < fp; i++)
    {
        Int curFsRowIndex = i; //current fully summed row index
        Int curFsRow = frowList [curFsRowIndex];
        PRLEVEL (1, ("%% curFsRow =%ld\n", curFsRow));
        tupleList *curRowTupleList = &RowList [curFsRow];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        ASSERT (numTuple <= m);
        paru_Tuple *listRowTuples = curRowTupleList->list;
        PRLEVEL (1, ("%% numTuple = %ld\n", numTuple));
        for (Int k = 0; k < numTuple; k++)
        {
            paru_Tuple curTpl = listRowTuples [k];
            Int e = curTpl.e;
            paru_Element *el = elementList[e];
            if (el == NULL) continue;

            Int curRowIndex = curTpl.f;
            Int nEl = el->ncols;
            //Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_rowIndex = (Int*) (el+1) + nEl; 

            if(el_rowIndex[curRowIndex] < 0 ) continue; 

            Int mEl = el->nrows;
            //Int *rowRelIndex = relRowInd (el);
            Int *rowRelIndex = (Int*) (el+1) + 2*nEl + mEl; 

            //Int *colRelIndex = relColInd (el);
            Int *colRelIndex =  (Int*)(el+1)+ nEl + mEl;
            Int *colIndex = (Int*)(el+1);


            PRLEVEL (1, ("%% curFsRowIndex =%ld\n", curFsRowIndex));
            ASSERT (el_rowIndex[curRowIndex] == curFsRow);
            ASSERT (curRowIndex < mEl);
            PRLEVEL (1, ("%% curColIndex =%ld\n", curRowIndex));

            //double *el_Num = numeric_pointer (el);
            double *el_Num =  (double*)((Int*) (el+1) + 2*nEl + 2*mEl); 
            PRLEVEL (1, ("%% element= %ld  nEl =%ld \n",e, nEl));


            assemble_row_hash (el_Num, uPart, mEl, nEl, fp, 
                    curRowIndex, curFsRowIndex, colIndex, colHash);


            //FLIP(el_rowIndex[curRowIndex]); //marking row assembled
            el_rowIndex[curRowIndex] = -1;
            rowRelIndex [curRowIndex] = -1;
            el->nrowsleft--;  
            if (el->nrowsleft == 0) 
            { //free el
                Int tot_size = sizeof(paru_Element) +
                    sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
                PRLEVEL (-1, ("%%inside Front: Free %ld\n",e));
                paru_free (1, tot_size, el, cc);
                elementList[e] = NULL;
            }
        }
    }

#ifndef NDEBUG  // Printing the  U part
    p = 0;
    PRLEVEL (p, ("%% U part Before TRSM: %ld x %ld\n", fp, colCount));
    PRLEVEL (p, ("%% U\t"));
    for (Int i = 0; i < colCount; i++)
        PRLEVEL (p, ("%ld\t\t", fcolList[i]));
    PRLEVEL (p, ("\n"));
    for (Int i = 0; i < fp; i++)
    {
        PRLEVEL (p, ("%% %ld\t",  frowList [i]));
        for (Int j = 0; j < colCount; j++)
            PRLEVEL (p, (" %2.5lf\t", uPart[j*fp+i]));
        PRLEVEL (p, ("\n"));
    }

#endif


    /**** 6 ****                     TRSM and DGEMM                         ***/ 

    paru_trsm(pivotalFront , uPart, fp, rowCount, colCount);

#ifdef COUNT_FLOPS
    paruMatInfo->flp_cnt_trsm += (double) (fp+1)*fp*colCount;
#ifndef NDEBUG  
    p = 0;
    PRLEVEL (p, ("\n%% FlopCount Trsm front %ld %ld ",fp, colCount  ));
    PRLEVEL (p, ("cnt = %lf\n ",   paruMatInfo->flp_cnt_trsm ));
#endif
#endif

#ifndef NDEBUG  // Printing the  U part
    p = -1;
    PRLEVEL (p, ("%% rowCount=%ld;\n",rowCount));
    PRLEVEL (p, ("%% U part After TRSM: %ld x %ld\n", fp, colCount));

    PRLEVEL (p, ("Ucols{%ld} = [",f+1));
    for (Int i = 0; i < colCount; i++)
        PRLEVEL (p, ("%ld ", fcolList[i]+1));
    PRLEVEL (p, ("];\n"));


    PRLEVEL (p, ("Urows{%ld} = [",f+1));
    for (Int i = 0; i < fp; i++)
        PRLEVEL (p, ("%ld ",  frowList [i]+1));
    PRLEVEL (p, ("];\n"));


    PRLEVEL (p, ("Us{%ld} = [",f+1));

    for (Int i = 0; i < fp; i++)
    {
        for (Int j = 0; j < colCount; j++)
            PRLEVEL (p, (" %.16g ", uPart[j*fp+i]));
        PRLEVEL (p, (";\n    "));
    }
    PRLEVEL (p, ("];\n"));
#endif


    paru_Element *curEl;
    PRLEVEL (1, ("%% rowCount=%ld, colCount=%ld, fp=%ld\n",
                rowCount, colCount, fp));
    PRLEVEL (1, ("%% curEl is %ld by %ld\n",rowCount-fp,colCount));
    if (fp < rowCount )
    { 
        curEl = elementList[eli] = paru_create_element (rowCount-fp,
                colCount, 0 ,cc); // allocating an un-initialized part of memory

        // While insided the DGEMM BETA == 0
        if ( curEl == NULL )
        {
            printf ("%% Out of memory when tried to allocate current CB %ld",
                    eli);
            return 1;
        }
        PRLEVEL (1, ("%% Created ele %ld in curEl =%p\n", eli, curEl));
    }
    else //EXIT point
    {   //NO rows for current contribution block
        delete curHeap;
        paruMatInfo->heapList[eli] = nullptr;
        PRLEVEL (1, ("%%(2)Heap freed inside front %p id=%ld\n",
                    curHeap, eli ));
        PRLEVEL (1, ("%% pivotalFront =%p\n",pivotalFront));
        return 0;
    }

    // Initializing curEl global indices
    //Int *el_colIndex = colIndex_pointer (curEl);
    Int *el_colIndex = (Int*)(curEl+1);
    curEl->lac = 0;
    Int *lacList = paruMatInfo->lacList;
    lacList[eli] = fcolList[0];
    for (Int i = 0; i < colCount; ++ i) 
        el_colIndex [i] = fcolList[i];
    Int *el_rowIndex = rowIndex_pointer (curEl);
    for (Int i = fp; i < rowCount; ++ i) 
    {
        Int locIndx = i-fp; 
        Int curRow = frowList [i];
        el_rowIndex [locIndx] = curRow;
        //Updating isRowInFront after the pivoting
        //attention it is the old rowMark not the updated rowMarkp + eli
        //If I decide to add rowMark I should change paru_pivotal
        //TODO decide about adding rowMark here later // + rowMark ; 
        isRowInFront [curRow] = locIndx ;
        PRLEVEL (1, ("%% el_rowIndex [%ld] =%ld\n",
                    locIndx, el_rowIndex [locIndx]));
    }

    //double *el_numbers = numeric_pointer (curEl);
    double *el_numbers = (double*)
        ((Int*)(curEl+1) + 2*colCount + 2*(rowCount-fp));

    paru_dgemm(pivotalFront, uPart, el_numbers, fp, rowCount, colCount);
    //printf ("%ld %ld %ld \n",rowCount-fp, colCount, fp);

#ifdef COUNT_FLOPS
    paruMatInfo->flp_cnt_dgemm += (double) 2*(rowCount -fp)*fp*colCount;

#ifndef NDEBUG  
    PRLEVEL (p, ("\n%% FlopCount Dgemm front %ld %ld %ld \n",rowCount -fp,
                fp, colCount  ));
    PRLEVEL (p, ("%ld %ld %ld \n",rowCount -fp,
                fp, colCount  ));

    PRLEVEL (p, ("cnt = %lf\n ",   paruMatInfo->flp_cnt_dgemm ));
#endif
#endif

#ifndef NDEBUG
    //Printing the contribution block after dgemm
    p = 1;
    PRLEVEL (p, ("\n%%After DGEMM:"));
    if (p <= 0)
        paru_print_element (paruMatInfo, eli);
#endif


    /**** 7 **** Count number of rows and columsn of prior CBs to asslemble ***/ 

    paruMatInfo->time_stamp[f]++; //invalidating all the marks
    PRLEVEL (-1, ("\n%%||||  Start Finalize %ld ||||\n", f));
    //paru_finalize (paruMatInfo,  f, start_fac, cc);
    paru_prior_assemble ( f, start_fac, 
            pivotal_elements, colHash, paruMatInfo, cc);
    PRLEVEL (-1, ("\n%%||||  Finish Finalize %ld ||||\n", f));


    ////////////////////////////////////////////////////////////////////////////

    // adding tuples for current front
    for (Int i = 0; i < colCount; ++ i) 
    {
        paru_Tuple colTuple;
        colTuple.e = eli;
        colTuple.f = i;
        if (paru_add_colTuple (ColList, fcolList[i], colTuple, cc) )
        {
            printf("%% Out of memory: add_colTuple \n");
            return 1;
        }
    }
    for (Int i = fp; i < rowCount; ++ i) 
    {
        Int locIndx = i-fp; 
        paru_Tuple rowTuple;
        rowTuple.e = eli;
        rowTuple.f = locIndx;
        if (paru_add_rowTuple (RowList, frowList[i], rowTuple, cc) )
        {
            printf("%% Out of memory: add_colTuple \n");
            return 1; 
        }

    }


#ifndef NDEBUG /* chekcing if isRowInFront is correct */
    rowMark = rowMarkp[eli];
    Int *Sleft = LUsym->Sleft;
    for (Int i = Sleft[col1]; i < Sleft[Super[f+1]]; i++)
        ASSERT ( isRowInFront [i] < rowMark);
#endif


    PRLEVEL (1, ("%%rowCount =%ld\n", rowCount));
    PRLEVEL (1, ("%%colCount =%ld\n", colCount));
    PRLEVEL (-1, ("fp =%ld;\n", fp));
    PRLEVEL (1, ("%%~~~~~~~Assemble Front %ld finished\n", f));
    return 0;
}
