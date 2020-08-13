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

    DEBUGLEVEL(-2);
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
    Int *panel_row = (Int*) paru_alloc ( (Int) ceil( (double)fp/panel_width) , 
            sizeof (Int), cc);
    if (panel_row == NULL )
    {
        printf ("%% Out of memory when tried to allocate for panel_row%ld",f);
        return 1;
    }
 
    Int *isRowInFront = Work->rowSize; 
    Int rowMark = Work->rowMark;
    PRLEVEL (-1, ("rowMark=%ld;\n", rowMark));

    if (rowMark < 0) 
    {  // in rare case of overflow
        memset (isRowInFront, -1, m*sizeof(Int));
        rowMark = Work->rowMark = 0;
    }


    Int fm = LUsym->Fm[f];     /* Upper bound number of rows of F */ 
    PRLEVEL (1, ("%% the size of fm is %ld\n",fm));
    Int *frowList = (Int*) paru_alloc (fm, sizeof (Int), cc);
    if (frowList == NULL )
    {
        printf ("%% Out of memory when tried to allocate for frowList %ld",f);
        paru_free ( (Int) ceil( (double)fp/panel_width) ,
                sizeof (Int), panel_row, cc);
        return 1;
    }
    paruMatInfo->frowList[f] = frowList;

#ifndef NDEBUG /* chekcing first part of Work to be zero */
    for (Int i = 0; i < m; i++)
    {  
        if ( isRowInFront [i] >= rowMark)
            PRLEVEL (1, ("%%rowMark = %ld, isRowInFront[%ld] = %ld\n", 
                        rowMark ,i,
                        isRowInFront [i]));
        ASSERT ( isRowInFront [i] < rowMark);
    }
#endif 


    std::set<Int>::iterator it;
#ifndef NDEBUG
    std::set<Int> stl_rowSet;
#endif 
    // Initializing relative index validation flag of current front
    paru_init_rel (paruMatInfo, f);
    Int time_f = paruMatInfo->time_stamp[f];

    PRLEVEL (0, ("%% Begin of Front %ld time_f = %ld\n", f, time_f));

    Int panel_num = 0; 
    paruMatInfo->frowCount[f] = 0;
    Int rowCount= 0;


    /*************** Making the link list of the immediate children ***********/
    paru_make_heap(paruMatInfo, f);
#ifndef NDEBUG  
    Int p = 1;
#endif
    Int * aChild = LUsym->aChild;
    Int * aChildp = LUsym->aChildp;
    Int *snM = LUsym->super2atree;

//    curEl->next = aChild[aChildp[eli]];      //first immediate child
//    curEl->prev= -1;
//    curEl->lad= aChild[aChildp[eli+1]-1];        //last immediate child
//    PRLEVEL (p, ("%% prev is %ld next is %ld\n", curEl->prev, curEl->next));

    Int eli = snM [f]; 
    Int prEl = eli;
    
//    paru_Element *chel;
//    paru_Element *ladel;
//    Int elidLadch;
//    PRLEVEL (p, ("%% Making the list for front %ld with id %ld\n",f,eli));

//TODO  add eli to th list wait until everything is done
    Int i = 0;
    for (i = aChildp[eli]; i <= aChildp[eli+1]-1; i++) 
    {
//        PRLEVEL (p, ("%% HEERREE \n"));
//        PRLEVEL (p, ("%% i = %ld and aChild[i]=%ld\n", i, aChild[i]));
//        Int elidch = aChild[i];  // element id of the child
//        PRLEVEL (p, ("%% elidch is %ld\n", elidch));
//        chel = elementList [elidch]; // immediate child element
//        ASSERT(chel != NULL);
//        elidLadch = chel->lad;// element if of last active descendent child
//        PRLEVEL (p, ("%% elidLadch is %ld\n", elidLadch));
//        //last active descendent of curent child
//        ladel = elementList[elidLadch]; 
//        ASSERT(ladel!= NULL);
//        ladel->next= aChild[i+1];
//        chel->prev = prEl;
//        PRLEVEL (p, ("%% prev is %ld lad is %ld\n", chel->prev, chel->lad));
//        prEl = elidch; 
//        
//
//        TODO: add all the active elements of aChild[i] to the current list
//        ?check to see if it is still active
    }
//    // last child it could be done inside the loop though it might be faster a
//    // little bit like this
//    chel = elementList [aChild[i+1]-1]; // immediate child element
//    ladel = elementList[elidLadch]; 
//    chel->prev = prEl;
//    ladel->next = -1;
//    PRLEVEL (p, ("%% prev is %ld lad is %ld\n", chel->prev, chel->lad));
//
#ifndef NDEBUG  
    p = 1;
#endif



    /**** 1 ******** finding set of rows in current front *********************/
    for (Int c = col1; c < col2; c++)
    {

        //saving number of rows of last panel if it passed the last col
        if( (c-col1) % panel_width == 0 && c != col1  )
        { 
            panel_row [panel_num++] = rowCount;
        }

        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listColTuples = curColTupleList->list;

#ifndef NDEBUG            
        Int p = 1;
        PRLEVEL (p, ("%%c =%ld  numTuple = %ld\n", c, numTuple));
        if (p <= 0 )
            paru_print_tupleList (ColList, c);
#endif
        for (Int i = 0; i < numTuple; i++)
        {
            paru_Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%%e =%ld  curColIndex = %ld\n", e, curColIndex));

            /*! TODO: Never negate e or f for now... keep it?	 */
            if(e < 0 || curColIndex < 0 ) continue;  //already deleted

            paru_Element *el = elementList[e];
            if (el == NULL) continue;

            //Int *el_colIndex = colIndex_pointer (el);
            Int *el_colIndex = (Int*)(el+1);

            PRLEVEL (1, ("%%point to col = %ld\n", el_colIndex[curColIndex]));
            /*! TODO: Keep the tuple or delete it?	 */
            if (el_colIndex [curColIndex]< 0 ) continue; // already assembled
            ASSERT (el_colIndex[curColIndex] == c);

            PRLEVEL (1, ("%%1:element= %ld elCol=%ld elCMark=%ld \n",
                        e, elCol[e], Work -> elCMark));

            if(el->rValid !=  time_f)
            {  // an element never seen before
                el->rValid = time_f;
#ifndef NDEBUG            
                if (el->rValid >  time_f )
                    PRLEVEL (1, ("%%time_f =%ld  rVal= %ld\n",
                                time_f , el->rValid));
#endif               
                ASSERT(el->rValid <= time_f);
            }
            else 
            { 
                elCol [e]--;    //keep track of number of cols
                PRLEVEL (1, ("%%  element= %ld is seen before \n",e));
                //counting prior element's columns
                continue;       // already have the set of rows
            }
            Int mEl = el->nrows;
            Int nEl = el->ncols;

            //Int *el_rowIndex = rowIndex_pointer (el); 
            Int *el_rowIndex = (Int*)(el+1)+nEl; 

            //Int *rowRelIndex = relRowInd (el);
            Int *rowRelIndex = (Int*)(el+1)+ 2*nEl + mEl;

            //Int *colRelIndex = relColInd (el);
            Int *colRelIndex =  (Int*)(el+1)+ nEl + mEl;


            colRelIndex [curTpl.f] = c - col1; //Initialzing relative index


            PRLEVEL (1, ("%%2:element= %ld  mEl =%ld \n",e, mEl));
            //Set union of all the rows of the current element
            for (Int rEl = 0; rEl < mEl; rEl++)
            {
                Int curRow = el_rowIndex [rEl]; 
                PRLEVEL (1, ("%%@@curRow =%ld rEl=%ld\n", curRow, rEl));
                if (curRow < 0 ) continue; // that row has already deleted
                PRLEVEL (1, ("%%**curRow =%ld rEl=%ld\n", curRow, rEl));
                ASSERT (curRow < m ) ;
#ifndef NDEBUG
                Int p = 1;
                if (p <= 0) paru_print_element (paruMatInfo, e);
                stl_rowSet.insert (curRow);
                PRLEVEL (1, ("%% %p ---> isRowInFront [%ld]=%ld\n", 
                           isRowInFront+curRow, curRow, isRowInFront[curRow]));
#endif

                if (isRowInFront[curRow] < rowMark )
                {  
                    // Adding curRow to the set
                    PRLEVEL (1, ("%%curRow =%ld rowCount=%ld\n", 
                                curRow, rowCount));
                    frowList [rowCount] = curRow;
                    rowRelIndex [rEl] = rowCount ;
                    PRLEVEL (1, ("%%1st: rowRelIndex[%ld] = %ld\n",
                                rEl, rowCount ));
                    isRowInFront [curRow] = rowMark + rowCount++; 
                }
                else
                {
                    rowRelIndex [rEl] = isRowInFront [curRow] - rowMark;
                    PRLEVEL (1, ("%%N1st: rowRelIndex[%ld] = %ld\n",
                                rEl, rowCount ));
                }
                ASSERT (rowCount <= m); 
            }
        }
    }


    panel_row [panel_num++] = rowCount;
    PRLEVEL (1, ("%%col1 =%ld col2 = %ld  panel_num= %ld\n", 
                col1, col2, panel_num));
    ASSERT (panel_num == (Int) ceil( (double)fp/panel_width) );

#ifndef NDEBUG /* Checking if pivotal rows are correct */
    p = 1;
    PRLEVEL (p, ("%%There are %ld rows in this front: \n %%", rowCount));
    for (Int i = 0; i < rowCount; i++)
        PRLEVEL (p, (" %ld", frowList [i]));
    PRLEVEL (p, ("\n"));
    Int stl_rowSize = stl_rowSet.size();
    if (rowCount != stl_rowSize)
    {
        PRLEVEL (p, ("%% STL %ld:\n",stl_rowSize));
        for (it = stl_rowSet.begin(); it != stl_rowSet.end(); it++)
            PRLEVEL (p, ("%% %ld", *it));
        PRLEVEL (p, ("\n%%My Set %ld:\n",rowCount));
        for (Int i = 0; i < rowCount; i++)
            PRLEVEL (p, ("%% %ld", frowList [i]));
        PRLEVEL (p, ("\n"));
    }
    ASSERT (rowCount == stl_rowSize );

#endif 



    double *pivotalFront = 
        (double*) paru_calloc (rowCount*fp, sizeof (double), cc);

    if (pivotalFront == NULL )
    {
        printf ("%% Out of memory when tried to allocate for pivotal part %ld",
                f);
        paru_free ( (Int) ceil( (double)fp/panel_width) ,
                sizeof (Int), panel_row, cc);
        return 1;
    }

    PRLEVEL (1, ("%% pivotalFront =%p \n", pivotalFront));
    PRLEVEL (1, ("%% fm=%ld rowCount=%ld \n", fm, rowCount));
    ASSERT ( fm >= rowCount );
    //freeing extra space for rows
    if (rowCount != fm)
    {
        Int sz = sizeof(Int)*fm; 
        frowList =
            (Int*) paru_realloc (rowCount, sizeof(Int), frowList, &sz, cc);
        paruMatInfo ->frowList[f] = frowList;
    }

    paru_fac *LUs =  paruMatInfo->partial_LUs;
    paruMatInfo->frowCount[f] = rowCount;

    LUs[f].m = rowCount;
    LUs[f].n = fp;
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
    for (Int c = col1; c < col2; c++)
    {
        tupleList *curTupleList = &ColList[c];
        Int numTuple = curTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listColTuples = curTupleList->list;
        PRLEVEL (1, ("%% c =%ld numTuple = %ld\n", c, numTuple));

        Int colIndexF = c - col1;  // relative column index

        for (Int i = 0; i < numTuple; i++)
        {

            paru_Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
            // Assembly of column curColIndex of e in colIndexF
            paru_Element *el = elementList[e];

            if (el == NULL) continue;

            Int curColIndex = curTpl.f;

            //Int *el_colIndex = colIndex_pointer (el);
            Int *el_colIndex = (Int*)(el+1);
            
            if (el_colIndex [curColIndex]< 0 ) continue;

            Int mEl = el->nrows;
            Int nEl = el->ncols;

            //Int *rowRelIndex = relRowInd (el);
            Int *rowRelIndex = (Int*)(el+1) + 2*nEl + mEl;

            //Int *colRelIndex = relColInd (el);
            Int *colRelIndex = (Int*)(el+1) + nEl + mEl;
            ASSERT (el_colIndex[curColIndex] == c);
            PRLEVEL (1, ("%% curColIndex =%ld\n", curColIndex));



            ASSERT (curColIndex < nEl);
            //double *el_Num = numeric_pointer (el);
            double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);

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
            
            //el_colIndex[curColIndex] = -1;
            el_colIndex[curColIndex] = flip (el_colIndex[curColIndex] );

            colRelIndex [curColIndex] = -1;
            el->ncolsleft--;     
            if (el->ncolsleft == 0)
            { //free el
                Int tot_size = sizeof(paru_Element) +
                    sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
                paru_free (1, tot_size, el, cc);
                elementList[e] = NULL;
            }
            else 
            {
                // updating least numbered column to point to an active column
                if (curColIndex <= el->lnc)
                {
                    while(el_colIndex[el->lnc] < 0)
                    {
                        PRLEVEL (1, ("\n%%el->lnc= %ld ",el->lnc));
                        PRLEVEL (1, ("el_colIndex[el->lnc]=%ld :\n",
                                    el_colIndex[el->lnc]));
                        el->lnc++;
                    }
                }
                PRLEVEL (1, ("%%Final: curColIndex=%ld ",curColIndex));
                PRLEVEL (1, ("%%ncols=%ld ",el->ncols));
                PRLEVEL (1, (" el->lnc %ld ",el->lnc));
                PRLEVEL (1, (" el_colIndex[el->lnc]=%ld\n",
                            el_colIndex[el->lnc]));
                ASSERT (el->lnc < el->ncols);
            }



#ifndef NDEBUG  // Printing the pivotal front
            p = 1;
            PRLEVEL (p, ("%%Inside the assembly loop \n"));
            PRLEVEL (p, ("%%x =  \t"));
            for (Int c = col1; c < col2; c++) 
            {
                PRLEVEL (p, ("%% %ld\t\t", c));
            }
            PRLEVEL (p, (" ;\n"));
            for (Int r = 0; r < rowCount; r++)
            {
                PRLEVEL (p, ("%% %ld\t", frowList [r]));
                for (Int c = col1; c < col2; c++)
                {
                    PRLEVEL (p, (" %2.5lf\t", 
                                pivotalFront [(c-col1)*rowCount + r]));
                }
                PRLEVEL (p, ("\n"));
            }
#endif

        }
    }

#ifndef NDEBUG  // Printing the pivotal front
    p = 1;
    PRLEVEL (p, ("%% Before pivoting\n"));
    PRLEVEL (p, ("%% x =  \t"));
    for (Int c = col1; c < col2; c++) 
        PRLEVEL (p, ("%ld\t\t", c));
    PRLEVEL (p, (" ;\n"));
    for (Int r = 0; r < rowCount; r++)
    {
        PRLEVEL (p, ("%% %ld\t", frowList [r]));
        for (Int c = col1; c < col2; c++)
            PRLEVEL (p, (" %2.5lf\t", pivotalFront [(c-col1)*rowCount + r]));
        PRLEVEL (p, ("\n"));
    }
#endif

    /**** 3 ********  factorizing the fully summed part of the matrix        ***
     *****  a set of pivot is found in this part that is crucial to assemble **/
    PRLEVEL (1, ("%% rowCount =%ld\n", rowCount));

#ifndef NDEBUG  // Printing the list of rows
    p = 1;
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
        paru_free ( (Int) ceil( (double)fp/panel_width) ,
                sizeof (Int), panel_row, cc);
        return 1;
    }
    Int next = -1;

    Int start_fac = paruMatInfo->time_stamp[f]; 
    PRLEVEL (1, ("%% start_fac= %ld\n",start_fac));

    Int fac = paru_factorize(pivotalFront, frowList, rowCount, f, panel_row, 
            &next, stl_colSet, paruMatInfo);
    time_f = paruMatInfo->time_stamp[f]; 
    PRLEVEL (1, ("%%After factorization time_f = %ld\n",time_f));

    /* To this point fully summed part of the front is computed and L and U    /  
     *  The next part is to find columns of nonfully summed then rows
     *  the rest of the matrix and doing TRSM and GEMM,                       */

    paru_free ( (Int) ceil( (double)fp/panel_width) ,  // Do not need this space
            sizeof (Int), panel_row, cc);

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


    //hasing from fcolList indices to column index 
    std::unordered_map <Int, Int> colHash; 

    //fcolList copy from the stl_colSet
    i = 0;
    for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
    {
        colHash.insert({*it , i});
        fcolList[i++] = *it;
    }

    std::vector<Int>** heapList = paruMatInfo->heapList;
    std::vector<Int>* curHeap = heapList[eli];
    // EXIT point HERE 
    if (colCount == 0)
    {  // there is no CB, Nothing to be done
        Work->rowMark +=  rowCount;
        paruMatInfo->fcolCount[f] = 0;
        delete curHeap;
        paruMatInfo->heapList[eli] = nullptr;
        PRLEVEL (1, ("%% pivotalFront =%p\n",pivotalFront));
        return 0;
    }


    /**** 5 ** assemble U part         Row by Row                          ****/ 

    //TODO the index is not correct
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

            //TODO: fix here 
            //           assemble_row (el_Num, uPart, mEl, nEl, fp, 
            //                   curRowIndex, curFsRowIndex, colRelIndex);

            assemble_row_hash (el_Num, uPart, mEl, nEl, fp, 
                    curRowIndex, curFsRowIndex, colIndex, colHash);


            //FLIP(el_rowIndex[curRowIndex]); //marking column as assembled
            el_rowIndex[curRowIndex] = -1;
            rowRelIndex [curRowIndex] = -1;
            el->nrowsleft--;  
            if (el->nrowsleft == 0)
            { //free el
                Int tot_size = sizeof(paru_Element) +
                    sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
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
    if (fp <= rowCount )
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
        PRLEVEL (1, ("%% curEl =%p\n", curEl));
    }



    // Initializing curEl global indices
    Int *el_colIndex = colIndex_pointer (curEl);
    for (Int i = 0; i < colCount; ++ i) 
        el_colIndex [i] = fcolList[i];
    Int *el_rowIndex = rowIndex_pointer (curEl);
    for (Int i = fp; i < rowCount; ++ i) 
    {
        Int locIndx = i-fp; 
        el_rowIndex [locIndx] = frowList[i];
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
    paru_finalize (paruMatInfo,  f, start_fac, cc);
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

    //fixing the current heap
    curHeap = heapList[eli];
    ASSERT (curHeap != nullptr);
#ifndef NDEBUG
    for (Int k = 0; k < curHeap->size(); k++)
    {
        Int elid = (*curHeap)[k];
        if (elementList[eli] == nullptr)
            curHeap->erase(curHeap->begin()+k);
    }

#endif
    curHeap->push_back(eli);


#ifndef NDEBUG
    //Printing the contribution block after prior blocks assembly
    p = 1;
    PRLEVEL (p, ("\n%%After prior blocks assembly:"));
    if (p <= 0)
        paru_print_element (paruMatInfo, eli);
#endif



    Work->rowMark +=  rowCount;
    rowMark = Work -> rowMark;

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
    return 0;
}
