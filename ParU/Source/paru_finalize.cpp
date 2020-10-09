/** =========================================================================  /
 * =======================  paru_finalize ===================================  /
 * ========================================================================== */
/*! @brief numerical assemble of prior fronts
 *    prior fronts
 * @author Aznaveh
 */

#include "Parallel_LU.hpp"

void paru_finalize (paru_matrix *paruMatInfo, Int f, Int start_fac,  
        cholmod_common *cc)
{

    DEBUGLEVEL(0);
#ifndef NDEBUG
    Int p = 1;
    // counters to check the status of tuples scanning
    static Int f1 = 0, f2 = 0, f3 = 0 ;
    // counters to find if it is children or their progeny that are fully or
    // partially assembled
    static Int child_FA = 0, noChild_FA = 0, 
               child_rowPA = 0, noChild_rowPA = 0,  
               elRow_child_rowPA = 0, elRow_noChild_rowPA = 0,  
               child_colPA= 0,  noChild_colPA= 0,
               elRow_child_colPA= 0,  elRow_noChild_colPA= 0;
#endif

    work_struct *Work =  paruMatInfo->Work;

    // Couning how many rows/cols of an element is seen
    Int *elRow = Work -> elRow; 
    Int *elCol = Work -> elCol;

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int el_ind = snM [f]; 
#ifndef NDEBUG
    Int *rM = LUsym->row2atree;
    Int *aParent = LUsym->aParent; //augmented tree size m+nf
#endif

    Int *Super = LUsym->Super;
    Int col1 = Super [f];     
    Int col2 = Super [f+1];
    Int fp = col2-col1;


    paru_Element **elementList = paruMatInfo->elementList;
    paru_Element *curFr = elementList[el_ind]; 
    Int curFrNrows = curFr->nrows;
    Int rowCount= curFrNrows + fp;
    Int colCount = curFr->ncols;

    Int *fcolList = paruMatInfo->fcolList[f] ;


    Int *first = LUsym->first;

    Int time_f = ++paruMatInfo->time_stamp[f]; //making all the markings invalid


    // double *cur_Numeric = numeric_pointer (curFr);
    double *cur_Numeric = 
        (double*)((Int*)(curFr+1) + 2*colCount + 2*curFrNrows);
    /**************************************************************************/
    /****************************1st pass: assemble columns********************/
    tupleList *ColList = paruMatInfo->ColList;
    for (Int k = 0; k < colCount; k++)
    {
        Int c = fcolList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG            
        p = 1;

        PRLEVEL (p, ("\n %%-------->  3rd: c =%ld  numTuple = %ld\n",
                    c, numTuple));
        if (p <= 0 )
        {
            paru_print_tupleList (ColList, c);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif

        for (Int i = 0; i < numTuple; i++)
        {
            paru_Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
#ifndef NDEBUG
            f1++;
#endif
            if ( e >= el_ind || e < first[el_ind])
            { 
                //Not any of descendents
                continue;
            }

            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curColIndex));

            ASSERT (e >= 0);
            ASSERT (e != el_ind);
            ASSERT (curColIndex >= 0);

            paru_Element *el = elementList[e];
            if (el == NULL) continue;


            //Int *el_colIndex = colIndex_pointer (el);
            Int *el_colIndex = (Int*)(el+1);
            
            if (el_colIndex [curColIndex] < 0 )
             //already assembled somewhere
                continue;  
            

            ASSERT (el_colIndex[curColIndex] == c);
            Int nEl = el->ncols;
            Int mEl = el->nrows;
            ASSERT (curColIndex < nEl);

            //double *el_Num = numeric_pointer (el);
            double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);
            
            PRLEVEL (1, ("%% elRow[%ld]=%ld currVal= %ld, curcVal=%ld ", 
                        e, elRow[e], el->rValid , el->cValid));
            PRLEVEL (1, ("%% time_f =%ld ", time_f));
            PRLEVEL (1, ("%% start_fac=%ld \n", start_fac));

            if (elRow [e] == 0) 
            {

                ASSERT ( el->rValid >= start_fac);

                //all the columns are in CB
                if (elCol[e] == 0 )
                {

                    ASSERT ( el->cValid > start_fac);
                    ASSERT ( el->cValid < time_f);

                    // Whole prior front assembly
                    // do complete assembly of e into current front, now
                    PRLEVEL (1, ("%% element %ld is going to be eliminated\n",
                                e));
                    paru_update_rel_ind_row (curFr, el, cc) ;
                    paru_update_rel_ind_col (paruMatInfo, f, curFr, el, cc) ;
                    //Int *rowRelIndex = relRowInd (el);
                    Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;
                    //Int *colRelIndex = relColInd (el);
                    Int *colRelIndex = (Int*)(el+1) + mEl + nEl;

#ifndef NDEBUG            
                    p = 1;
                    PRLEVEL (p, ("%% Full assembly from:"));
                    if (p <= 0 )
                        paru_print_element (paruMatInfo, e);
                    PRLEVEL (p, ("%% rowRelIndex:\n%%"));
                    for (int i=0 ; i< el->nrows; i++)
                        PRLEVEL (p, (" %ld",rowRelIndex[i]));
                    PRLEVEL (p, ("\n%% colRelIndex:\n%%"));
                    for (int i=0 ; i< el->ncols; i++)
                        PRLEVEL (p, (" %ld",colRelIndex[i]));
                    PRLEVEL (p, ("\n"));

#endif


                    assemble_all (el_Num, cur_Numeric, mEl, nEl, curFrNrows,
                            el->nrowsleft, el->ncolsleft, rowRelIndex, 
                            colRelIndex);
                    // delete e
                    Int tot_size = sizeof(paru_Element) +
                        sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
                    paru_free (1, tot_size, el, cc);


#ifndef NDEBUG
                    if (aParent [e] = el_ind )
                        child_FA++ ;
                    else
                        noChild_FA++ ;
                    p = 1;
                    if (p <= 0) 
                    {
                        
                        PRLEVEL (p, ("\n%%After all assembly of %ld:",e));

                        PRLEVEL (p, ("\n%%Source%ld:",e));
                        paru_print_element (paruMatInfo, e);
 
                        PRLEVEL (p, ("\n%%Destin%ld:",el_ind));
                        paru_print_element (paruMatInfo, el_ind);
                    }

#endif
                    elementList[e] = NULL;
                    continue;
                }


                if(el->rValid !=  time_f)
                { 
                    /*  Update rowRelIndex just once	 */

#ifndef NDEBUG
                    p = 1;
                    PRLEVEL (p, ("%% update row relative element %ld\n", e ));
                    //Printing the contribution block prior index update 
                    if (p <= 0) 
                    {
                        PRLEVEL (p, ("\n%%Before index update %ld:",e));
                        paru_print_element (paruMatInfo, e);
                    }
#endif
                    paru_update_rel_ind_row (curFr, el, cc) ;
#ifndef NDEBUG            
                    Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;
                    for(Int i=0; i < el->nrows; i++)
                    {
                        PRLEVEL (1, ("%% rowRelIndex[%ld] =%ld\t", i,
                                    rowRelIndex [i]));
                        ASSERT(rowRelIndex [i] < curFrNrows);
                        PRLEVEL (1,("\n"));
                    }
#endif

                    el->rValid =  time_f;
                }
#ifndef NDEBUG
                //Printing the contribution block before 
                //   prior blocks assembly
                p = 1;
                if (p <= 0)
                {
                    PRLEVEL (p, ("\n%%Before column assembly of %ld:\n",e));
                    paru_print_element (paruMatInfo, el_ind);
                    paru_print_element (paruMatInfo, e);
                }

                PRLEVEL (p, ("%%colCount=%ld k=%ld", colCount, k));
                PRLEVEL (p, ("%%curFr->nrows=%ld ", curFrNrows));
                PRLEVEL (p, ("%% cur_Numeric=%2.4lf\n",
                            *(cur_Numeric+k*curFrNrows)));
#endif
                //Int *rowRelIndex = relRowInd (el);
                Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;

                //Int *colRelIndex    = relColInd (el);
                Int *colRelIndex = (Int*)(el+1) + mEl + nEl;

                assemble_col (el_Num+curColIndex*mEl,
                        cur_Numeric+k*curFrNrows, mEl, rowRelIndex);
                colRelIndex [curColIndex] = -1;
                //el_colIndex [curColIndex] = -1;
                el_colIndex[curColIndex] = flip (el_colIndex[curColIndex] );
                el->ncolsleft --;

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
                    if (curColIndex <= el->lac)
                    {
                        while(el_colIndex[el->lac] < 0)
                        {
                            PRLEVEL (1, ("\n%%el->lac= %ld ",el->lac));
                            PRLEVEL (1, ("el_colIndex[el->lac]=%ld :\n",
                                        el_colIndex[el->lac]));
                            el->lac++;
                        }
                    }
                    PRLEVEL (1, ("%%Final: curColIndex=%ld ",curColIndex));
                    PRLEVEL (1, ("%%ncols=%ld ",el->ncols));
                    PRLEVEL (1, (" el->lac %ld ",el->lac));
                    PRLEVEL (1, (" el_colIndex[el->lac]=%ld\n",
                                el_colIndex[el->lac]));
                    ASSERT (el->lac < el->ncols);
                }




#ifndef NDEBUG
                Int *el_rowIndex = rowIndex_pointer (el);
                //Printing the contribution block after 
                //  prior blocks assembly
                p = 1;
                PRLEVEL (p, ("\n%%After column assembly of %ld:\n",e));
                if (p <= 0)
                {
                    paru_print_element (paruMatInfo, el_ind);
                    paru_print_element (paruMatInfo, e);
                }

                if (aParent [e] = el_ind )
                {
                    if (mEl == 1 && rM[el_rowIndex[0]] == e ) //elRow
                        elRow_child_colPA++;
                    else
                        child_colPA++ ;
                }
                else
                {
                    if (mEl == 1 && rM[el_rowIndex[0]] == e ) //elRow
                        elRow_noChild_colPA++ ;
                    else
                        noChild_colPA++ ;
                }
#endif

            } 
        }

#ifndef NDEBUG
        p = 1;
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif

    }


    /**************************************************************************/

    /*******************       2nd pass: assemble rows            *************/

    tupleList *RowList = paruMatInfo->RowList;
    for (Int k = fp; k < rowCount; k++)
    {
        Int *frowList = paruMatInfo->frowList[f] ;
        Int r = frowList [k];
        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG            
        p = 1;
        PRLEVEL (1, ("\n %%------->  4th: r =%ld  numTuple = %ld\n",
                    r, numTuple));
        if (p <= 0 )
        {
            paru_print_tupleList (RowList, r);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif

        // for (Int i = 0; i < numTuple; i++)
        // {
        //     paru_Tuple curTpl = listRowTuples [i];
        Int pdst = 0, psrc;
        for (psrc = 0; psrc < numTuple; psrc ++)
        {
            paru_Tuple curTpl = listRowTuples [psrc];
            Int e = curTpl.e;
#ifndef NDEBUG
            f2++;
#endif

            if ( e >= el_ind || e < first[el_ind])
            {
                //Not any of descendents; clear elRow if changed not to mangle
                //with other's computation
                elRow [e] = -1;
                listRowTuples [pdst++] = curTpl; //keeping the tuple
                continue;
            }

            paru_Element *el = elementList[e];
            if (el == NULL) continue;

            Int curRowIndex = curTpl.f;

            Int nEl = el->ncols;
            //Int *el_rowIndex = rowIndex_pointer (el);
            Int *el_rowIndex = (Int*) (el+1) + nEl;

            if (el_rowIndex [curRowIndex] < 0 )
                // it will be deleted here
                continue;  

            PRLEVEL (1, ("%% elCol[%ld]=%ld ",e, elCol[e]));


            if (elCol [e] == 0)
            { 
                //all the rows are in CB

                ASSERT (el_rowIndex[curRowIndex] == r);
                ASSERT (elRow[e] != 0);

#ifndef NDEBUG            
                p = 1;
                PRLEVEL (1, ("%% Before row assembly: \n" ));
                PRLEVEL (1, ("%% elRow[%ld]=%ld ",e, elRow[e]));
                if (p <= 0 )
                {
                    paru_print_element (paruMatInfo, e);
                    paru_print_element (paruMatInfo, el_ind);
                }
#endif
                Int mEl = el->nrows;
                //Int *colRelIndex = relColInd (el);
                Int *colRelIndex = (Int*)(el+1) + nEl + mEl;

                if(el->cValid !=  time_f)
                {
                    /* Update colRelIndex	 */
                    PRLEVEL (1, ("%% update column relative index %ld\n"
                                ,e ));
                    paru_update_rel_ind_col (paruMatInfo, f, curFr, el, cc) ;
#ifndef NDEBUG            
                    for(Int i=0 ; i <el->ncols ; i++)
                    {
                        PRLEVEL (1, ("%% colRelIndex[%ld] =%ld\t", i,
                                    colRelIndex [i]));
                        ASSERT(colRelIndex [i] < curFr->ncols);
                        PRLEVEL (1,("\n"));
                    }
#endif
                    el->cValid =  time_f;
                }

                //double *el_Num = numeric_pointer (el);
                double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);

                assemble_row (el_Num, cur_Numeric, mEl, nEl, 
                        curFrNrows, curRowIndex , k-fp, colRelIndex );
#ifndef NDEBUG            
                p = 1;
                PRLEVEL (p, ("%% after row assembly: \n" ));
                if (p <= 0 )
                {
                    paru_print_element (paruMatInfo, e);
                    paru_print_element (paruMatInfo, el_ind);
                }

                if (aParent [e] = el_ind )
                {
                    if (mEl == 1 && rM[el_rowIndex[0]] == e ) //elRow
                        elRow_child_rowPA++ ;
                    else
                        child_rowPA++ ;
                }
                else
                {
                    if (mEl == 1 && rM[el_rowIndex[0]] == e ) //elRow
                        elRow_child_rowPA++ ;
                    else
                        noChild_rowPA++ ;
                }
#endif
                //Int *rowRelIndex = relRowInd (el);
                Int *rowRelIndex = (Int*)(el+1) + 2*nEl+ mEl;

                rowRelIndex [curRowIndex] = -1;
                el_rowIndex [curRowIndex] = -1;
                el->nrowsleft --;
                if (el->nrowsleft == 0)
                { //free el
                    Int tot_size = sizeof(paru_Element) +
                        sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
                    paru_free (1, tot_size, el, cc);
                    elementList[e] = NULL;
                }


            } 
            else
                listRowTuples [pdst++] = curTpl; //keeping the tuple
            PRLEVEL (1, ("%% RESET elRow[%ld]=%ld \n",e,  elRow[e]));
            elRow [e] = -1;
            PRLEVEL (1, ("%% After SET elRow[%ld]=%ld \n",e,  elRow[e]));
        }
        curRowTupleList->numTuple = pdst;

    }


    /********************* 3rd path: clearing column tuples and uncheck *******/
    for (Int k = 0; k < colCount; k++)
    {
        Int c = fcolList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        paru_Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG            
        p = 1;
        PRLEVEL (p, ("\n %%-------->  5th: c =%ld  numTuple = %ld\n",
                    c, numTuple));
        if (p <= 0 )
        {
            paru_print_tupleList (ColList, c);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif
        Int pdst = 0, psrc;

        for (psrc = 0; psrc < numTuple; psrc ++)
        {
            paru_Tuple curTpl = listColTuples [psrc];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curColIndex));

#ifndef NDEBUG
            f3++;
#endif
            paru_Element *el = elementList[e];
            if (el == NULL)
            {
                PRLEVEL (1, ("%% El==NULL\n"));
                continue;
            }

            Int *el_colIndex = colIndex_pointer (el);

            if (el_colIndex [curColIndex] < 0 )
            { //it will be deleted here
                PRLEVEL (1, ("%% psrc=%ld\n", psrc));
                continue;  
            }
            else 
                listColTuples [pdst++] = curTpl; //keeping the tuple

            PRLEVEL (1, ("%% time_f =%ld \n", time_f));

        }
        curColTupleList->numTuple = pdst;
        PRLEVEL (1, ("%% new number of tuple=%ld\n", pdst));
#ifndef NDEBUG
        p = 1;
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif

    }
    /**************************************************************************/


#ifndef NDEBUG            

    PRLEVEL (p, ("%% Finalized counters f1=%ld f2=%ld f3=%ld "
                " sum=%ld\n", f1, f2, f3, f1+f2+f3));

    p = 1;
    //if(f == LUsym-> nf -2)
    PRLEVEL (p, ("%% child_FA=%ld noChild_FA=%ld"
                "  child_colPA=%ld noChild_colPA=%ld\n%%" 
                "  elRow_child_colPA=%ld elRow_noChild_colPA=%ld\n%%" 
                "  child_rowPA=%ld noChild_rowPA=%ld\n%%" 
                "  elRow_child_rowPA=%ld elRow_noChild_rowPA=%ld\n%%" ,
                child_FA, noChild_FA, child_colPA, noChild_colPA, 
                elRow_child_colPA, elRow_noChild_colPA, 
                child_rowPA, noChild_rowPA,
                elRow_child_rowPA, elRow_noChild_rowPA 
                ));

#endif

    // free the sorting space if allocated
    paru_free ( 2*curFrNrows, sizeof(Int), curFr->rWork, cc); 
    curFr->rWork = NULL;

}
