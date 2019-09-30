/** =========================================================================  /
 * =======================  paru_finalize ===================================  /
 * ========================================================================== */
/*! @brief numerical assemble of prior fronts
 *    prior fronts
 * @author Aznaveh
 */

#include "Parallel_LU.hpp"

void paru_finalize (paru_matrix *paruMatInfo, Int f, cholmod_common *cc){

    DEBUGLEVEL(1);
#ifndef NDEBUG
    Int p = 0;
    // counters to check the status of tuples scanning
    static Int f1 = 0, f2 = 0, f3 = 0, f4 = 0;
#endif

    work_struct *Work =  paruMatInfo->Work;

    Int *isRowInFront = Work->rowSize; 
    Int rowMark = Work->rowMark;

    // Couning how many rows/cols of an element is seen
    Int *elRow = Work -> elRow; 
    Int *elCol = Work -> elCol;

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int el_ind = snM [f]; 

    Int *Super = LUsym->Super;
    Int col1 = Super [f];     
    Int col2 = Super [f+1];
    Int fp = col2-col1;


    paru_Element **elementList = paruMatInfo->elementList;
    paru_Element *curFr = elementList[el_ind]; 
    Int rowCount= curFr->nrows + fp;
    Int colCount = curFr->ncols;

    Int *fcolList = paruMatInfo->fcolList[f] ;

    Int *isColInCBcolSet = Work -> colSize;
    Int colMark = Work -> colMark;

    Int *first = LUsym->first;
    Int *row_degree_bound = paruMatInfo->row_degree_bound;


    Int time_f = ++paruMatInfo->time_stamp[f]; //making all the markings invalid 


    double *cur_Numeric = numeric_pointer (curFr);
    Int new_row_degree_bound_for_r ;
    Int *frowList = paruMatInfo->frowList[f] ;
    /**************************************************************************/


    /****************************1st pass: assemble columns********************/
    tupleList *ColList = paruMatInfo->ColList;
    for (Int k = 0; k < colCount; k++){
        Int c = fcolList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG            
        p = 1;
        
        PRLEVEL (p, ("\n %%-------->  3rd: c =%ld  numTuple = %ld\n",
                    c, numTuple));
        if (p <= 0 ){
            paru_print_tupleList (ColList, c);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif

        for (Int i = 0; i < numTuple; i++){
            paru_Tuple curTpl = listColTuples [i];
            Int e = curTpl.e;
#ifndef NDEBUG
            f1++;
#endif
            if ( e >= el_ind || e < first[el_ind]){ 
                //Not any of descendents
                continue;
            }

            Int curColIndex = curTpl.f;
            PRLEVEL (0, ("%% element= %ld  f =%ld \n",e, curColIndex));

            ASSERT (e >= 0);
            ASSERT (e != el_ind);
            ASSERT (curColIndex >= 0);

            paru_Element *el = elementList[e];
            if (el == NULL) continue;
            Int mEl = el->nrows;
            Int nEl = el->ncols;

            Int *el_colIndex = colIndex_pointer (el);
            Int *el_rowIndex = rowIndex_pointer (el);
            Int *rowRelIndex = relRowInd (el);
            Int *colRelIndex    = relColInd (el);

            if (el_colIndex [curColIndex] < 0 ){ //already assembled somewhere
                continue;  
            }

            ASSERT (el_colIndex[curColIndex] == c);
            ASSERT (curColIndex < nEl);

            double *el_Num = numeric_pointer (el);
            PRLEVEL (1, ("%% elRow[%ld]=%ld currVal= %ld ", 
                        e, elRow[e], el->rValid ));
            PRLEVEL (1, ("%% time_f =%ld \n", time_f));

            if (elRow [e] == 0 ){         //all the columns are in CB


                if (elCol[e] == 0) {    // Whole prior front assembly
                    // do complete assembly of e into current front, now
                    PRLEVEL (0, ("%% element %ld is going to be eliminated\n",
                                e));
                    paru_update_rel_ind (curFr, el, 'r',cc) ;
                    paru_update_rel_ind (curFr, el, 'c', cc) ;

                    Int *rowRelIndex = relRowInd (el);
                    Int *colRelIndex = relColInd (el);

                    Int *el_rowIndex = rowIndex_pointer (el);
                    Int *el_colIndex = colIndex_pointer (el);

                    double *el_Num = numeric_pointer (el);
                    assemble_all (el_Num, cur_Numeric,
                            el->nrows, el->ncols, curFr->nrows,
                            rowRelIndex, colRelIndex);
                    // delete e
                    Int tot_size = sizeof(paru_Element) +
                        sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
                    paru_free (1, tot_size, el, cc);
                    elementList[e] = NULL;
                    PRLEVEL (0, ("%% NULLIFIED\n"));
                    continue;
                }

                ASSERT ( e < el_ind && e >= first[el_ind]);

                if(el->rValid !=  time_f){ /*  Update rowRelIndex	 */


#ifndef NDEBUG
                    p = 1;
                    PRLEVEL (p, ("%% update row relative element%ld\n", e ));
                    //Printing the contribution block prior index update 
                    if (p <= 0) {
                        PRLEVEL (p, ("\n%%Before index update %ld:",e));
                        paru_print_element (paruMatInfo, e);
                    }
#endif
                    paru_update_rel_ind (curFr, el, 'r',cc) ;
#ifndef NDEBUG            
                    for(Int i=0; i < el->nrows; i++){
                        PRLEVEL (1, ("%% rowRelIndex[%ld] =%ld\t", i,
                                    rowRelIndex [i]));
                        ASSERT(rowRelIndex [i] < curFr->nrows);
                        PRLEVEL (1,("\n"));
                    }
#endif

                    el->rValid =  time_f;
                }
#ifndef NDEBUG
                //Printing the contribution block before 
                //   prior blocks assembly
                p = 1;
                if (p <= 0){
                    PRLEVEL (p, ("\n%%Before column assembly of %ld:",e));
                    paru_print_element (paruMatInfo, el_ind);
                    paru_print_element (paruMatInfo, e);
                }

                PRLEVEL (p, ("%%colCount=%ld k=%ld", colCount, k));
                PRLEVEL (p, ("%%curFr->nrows=%ld ", curFr->nrows));
                PRLEVEL (p, ("%% cur_Numeric=%2.4lf\n",
                            *(cur_Numeric+k*curFr->nrows)));
#endif

                assemble_col (el_Num+curColIndex*mEl,
                        cur_Numeric+k*curFr->nrows, mEl, rowRelIndex);
                colRelIndex [curColIndex] = -1;
                el_colIndex [curColIndex] = -1;
                el->ncolsleft --;
#ifndef NDEBUG
                //Printing the contribution block after 
                //  prior blocks assembly
                p = 1;
                PRLEVEL (p, ("\n%%After column assembly of %ld:",e));
                if (p <= 0){
                    paru_print_element (paruMatInfo, el_ind);
                    paru_print_element (paruMatInfo, e);
                }
#endif

            } 
        }

#ifndef NDEBUG
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif

    }


    /**************************************************************************/

    /*******************       2nd pass: assemble rows            *************/
    tupleList *RowList = paruMatInfo->RowList;
    for (Int k = fp; k < rowCount; k++){
        Int r = frowList [k];
        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG            
        p = 1;
        PRLEVEL (1, ("\n %%------->  4th: r =%ld  numTuple = %ld\n",
                    r, numTuple));
        if (p <= 0 ){
            paru_print_tupleList (RowList, r);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif

        for (Int i = 0; i < numTuple; i++){
            paru_Tuple curTpl = listRowTuples [i];
            Int e = curTpl.e;
#ifndef NDEBUG
            f2++;
#endif
 
            if ( e >= el_ind || e < first[el_ind]){
                //Not any of descendents
                continue;
            }
            //TODO: these assertion doesn't work: why?
            // ASSERT ( e >= el_ind );
            // ASSERT ( e < first[el_ind] );

            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",
                        e, curColIndex));


            Int curRowIndex = curTpl.f;
            ASSERT (e >= 0);
            ASSERT (curRowIndex >= 0);

            paru_Element *el = elementList[e];
            if (el == NULL) continue;
            Int mEl = el->nrows;
            Int nEl = el->ncols;

            Int *el_colIndex = colIndex_pointer (el);
            Int *el_rowIndex = rowIndex_pointer (el);
            Int *rowRelIndex = relRowInd (el);
            Int *colRelIndex    = relColInd (el);


            if (el_rowIndex [curRowIndex] < 0 ){ 
                // it will be deleted here
                continue;  
            }
            if (el_rowIndex[curRowIndex] < 0)  
                continue; //not valid
            PRLEVEL (1, ("%% el_rowIndex [%ld] =%ld\n", 
                        curRowIndex, el_rowIndex [curRowIndex]));
            ASSERT (el_rowIndex[curRowIndex] == r);
            ASSERT (curRowIndex < mEl);

            double *el_Num = numeric_pointer (el);
            PRLEVEL (1, ("%% elCol[%ld]=%ld ",e, elCol[e]));

#if 1
            if (elCol [e] == 0 ){ //all the rows are in CB

                if (elRow[e] == 0){
                    continue; //already assembled
                }
                ASSERT ( e < el_ind && e >= first[el_ind]);

#ifndef NDEBUG            
                p = 1;
                PRLEVEL (1, ("%% Before row assembly: \n" ));
                if (p <= 0 ){
                    paru_print_element (paruMatInfo, e);
                    paru_print_element (paruMatInfo, el_ind);
                }
#endif
                if(el->cValid !=  time_f){
                    /* Update colRelIndex	 */
                    PRLEVEL (1, ("%% update column relative index %ld\n"
                                ,e ));
                    paru_update_rel_ind (curFr, el, 'c', cc) ;
#ifndef NDEBUG            
                    for(Int i=0 ; i <el->ncols ; i++){
                        PRLEVEL (1, ("%% colRelIndex[%ld] =%ld\t", i,
                                    colRelIndex [i]));
                        ASSERT(colRelIndex [i] < curFr->ncols);
                        PRLEVEL (1,("\n"));
                    }
#endif
                    el->cValid =  time_f;
                }
                assemble_row (el_Num, cur_Numeric, mEl, nEl, 
                        curFr->nrows, curRowIndex , k-fp, colRelIndex );
#ifndef NDEBUG            
                p = 1;
                PRLEVEL (1, ("%% after row assembly: \n" ));
                if (p <= 0 ){
                    paru_print_element (paruMatInfo, e);
                    paru_print_element (paruMatInfo, el_ind);
                }
#endif
                rowRelIndex [curRowIndex] = -1;
                el_rowIndex [curRowIndex] = -1;
                el->nrowsleft --;
            } 
#endif
            elRow [e] = -1;
        }

    }


    /********************* 3rd path: clearing column tuples and uncheck *******/
    for (Int k = 0; k < colCount; k++){
        Int c = fcolList [k];   //non pivotal column list
        tupleList *curColTupleList = &ColList[c];
        Int numTuple = curColTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listColTuples = curColTupleList->list;
#ifndef NDEBUG            
        p = 0;
        PRLEVEL (p, ("\n %%-------->  5th: c =%ld  numTuple = %ld\n",
                    c, numTuple));
        if (p <= 0 ){
            paru_print_tupleList (ColList, c);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif
        Int pdst = 0, psrc;

        for (psrc = 0; psrc < numTuple; psrc ++){
            paru_Tuple curTpl = listColTuples [psrc];
            Int e = curTpl.e;
            Int curColIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curColIndex));

#ifndef NDEBUG
            f3++;
#endif
            paru_Element *el = elementList[e];
            if (el == NULL){
                PRLEVEL (1, ("%% El==NULL\n"));
                continue;
            }
            Int mEl = el->nrows;
            Int nEl = el->ncols;

            Int *el_colIndex = colIndex_pointer (el);
            // elCol[e]= -1;

            if (el_colIndex [curColIndex] < 0 ){ //it will be deleted here
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
        p = 0;
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif

    }
    /**************************************************************************/

#if 0
    /************** 4th path: clearing row tuples and unhceck elRow************/
    for (Int k = fp; k < rowCount; k++){
        Int r = frowList [k];
        tupleList *curRowTupleList = &RowList[r];
        Int numTuple = curRowTupleList->numTuple;
        ASSERT (numTuple >= 0);
        paru_Tuple *listRowTuples = curRowTupleList->list;
#ifndef NDEBUG            
        p = 1;
        PRLEVEL (p, ("\n %%-------->  6th: r =%ld  numTuple = %ld\n",
                    r, numTuple));
        if (p <= 0 ){
            paru_print_tupleList (RowList, r);
            paru_print_element (paruMatInfo, el_ind);
        }
#endif
        Int pdst = 0, psrc;

        for (psrc = 0; psrc < numTuple; psrc ++){
            paru_Tuple curTpl = listRowTuples [psrc];
            Int e = curTpl.e;
            Int curRowIndex = curTpl.f;
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curRowIndex));
#ifndef NDEBUG
            f4++;
#endif


            paru_Element *el = elementList[e];
            if (el == NULL) continue;
            Int mEl = el->nrows;
            Int nEl = el->ncols;

            Int *el_rowIndex = rowIndex_pointer (el);
            elRow[e]= -1;

            if (el_rowIndex [curRowIndex] < 0 ){ //it will be deleted here
                PRLEVEL (1, ("%% psrc=%ld\n", psrc));
                continue;  
            }
            else 
                listRowTuples [pdst++] = curTpl; //keeping the tuple
            PRLEVEL (1, ("%% time_f =%ld \n", time_f));

        }
        curRowTupleList->numTuple = pdst;
        PRLEVEL (1, ("%% new number of tuples=%ld\n", pdst));


#ifndef NDEBUG            
        p = 0;
        if (p <= 0 )
            paru_print_tupleList (RowList, r);

        PRLEVEL (p, ("%% Finalized counters f1=%ld f2=%ld f3=%ld f4=%ld"
                    " sum=%ld\n", f1, f2, f3, f4, f1+f2+f3+f4));
#endif
    }
#endif

    // free the sorting space if allocated
    paru_free ( 2*curFr->nrows, sizeof(Int), curFr->rWork, cc); 
    curFr->rWork = NULL;
    paru_free ( 2*curFr->ncols, sizeof(Int), curFr->cWork, cc); 
    curFr->cWork = NULL;

}
