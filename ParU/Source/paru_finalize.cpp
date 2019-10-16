/** =========================================================================  /
 * =======================  paru_finalize ===================================  /
 * ========================================================================== */
/*! @brief numerical assemble of prior fronts
 *    prior fronts
 * @author Aznaveh
 */

#include "Parallel_LU.hpp"

void paru_finalize (paru_matrix *paruMatInfo, Int f, Int start_fac, 
        cholmod_common *cc){

    DEBUGLEVEL(1);
#ifndef NDEBUG
    Int p = 1;
    // counters to check the status of tuples scanning
    static Int f1 = 0, f2 = 0, f3 = 0, f4 = 0;
    // counters to find if it is children or their progeny that are fully or
    // partially assembled
    static Int child_FA = 0, noChild_FA = 0, 
               child_rowPA = 0, noChild_rowPA = 0,  
               elRow_child_rowPA = 0, elRow_noChild_rowPA = 0,  
               child_colPA= 0,  noChild_colPA= 0,
               elRow_child_colPA= 0,  elRow_noChild_colPA= 0;
#endif

    work_struct *Work =  paruMatInfo->Work;

    Int *isRowInFront = Work->rowSize; 
    Int rowMark = Work->rowMark;

    // Couning how many rows/cols of an element is seen
    Int *elRow = Work -> elRow; 
    Int *elCol = Work -> elCol;

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int *rM = LUsym->row2atree;
    Int *aParent = LUsym->aParent; //augmented tree size m+nf
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

#ifndef NDEBUG
#if 0
    PRLEVEL (p, ("%%*** List of children of %ld\n%%", el_ind));
    paru_Element *el = curFr;
    Int next = curFr-> next;
    while (next >= 0) {
        el->next = -1; 
        el = elementList[next];
        PRLEVEL (p, (" %ld -", next));
        ASSERT (elRow[next] == 0 && elCol [next] == 0);

        Int mEl = el->nrows;
        Int nEl = el->ncols;

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
        elementList[next] = NULL;
        PRLEVEL (0, ("%% $$ NULLIFIED\n"));
        elRow[next] = -1;
        elCol[next] = -1;

        next = elementList[next]->next;

        continue;

    }
#endif
#endif



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
            PRLEVEL (1, ("%% element= %ld  f =%ld \n",e, curColIndex));

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
            PRLEVEL (1, ("%% elRow[%ld]=%ld currVal= %ld, curcVal=%ld ", 
                        e, elRow[e], el->rValid , el->cValid));
            PRLEVEL (1, ("%% time_f =%ld ", time_f));
            PRLEVEL (1, ("%% start_fac=%ld \n", start_fac));

            if (elRow [e] == 0) {
               
                ASSERT ( el->rValid > start_fac);
                //ASSERT ( el->rValid < time_f);

                //all the columns are in CB
                if (elCol[e] == 0 ){
               
                ASSERT ( el->cValid > start_fac);
                ASSERT ( el->cValid < time_f);

                    // Whole prior front assembly
                    // do complete assembly of e into current front, now
                    PRLEVEL (1, ("%% element %ld is going to be eliminated\n",
                                e));
                    paru_update_rel_ind (curFr, el, 'r',cc) ;
                    paru_update_rel_ind (curFr, el, 'c', cc) ;

#ifndef NDEBUG            
                    p = 1;
                    PRLEVEL (p, ("%% Full assembly from: \n"));
                    if (p <= 0 ){
                        paru_print_element (paruMatInfo, e);
                    }
#endif

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
#ifndef NDEBUG
                    if (aParent [e] = el_ind )
                        child_FA++ ;
                    else
                        noChild_FA++ ;

#endif
                    elementList[e] = NULL;
                    PRLEVEL (1, ("%% NULLIFIED\n"));
                    continue;
                }

                ASSERT ( e < el_ind && e >= first[el_ind]);

                if(el->rValid !=  time_f){ /*  Update rowRelIndex	 */


#ifndef NDEBUG
                    p = 0;
                    PRLEVEL (p, ("%% update row relative element %ld\n", e ));
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
                    PRLEVEL (p, ("\n%%Before column assembly of %ld:\n",e));
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

                if (el->ncolsleft == 0){ //free el
                    Int tot_size = sizeof(paru_Element) +
                        sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
                    paru_free (1, tot_size, el, cc);
                    elementList[e] = NULL;
                }


#ifndef NDEBUG
                //Printing the contribution block after 
                //  prior blocks assembly
                p = 0;
                PRLEVEL (p, ("\n%%After column assembly of %ld:\n",e));
                if (p <= 0){
                    paru_print_element (paruMatInfo, el_ind);
                    paru_print_element (paruMatInfo, e);
                }

                if (aParent [e] = el_ind ){
                    if (mEl == 1 && rM[el_rowIndex[0]] == e ) //elRow
                        elRow_child_rowPA++;
                    else
                        child_rowPA++ ;
                }
                else{
                    if (mEl == 1 && rM[el_rowIndex[0]] == e ) //elRow
                        elRow_noChild_rowPA++ ;
                    else
                        noChild_rowPA++ ;
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
                //Not any of descendents; clear elRow if changed not to mangle
                //with other's computation
                elRow [e] = -1;
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

//#if 0
            if (elCol [e] == 0){ // && el->cValid > start_fac){
                //all the rows are in CB

                if (elRow[e] == 0){
                    PRLEVEL (1, ("%% %ld is already assembled \n",e ));
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
                PRLEVEL (p, ("%% after row assembly: \n" ));
                if (p <= 0 ){
                    paru_print_element (paruMatInfo, e);
                    paru_print_element (paruMatInfo, el_ind);
                }

                if (aParent [e] = el_ind ){
                    if (mEl == 1 && rM[el_rowIndex[0]] == e ) //elRow
                        elRow_child_colPA++ ;
                    else
                        child_colPA++ ;
                }
                else{
                    if (mEl == 1 && rM[el_rowIndex[0]] == e ) //elRow
                        elRow_child_colPA++ ;
                    else
                        noChild_colPA++ ;
                }
#endif
                rowRelIndex [curRowIndex] = -1;
                el_rowIndex [curRowIndex] = -1;
                el->nrowsleft --;
                if (el->nrowsleft == 0){ //free el
                    Int tot_size = sizeof(paru_Element) +
                        sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
                    paru_free (1, tot_size, el, cc);
                    elementList[e] = NULL;
                }


            } 
            //#endif
            PRLEVEL (1, ("%% RESET elRow[%ld]=%ld \n",e,  elRow[e]));
            elRow [e] = -1;
            PRLEVEL (1, ("%% After SET elRow[%ld]=%ld \n",e,  elRow[e]));
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
            p = 1;
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
            //el->cValid = -1;

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
        p = 1;
        if (p <= 0)
            paru_print_tupleList (ColList, c);
#endif

    }
    /**************************************************************************/


#ifndef NDEBUG            

    PRLEVEL (p, ("%% Finalized counters f1=%ld f2=%ld f3=%ld f4=%ld"
                " sum=%ld\n", f1, f2, f3, f4, f1+f2+f3+f4));

    p = 0;
    //if(f == LUsym-> nf -2)
        PRLEVEL (p, ("%% child_FA=%ld noChild_FA=%ld"
                    "  child_colPA=%ld noChild_colPA=%ld\n" 
                    "  elRow_child_colPA=%ld elRow_noChild_colPA=%ld\n" 
                    "  child_rowPA=%ld noChild_rowPA=%ld\n" 
                    "  elRow_child_rowPA=%ld elRow_noChild_rowPA=%ld\n" ,
                    child_FA, noChild_FA, child_colPA, noChild_colPA, 
                    elRow_child_colPA, elRow_noChild_colPA, 
                    child_rowPA, noChild_rowPA,
                    elRow_child_rowPA, elRow_noChild_rowPA 
                    ));

#endif

    // free the sorting space if allocated
    paru_free ( 2*curFr->nrows, sizeof(Int), curFr->rWork, cc); 
    curFr->rWork = NULL;
    paru_free ( 2*curFr->ncols, sizeof(Int), curFr->cWork, cc); 
    curFr->cWork = NULL;

}
