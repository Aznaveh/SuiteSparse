/** =========================================================================  /
 * =======================  paru_prior_assemble =============================  /
 * ========================================================================== */
/*! @brief numerical assemble of prior fronts
 *    
 * @author Aznaveh
 */

#include "Parallel_LU.hpp"

void paru_prior_assemble ( Int f, Int start_fac,  
        std::vector<Int> &pivotal_elements,
        paru_matrix *paruMatInfo,
        cholmod_common *cc)
{
    DEBUGLEVEL(1);

    work_struct *Work =  paruMatInfo->Work;
    Int *elRow = Work -> elRow; 
    Int *elCol = Work -> elCol;
    
    paru_Element **elementList = paruMatInfo->elementList;
    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
 
    Int el_ind = snM [f]; 
    paru_Element *curEl = elementList[el_ind]; 

    Int *Super = LUsym->Super;
    Int col1 = Super [f];     
    Int col2 = Super [f+1];
    Int fp = col2-col1;

    Int curElNrows = curEl->nrows;
    Int rowCount= curElNrows + fp;



    Int colCount = curEl->ncols;
    // double *cur_Numeric = numeric_pointer (curFr);
    double *cur_Numeric = 
        (double*)((Int*)(curEl+1) + 2*colCount + 2*curElNrows);

    Int *fcolList = paruMatInfo->fcolList[f] ;

    Int pMark = start_fac;

    PRLEVEL (1, ("%%Inside prior\n"));
        PRLEVEL (1, ("%% pivotal size is %ld ", pivotal_elements.size()));
    Int ii = 0;
    for(Int i = 0 ; i < pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        paru_Element *el = elementList[e];
        PRLEVEL (1, ("%% element= %ld  \n",e));
        if ( el == NULL) 
        {
            PRLEVEL (3, ("%% element= %ld is NULL ii=%ld \n",e, ii));
            continue;
        }
        PRLEVEL (1, ("%%elRow[%ld]=%ld \n",e, elRow[e]));
        ASSERT (elRow[e] == 0);
        if (el->rValid == pMark || elCol[e] == 0)
            // it can be eliminated fully
            // both a pivotal column and pivotal row
        {//TODO asselmble all and delete from the pivotal_elements

            paru_update_rel_ind_row (curEl, el, cc) ;
            paru_update_rel_ind_col (paruMatInfo, f, curEl, el, cc) ;

            Int nEl = el->ncols;
            Int mEl = el->nrows;

            //Int *rowRelIndex = relRowInd (el);
            Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;
            //Int *colRelIndex = relColInd (el);
            Int *colRelIndex = (Int*)(el+1) + mEl + nEl;

            //double *el_Num = numeric_pointer (el);
            double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);

            assemble_all (el_Num, cur_Numeric, mEl, nEl, curElNrows,
                    el->nrowsleft, el->ncolsleft, rowRelIndex, 
                    colRelIndex);
            // delete e
            Int tot_size = sizeof(paru_Element) +
                sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
            paru_free (1, tot_size, el, cc);
            PRLEVEL (1, ("%%Prior assembly Free %ld  %p size %ld\n",
                        e, el, tot_size));
            elementList[e] = NULL;
            continue; 
        }

//        if (el->rValid > start_fac && elCol[e] < el->ncolsleft)
//        { 
//            //TODO search for columns and asslemble them
//            paru_update_rel_ind_row (curEl, el, cc) ;
//        }
        pivotal_elements [ii++] = pivotal_elements [i];
    }

    if ( ii+1 < pivotal_elements.size())
    {
        PRLEVEL (1, ("%% Prior: size was %ld ", pivotal_elements.size()));
        PRLEVEL (1, ("%% and now is %ld\n ", ii+1));
        pivotal_elements.resize(ii+1);
    }



    Int eli = snM [f]; 
    std::vector<Int>** heapList = paruMatInfo->heapList;
    std::vector<Int>* curHeap = heapList[eli];

//    for(Int i = curHeap->size()-1; i > 0 ; i--)
//    {
//        Int e = (*curHeap)[i];
//        paru_Element *el = elementList[e];
//
//        if (el->rValid >= pMark && elRow[e] == 0)
//        { // all the rows are inside current fron
//            if (el->cValid >= pMark &&elCol[e] < el->ncolsleft)
//            {//asslemble all and delete from the heap
//                paru_update_rel_ind_row (curEl, el, cc) ;
//            }
//        }
//    }
//


}

