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
        std::unordered_map <Int, Int> colHash, 
        paru_matrix *paruMatInfo,
        cholmod_common *cc)
{
    DEBUGLEVEL(1);
#ifndef NDEBUG  
    Int p = 1;
#endif


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

    PRLEVEL (p, ("%%Inside prior\n"));
    PRLEVEL (p, ("%% pivotal size is %ld ", pivotal_elements.size()));
    Int ii = 0;
    //TODO the relative indecis are valid
    for(Int i = 0 ; i < pivotal_elements.size(); i++)
    {
        Int e = pivotal_elements[i];
        paru_Element *el = elementList[e];
        PRLEVEL (p, ("%% element= %ld  \n",e));
        if ( el == NULL) 
        {
            PRLEVEL (p, ("%% element= %ld is NULL ii=%ld \n",e, ii));
            continue;
        }
        PRLEVEL (p, ("%%elRow[%ld]=%ld \n",e, elRow[e]));
        ASSERT (elRow[e] == 0);
        if (el->rValid == pMark || elCol[e] == 0)
            // it can be eliminated fully
            // both a pivotal column and pivotal row
        {//TODO asselmble all and delete from the pivotal_elements
            paru_update_rel_ind_row (curEl, el, cc) ;
            PRLEVEL (p, ("%%pivotal element rel col update %ld \n", e));
            paru_update_rel_ind_col (paruMatInfo, f, curEl, el, cc) ;

            Int nEl = el->ncols;
            Int mEl = el->nrows;

            //Int *rowRelIndex = relRowInd (el);
            Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;
            //Int *colRelIndex = relColInd (el);
            Int *colRelIndex = (Int*)(el+1) + mEl + nEl;

            //double *el_Num = numeric_pointer (el);
            double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);

            PRLEVEL (p, ("%%assembling %ld in %ld\n", e, el_ind));
            PRLEVEL (p, ("%% size %ld x %ld\n", mEl, nEl));
            paru_eliminate (e, f, colHash, paruMatInfo);
//            assemble_all (el_Num, cur_Numeric, mEl, nEl, curElNrows,
//                    el->nrowsleft, el->ncolsleft, rowRelIndex, 
//                    colRelIndex);
            PRLEVEL (p, ("%%assembling %ld in %ld done\n", e, el_ind));
            // delete e
            Int tot_size = sizeof(paru_Element) +
                sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
            paru_free (1, tot_size, el, cc);
            PRLEVEL (p, ("%%Prior assembly Free %ld  %p size %ld\n",
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

    if ( ii < pivotal_elements.size())
    {
        PRLEVEL (p, ("%% Prior: size was %ld ", pivotal_elements.size()));
        PRLEVEL (p, (" and now is %ld\n ", ii));
        pivotal_elements.resize(ii);
    }


    Int eli = snM [f]; 
    std::vector<Int>** heapList = paruMatInfo->heapList;
    std::vector<Int>* curHeap = heapList[eli];

#ifndef NDEBUG  
    p = 1;
#endif


    Int *lacList = paruMatInfo->lacList;
#ifndef NDEBUG  
    PRLEVEL (p, ("%% current heap:\n %%"));
    for(Int k = 0 ; k < curHeap->size(); k++)
    {
        Int ee = (*curHeap)[k];
        paru_Element *ell = elementList[ee];
        PRLEVEL (p, ("%ld-%ld", k, ee));
        if (ell != NULL)
        {PRLEVEL (p, ("(%ld) ", lacList[ee] ));}
        else
        { PRLEVEL (p, ("(*%ld) ", lacList[ee] ));}
    }
    PRLEVEL (p, ("\n"));
#endif


    //if ( curHeap->size() == 0 ) return;
    if ( curHeap->empty() ) return;

    auto greater = [&lacList](Int a, Int b) {return lacList[a] > lacList[b]; };

    // chekc the root
    while ( elementList [ (*curHeap)[0] ] == NULL && curHeap->size() > 0 )
    {
        std::pop_heap (curHeap->begin(), curHeap->end(), greater );
        curHeap->pop_back();
    }

    //if ( curHeap->size() == 0 ) return;
    if ( curHeap->empty() ) return;

    for(Int i = curHeap->size()-1; i > 0 ; i--)
    {
        Int e = (*curHeap)[i];
        paru_Element *el = elementList[e];

        if ( el == NULL) 
        { //delete el from the heap
            remove_heap (i, lacList, (*curHeap));
#ifndef NDEBUG  
            PRLEVEL (p, ("%% Intermediate heap %ld deleted:\n %%", e))
            for(Int k = 0 ; k < curHeap->size(); k++)
            {
                Int ee = (*curHeap)[k];
                paru_Element *ell = elementList[ee];
                PRLEVEL (p, ("%ld-%ld", k, ee));
                if (ell != NULL)
                {PRLEVEL (p, ("(%ld) ", lacList[ee] ));}
                else
                { PRLEVEL (p, ("(*%ld) ", lacList[ee] ));}
            }
            PRLEVEL (p, ("\n"));
#endif
            continue;
        }

        if (elRow [e] == 0 && elCol [e] == 0 && el->rValid >= pMark)
        {
            PRLEVEL (-1, ("%% Inside the heap %ld deleted:\n %%", e))
            paru_update_rel_ind_row (curEl, el, cc) ;
            PRLEVEL (p, ("%%heap element rel col update %ld \n", e));
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
            PRLEVEL (p, ("%%Prior assembly Free %ld  %p size %ld\n",
                        e, el, tot_size));
            elementList[e] = NULL;
            remove_heap (i, lacList, (*curHeap));
            continue; 
        }

        //if (el->rValid >= pMark && elRow[e] == 0)
        //{ // all the rows are inside current fron
        //    // assembling several columns
        //    if (el->cValid >= pMark &&elCol[e] < el->ncolsleft)
        //    {//asslemble all and delete from the heap
        //        paru_update_rel_ind_row (curEl, el, cc) ;
        //        //TODO: use exactly this withoug pop
        //        //std::pop_heap(elHeap->begin()+i, elHeap->end(), greater);
        //    }
        //}
    }

    ////////////////////LOOKING AT HEAP AS A LIST //////////////////////////////
//        ii = 0;
//        for(Int i = 0 ; i < curHeap->size(); i++)
//        {
//            Int e = (*curHeap)[i];
//            paru_Element *el = elementList[e];
//    
//            if ( el == NULL )
//            {
//                continue;
//            }
//            (*curHeap)[ii++] = (*curHeap)[i];
//        }
//    
//        if ( ii < curHeap->size())
//        {
//            PRLEVEL (p, ("%% Prior: size was %ld ", curHeap->size()));
//            PRLEVEL (p, (" and now is %ld\n ", ii));
//            curHeap->resize(ii);
//        }
//
//
//   ///TODO: this line should be deleted after prior is corrected
//    std::make_heap(curHeap->begin(), curHeap->end(), greater ); 


#ifndef NDEBUG  
    PRLEVEL (p, ("%% After assembly current heap:\n %%"));
    for(Int k = 0 ; k < curHeap->size(); k++)
    {
        Int ee = (*curHeap)[k];
        paru_Element *ell = elementList[ee];
        PRLEVEL (p, ("%ld-%ld", k, ee));
        if (ell != NULL)
        {PRLEVEL (p, ("(%ld) ", lacList[ee] ));}
        else
        { PRLEVEL (p, ("(*%ld) ", lacList[ee] ));}

    }
    PRLEVEL (p, ("\n"));

    //chekcing the heap
    for(Int i = curHeap->size()-1 ; i > 0; i--)
    {
        Int elid = (*curHeap)[i];
        Int pelid = (*curHeap)[(i-1)/2]; //parent id
        if ( lacList[pelid] > lacList[elid] )
        {
            PRLEVEL (p, ("%ld-%ld(%ld) <", (i-1)/2, pelid, lacList[pelid] ));
            PRLEVEL (p, ("%ld-%ld(%ld) \n", i, elid, lacList[elid] ));
        }
        ASSERT ( lacList[pelid] <= lacList[elid]);
    }

#endif

    // free the sorting space if allocated
    paru_free ( 2*curElNrows, sizeof(Int), curEl->rWork, cc); 
    curEl->rWork = NULL;
}
