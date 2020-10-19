/** =========================================================================  /
 * =======================  paru_eliminate ===================================  /
 * ========================================================================== */
/*! @brief  finding the  columns of prior element and eliminate it to the
 * current front 
 *
 *  @author Aznaveh
 */
#include "Parallel_LU.hpp"
#define C 4

void paru_eliminate ( Int e, Int f, 
        std::unordered_map <Int, Int> rowHash, 
        std::unordered_map <Int, Int> colHash, 
        paru_matrix *paruMatInfo)

{
    DEBUGLEVEL(0);

    //TODO bring the col list
    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM [f]; 


    paru_Element **elementList = paruMatInfo->elementList;
    paru_Element *el = elementList[e];

    Int nEl = el->ncols;
    Int *el_colIndex = (Int*)(el+1);

    //Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;

    //Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int*) (el+1) + nEl; 

    //double *el_Num = numeric_pointer (el);
    double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);
    // current elemnt numerical pointer
    double *curEl_Num = (double*)((Int*)(eli+1) + 2*nEl+ 2*mEl);

    //find the intersection between columns of e and 
    Int intersection = 0;
    // conditions for early stop

    ASSERT (el_colIndex[el->lac]  <= *stl_colSet.end())
        ASSERT (el_colIndex[nEl-1] <= 0 || 
                *stl_colSet.begin() <= el_colIndex[nEl-1])

        PRLEVEL (p, ("%% newColSet.size = %ld\n",stl_colSet.size() ));
    PRLEVEL (p, ("%% nEl = %ld\n",nEl));
    std::set<Int>::iterator it;

    if ( el->ncolsleft == 1)
    {
        //TODO linear time go throuh el->lac and assemble it and done
        Int sc = el->lac; //source column
        Int colind = el_colIndex [sc];
        Int dc = colHash [colind];
        for (Int i = 0; i < el->nrowsleft; i++) 
        {
            Int rowIndi = el_rowIndex[i];
            Int ri = rowHash[rowIndi];
            if (ri >= 0 )
            {
            }
        }
    }
    else
    {
        //TODO traverse once throuh row indices update relative indices and save
        //the active rows then just find active columns and assemble them one
        //bye one
        Int tempRel[el->nrowsleft]; 
        // C99 finding the row indices fisrt 

        //finding the active row indices
        Int jj = 0;
        for (Int j = 0; j < el->nrows; j++)
        {
            if (el_rowIndex[j] > 0)
            {
                tempRel[j] = jj++; 
                if (jj == el->nrowsleft) break;
            }
        }

        ASSERT (C*stl_colSet.size() >= el->ncolsleft);

        //        if (C*stl_colSet.size() < el->ncolsleft) //TODO this never happens here
        //            // if size el >> stl_colSet 
        //            //   binary search each of elements in stl_colSet in el
        //            //   log(nEl)*stl_colSet.size()
        //        {
        //            PRLEVEL (p, ("%%el >> stl_colSet\n"));
        //            for (it = stl_colSet.begin(); it != stl_colSet.end(); it++)
        //            {
        //                Int c = *it;
        //                Int col = bin_srch_col (el_colIndex, el->lac, nEl, c);
        //                PRLEVEL (p, ("%%intersection=%ld", intersection));
        //                PRLEVEL (p, ("%%after bsearch for c=%ld col=%ld \n",c, col ));
        //                if ( col != -1 && el_colIndex[col] == c) 
        //                { 
        //                    intersection++;
        //                    PRLEVEL (p, ("%%##1: c=%ld ", c));
        //                    PRLEVEL (p, ("%%intersection=%ld\n", intersection));
        //                };
        //            }
        //        }
        if (stl_colSet.size() > C*(el->ncolsleft) )
        {   //  else if stl_colSet >> el
            PRLEVEL (p, ("%%el << stl_colSet\n"));
            //      binary search each of elements in el in stl_colSet 
            //      log(stl_colSet.size())*nEl
            Int ncolsseen = el->ncolsleft;
            for (Int c = el->lac; c < nEl; c++)
            {
                Int col = el_colIndex[c];
                if (col < 0) continue;
                ncolsseen--;
                if ( stl_colSet.find(col) != stl_colSet.end() )
                {
                    intersection++; 
                    PRLEVEL (p, ("%%2: col=%ld ", col));
                    PRLEVEL (p, ("%%intersection=%ld\n", intersection));
                };
                if (ncolsseen == 0) break;
            }
        }
        else
        { //Merge style m+n
            PRLEVEL (p, ("%%Merge style\n"));
            it = stl_colSet.begin(); 
            Int c = el->lac;
            Int ncolsseen = el->ncolsleft;
            while (it != stl_colSet.end() && ncolsseen > 0)
            { 
                while (el_colIndex[c] < 0 && c < nEl) ++c; //skip dead columns
                if (c >= nEl)
                    break;

                if (*it < el_colIndex[c])
                    it++;
                else if (  el_colIndex[c] < *it)
                {
                    c++; 
                    ncolsseen--;
                }
                else
                    // *it == col
                { 
                    intersection++; 
                    PRLEVEL (p, ("%%3: c=%ld ", c));
                    PRLEVEL (p, ("%%col= %ld", el_colIndex[c]));
                    PRLEVEL (p, ("%%intersection=%ld\n", intersection));
                    it++; c++; 
                    ncolsseen--;
                };
            }
        }
    }
    PRLEVEL (p, (" e = %ld intersection= %ld\n",e, intersection));
}

