/** =========================================================================  /
 * =======================  paru_intersection ================================  /
 * ========================================================================== */
/*! @brief  finding the number of intersection between new columns and piror
 * block columns
 *
 *  @author Aznaveh
 */
#include "paru_internal.hpp"
#define C 4

int paru_intersection ( Int e, paru_Element **elementList, 
        std::set<Int> &stl_newColSet)
{
    DEBUGLEVEL(0);
#ifndef NDEBUG        
    Int p = 1;
    PRLEVEL (p, ("%%stl_newColSet:\n%%"));
    for (std::set<Int>::iterator it = stl_newColSet.begin(); 
            it != stl_newColSet.end(); it++)
    {
        PRLEVEL (p, (" %ld",*it));
    }
    PRLEVEL (p, ("\n"));

#endif

    paru_Element *el = elementList[e];

    Int nEl = el->ncols;
    Int *el_colIndex = (Int*)(el+1);

    //find the intersection between columns of e and 
    Int intersection = 0;
    // conditions for early stop
    if (el_colIndex[el->lac]  > *stl_newColSet.end())
        return 0;
    if (el_colIndex[nEl-1] > 0 && *stl_newColSet.begin() > el_colIndex[nEl-1])
        return 0;

    PRLEVEL (p, ("%% newColSet.size = %ld\n",stl_newColSet.size() ));
    PRLEVEL (p, ("%% nEl = %ld\n",nEl));
    std::set<Int>::iterator it;
    if ( (Int) (C*stl_newColSet.size()) < nEl-el->lac )
        // if size el >> stl_newColSet 
        //   binary search each of elements in stl_newColSet in el
        //   log(nEl)*stl_newColSet.size()
    {
        PRLEVEL (p, ("%%el >> stl_newColSet\n"));
        for (it = stl_newColSet.begin(); it != stl_newColSet.end(); it++)
        {
            Int c = *it;
            Int col = bin_srch_col (el_colIndex, el->lac, nEl-1, c);
            PRLEVEL (p, ("%%intersection=%ld", intersection));
            PRLEVEL (p, ("%%after bsearch for c=%ld col=%ld \n",c, col ));
            if ( col != -1 && el_colIndex[col] == c) 
            { 
                intersection++;
                PRLEVEL (p, ("%%##1: c=%ld ", c));
                PRLEVEL (p, ("%%intersection=%ld\n", intersection));
            };
        }
    }
    else if ( (Int) stl_newColSet.size() > C*(nEl-el->lac) )
    {   //  else if stl_newColSet >> el
        PRLEVEL (p, ("%%el << stl_newColSet\n"));
        //      binary search each of elements in el in stl_newColSet 
        //      log(stl_newColSet.size())*nEl
        Int ncolsseen = el->ncolsleft;
        for (Int c = el->lac; c < nEl; c++)
        {
            Int col = el_colIndex[c];
            if (col < 0) continue;
            ncolsseen--;
            if ( stl_newColSet.find(col) != stl_newColSet.end() )
            {
                intersection++; 
                PRLEVEL (p, ("%%2: col=%ld ", col));
                PRLEVEL (p, ("%%intersection=%ld\n", intersection));
            };
            if (ncolsseen == 0) return intersection;
        }
    }
    else
    { //Merge style m+n
        PRLEVEL (p, ("%%Merge style\n"));
        it = stl_newColSet.begin(); 
        Int c = el->lac;
        Int ncolsseen = el->ncolsleft;
        while (it != stl_newColSet.end() && ncolsseen > 0)
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
    PRLEVEL (p, (" e = %ld intersection= %ld\n",e, intersection));
    return intersection;
}
