/** =========================================================================  /
 * =======================  paru_intersection ================================  /
 * ========================================================================== */
/*! @brief  finding the number of intersection between new columns and piror
 * block columns
 *
 *  @author Aznaveh
 */
#include "Parallel_LU.hpp"
#define C 4

int paru_intersection ( Int e, paru_Element **elementList, 
        std::set<Int> &stl_newColSet)
{
    DEBUGLEVEL(1);
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
    if (el_colIndex[0]  > *stl_newColSet.end())
        return 0;
    if (*stl_newColSet.begin() > el_colIndex[nEl-1])
        return 0;

    PRLEVEL (p, ("%% newColSet.size = %ld\n",stl_newColSet.size() ));
    PRLEVEL (p, ("%% nEl = %ld\n",nEl));
    std::set<Int>::iterator it;
    if (C*stl_newColSet.size() < nEl )
        // if size el >> stl_newColSet 
        //   binary search each of elements in stl_newColSet in el
        //   log(nEl)*stl_newColSet.size()
    {
        PRLEVEL (p, ("%%el >> stl_newColSet\n"));
        for (it = stl_newColSet.begin(); 
                it != stl_newColSet.end(); it++)
        {
            Int c = *it;
            Int col = bin_srch_col (el_colIndex, el->lnc, nEl, c);
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
    else if (stl_newColSet.size() > C*nEl )
    {
        PRLEVEL (p, ("%%el << stl_newColSet\n"));
        //  else if stl_newColSet >> el
        //      binary search each of elements in el in stl_newColSet 
        //      log(stl_newColSet.size())*nEl
        for (Int c = el->lnc; c < nEl; c++)
        {
            Int col = el_colIndex[c];
            if (col < 0) continue;
            if ( stl_newColSet.find(col) != stl_newColSet.end() )
            {
                intersection++; 
                PRLEVEL (p, ("%%2: col=%ld ", col));
                PRLEVEL (p, ("%%intersection=%ld\n", intersection));
            };
        }
    }
    else
    { //Merge style m+n
        
        
        PRLEVEL (p, ("%%Merge style\n"));
        it = stl_newColSet.begin(); 
        Int c = el->lnc;
        while (it != stl_newColSet.end() && c < nEl)
        { 
            while (el_colIndex[c] < 0) ++c; //skip dead columns
            if (c >= nEl)
                break;

            if (*it < el_colIndex[c])
                it++;
            else if (  el_colIndex[c] < *it)
                c++; 
            else
                // *it == col
            { 
                intersection++; 
                PRLEVEL (p, ("%%3: c=%ld ", c));
                PRLEVEL (p, ("%%col= %ld", el_colIndex[c]));
                PRLEVEL (p, ("%%intersection=%ld\n", intersection));
                it++; c++; 
            };
        }
    }
    PRLEVEL (p, (" e = %ld intersection= %ld\n",e, intersection));
    return intersection;
}
