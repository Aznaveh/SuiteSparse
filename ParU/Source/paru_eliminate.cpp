/** =========================================================================  /
 *  ======================  paru_eliminate ==================================  /
 *  ========================================================================= */
/*! @brief  finding the  columns of prior element and fully eliminate it and add
 * it to the current front 
 *
 *  @author Aznaveh
 */
#include "Parallel_LU.hpp"

void paru_eliminate_all ( Int e, Int f, 
        std::vector <Int> colHash, 
        paru_matrix *paruMatInfo,
        cholmod_common *cc)

{
    DEBUGLEVEL(1);
#ifndef NDEBUG  
    Int p = 1;
#endif

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM [f]; 
    PRLEVEL (p, ("%% Eliminat all of %ld in %ld\n", e, eli));


    paru_Element **elementList = paruMatInfo->elementList;

    paru_Element *el = elementList[e];
    paru_Element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    //Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int*)(el+1);

    //Int *rowRelIndex = relRowInd (el);
    //Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;

    if (el->cValid != paruMatInfo->time_stamp[f] )
        paru_update_rel_ind_col ( e, f, colHash, paruMatInfo) ;

    // Int *colRelIndex = relColInd (paru_Element *el);
    Int *colRelIndex = (Int*)(el+1) + mEl+ nEl;

    //Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int*) (el+1) + nEl; 

    //double *el_Num = numeric_pointer (el);
    double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);
    // current elemnt numerical pointer
    //double *el_Num = numeric_pointer (curEl);
    double *curEl_Num = (double*)((Int*)
            (curEl+1) + 2*curEl->nrowsleft+ 2*curEl->ncolsleft);

    work_struct *Work =  paruMatInfo->Work;
    Int *isRowInFront = Work->rowSize; 
    Int *rowMarkp = Work->rowMark;
    Int rowMark = rowMarkp[eli];

    Int *fcolList = paruMatInfo->fcolList[f];
    paru_fac *Us =  paruMatInfo->partial_Us;
    Int colCount = Us[f].n;

    ASSERT (el_colIndex[el->lac]  <= fcolList[colCount-1] );
    ASSERT (el_colIndex[nEl-1] <= 0 || fcolList[0] <= el_colIndex[nEl-1]);

    PRLEVEL (p, ("%% newColSet.size = %ld\n", colCount ));
    PRLEVEL (p, ("%% nEl = %ld\n",nEl));

    if ( el->ncolsleft == 1)
    {
        PRLEVEL (p, ("%% 1 col left\n %%"));
        double *sC = el_Num + mEl*el->lac; //source column pointer
        Int colInd = el_colIndex [el->lac];
        PRLEVEL (1, ("%% colInd =%ld \n", colInd));
        ASSERT (colInd >= 0);
        // Int fcolcolind = paru_find_hash (colInd, colHash, fcolList);
        Int fcolcolind = colRelIndex [el->lac];
        double *dC = curEl_Num + fcolcolind*curEl->nrows;
        Int nrowsSeen = el->nrowsleft;
        for (Int i = 0; i < mEl ; i++) 
        {
            Int rowInd = el_rowIndex[i];
            PRLEVEL (1, ("%% rowInd =%ld \n", rowInd));
            if (rowInd >= 0 )
            {
                Int ri = isRowInFront [rowInd];
                PRLEVEL (1, ("%% ri = %ld \n", ri));
                PRLEVEL (1, ("%% sC [%ld] =%2.5lf \n", i, sC [i]));
                PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", ri, dC [ri]));
                dC [ri ] += sC[i];
                PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", i, dC [ri]));
                if (-- nrowsSeen == 0) break;
            }
        }

    }
    else
    {
        PRLEVEL (p, ("%% more than 1 col left\n %%"));

        // save the structure of the rows once at first
        Int tempRow[el->nrowsleft]; //C99 
        Int ii = 0;
        for (Int i = 0; i < mEl ; i++) 
        {
            Int rowInd = el_rowIndex[i];
            PRLEVEL (1, ("%% rowInd =%ld \n", rowInd));
            if (rowInd >= 0 )
            {
                tempRow[ii++] = i;
                if (ii == el->nrowsleft) break;
            }
        }


        for (Int j = el->lac; j < nEl ; j++) 
        {

            PRLEVEL (1, ("%% j =%ld \n", j));
            double *sC = el_Num + mEl*j; //source column pointer
            Int colInd = el_colIndex [j];
            PRLEVEL (1, ("%% colInd =%ld \n", colInd));
            if (colInd < 0) continue;
            //Int fcolcolind = paru_find_hash (colInd, colHash, fcolList);
            Int fcolcolind = colRelIndex [j];

            double *dC = curEl_Num + fcolcolind*curEl->nrows;

            for (Int ii = 0; ii < el->nrowsleft; ii++) 
            {
                Int i = tempRow[ii];
                Int rowInd = el_rowIndex[i];
                Int ri = isRowInFront [rowInd];

                PRLEVEL (1, ("%% ri = %ld \n", ri));
                PRLEVEL (1, ("%% sC [%ld] =%2.5lf \n", i, sC [i]));
                PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", ri, dC [ri]));
                dC [ri ] += sC[i];
                PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", i, dC [ri]));

            }

            if (-- el->ncolsleft == 0) break;
            PRLEVEL (1, ("\n"));
        }
    }

    Int tot_size = sizeof(paru_Element) +
        sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
    paru_free (1, tot_size, el, cc);
    PRLEVEL (p, ("%%Eliminate assembly Free %ld  %p size %ld\n",
                e, el, tot_size));
    elementList[e] = NULL;
}

// try to find columns and assemble them to current front. After the first
// column that is not in current front it gets a toll for each column doesn't
// fit

void paru_eliminate_cols ( Int e, Int f, 
        std::vector <Int> colHash, 
        paru_matrix *paruMatInfo,
        cholmod_common *cc)

{

    DEBUGLEVEL(0);
#ifndef NDEBUG  
    Int p = 1;
    Int c = 0; //number of columns assembled
#endif
    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM [f]; 

    PRLEVEL (p, ("%% Eliminat some cols of %ld in %ld\n", e, eli));
#ifndef NDEBUG
    p = 1;

    PRLEVEL (p, ("%% %ld :\n", eli));
    if (p <= 0) paru_print_element (paruMatInfo, eli);

    PRLEVEL (p, ("%% %ld :\n", e));
    if (p <= 0) paru_print_element (paruMatInfo, e);
#endif


    paru_Element **elementList = paruMatInfo->elementList;

    paru_Element *el = elementList[e];
    paru_Element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    //Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int*)(el+1);

    //Int *rowRelIndex = relRowInd (el);
    //Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;

    //Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int*) (el+1) + nEl; 

    //double *el_Num = numeric_pointer (el);
    double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);
    // current elemnt numerical pointer
    //double *el_Num = numeric_pointer (curEl);
    double *curEl_Num = (double*)((Int*)
            (curEl+1) + 2*curEl->nrowsleft+ 2*curEl->ncolsleft);

    work_struct *Work =  paruMatInfo->Work;
    Int *isRowInFront = Work->rowSize; 
    Int *rowMarkp = Work->rowMark;
    Int rowMark = rowMarkp[eli];

    Int *fcolList = paruMatInfo->fcolList[f];
    paru_fac *Us =  paruMatInfo->partial_Us;
    Int colCount = Us[f].n;

    //    ASSERT (el_colIndex[el->lac]  <= fcolList[colCount-1] );
    //    ASSERT (el_colIndex[nEl-1] <= 0 || fcolList[0] <= el_colIndex[nEl-1]);


    Int tempRow[el->nrowsleft]; //C99 
    Int tempRow_ready = 0;
    Int toll = 8; //number of times it continue when do not find anything


    //TOLL FREE zone
    while (paru_find_hash (el_colIndex[el->lac], colHash, fcolList)!= -1 )
    {
        PRLEVEL (p, ("%% Toll free\n"));
        if (tempRow_ready == 0 )
        {
            // save the structure of the rows once at first
            Int ii = 0;
            for (Int i = 0; i < mEl ; i++) 
            {
                Int rowInd = el_rowIndex[i];
                PRLEVEL (1, ("%% rowInd =%ld \n", rowInd));
                if (rowInd >= 0 )
                {
                    tempRow[ii++] = i;
                    if (ii == el->nrowsleft) break;
                }
            } 
            tempRow_ready = 1 ;
        }

        Int colInd = el_colIndex [el->lac];
        Int fcolcolind = paru_find_hash (colInd, colHash, fcolList);

        PRLEVEL (1, ("%% el->lac =%ld \n", el->lac));
        double *sC = el_Num + mEl*el->lac; //source column pointer
        PRLEVEL (1, ("%% colInd =%ld \n", colInd));
        ASSERT (colInd >= 0);

        double *dC = curEl_Num + fcolcolind*curEl->nrows;

        for (Int ii = 0; ii < el->nrowsleft; ii++) 
        {
            Int i = tempRow[ii];
            Int rowInd = el_rowIndex[i];
            Int ri = isRowInFront [rowInd];

            PRLEVEL (1, ("%% ri = %ld \n", ri));
            PRLEVEL (1, ("%% sC [%ld] =%2.5lf \n", i, sC [i]));
            PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", ri, dC [ri]));
            dC [ri ] += sC[i];
            PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", i, dC [ri]));

        }
#ifndef NDEBUG  
        c++;
#endif
        el_colIndex[el->lac] = flip (el_colIndex[el->lac] );
        if (-- el->ncolsleft == 0) break;
        while (el_colIndex[++el->lac] < 0 && el->lac < el->ncols);
    }
    // el->lac won't get updated after this
    Int *lacList = paruMatInfo->lacList;
    lacList[e] = el_colIndex[el->lac];


    //TOLL Zone

    for (Int j = el->lac+1; j < nEl && el->ncolsleft > 0 && toll > 0 ; j++) 
    {
        PRLEVEL (p, ("%% Toll zone\n"));
        toll --;
        if (tempRow_ready == 0 )
        {
            // save the structure of the rows once at first
            Int ii = 0;
            for (Int i = 0; i < mEl ; i++) 
            {
                Int rowInd = el_rowIndex[i];
                PRLEVEL (1, ("%% rowInd =%ld \n", rowInd));
                if (rowInd >= 0 )
                {
                    tempRow[ii++] = i;
                    if (ii == el->nrowsleft) break;
                }
            } 
            tempRow_ready = 1 ;
        }


        PRLEVEL (1, ("%% j =%ld \n", j));
        double *sC = el_Num + mEl*j; //source column pointer
        Int colInd = el_colIndex [j];
        PRLEVEL (1, ("%% colInd =%ld \n", colInd));
        if (colInd < 0) continue;
        Int fcolcolind = paru_find_hash (colInd, colHash, fcolList);
        if (fcolcolind == -1 )
        {// not found
            continue;
        }
        toll++; //if found 
        double *dC = curEl_Num + fcolcolind*curEl->nrows;

        for (Int ii = 0; ii < el->nrowsleft; ii++)
        {
            Int i = tempRow[ii];
            Int rowInd = el_rowIndex[i];
            Int ri = isRowInFront [rowInd];

            PRLEVEL (1, ("%% ri = %ld \n", ri));
            PRLEVEL (1, ("%% sC [%ld] =%2.5lf \n", i, sC [i]));
            PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", ri, dC [ri]));
            dC [ri ] += sC[i];
            PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", i, dC [ri]));

        }
#ifndef NDEBUG  
        c++;
#endif
        el_colIndex[j] = flip (el_colIndex[j] );
        if (-- el->ncolsleft == 0) break;
    }

    PRLEVEL (1, ("%%  %ld has found and assembled, ncolsleft %ld\n", 
                c, el->ncolsleft));

    if (el->ncolsleft == 0 )
    {
        Int tot_size = sizeof(paru_Element) +
            sizeof(Int)*(2*(mEl+nEl)) + sizeof(double)*nEl*mEl;
        paru_free (1, tot_size, el, cc);
        PRLEVEL (p, ("%%Some cols Free ALL %ld  %p size %ld\n",
                    e, el, tot_size));
        elementList[e] = NULL;
    }
}
