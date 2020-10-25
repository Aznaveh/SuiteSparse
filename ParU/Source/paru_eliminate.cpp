/** =========================================================================  /
 *  ======================  paru_eliminate ==================================  /
 *  ========================================================================= */
/*! @brief  finding the  columns of prior element and fully eliminate it and add
 * it to the current front 
 *
 *  @author Aznaveh
 */
#include "Parallel_LU.hpp"

void paru_eliminate ( Int e, Int f, 
        std::unordered_map <Int, Int> colHash, 
        paru_matrix *paruMatInfo,
        cholmod_common *cc)

{
    DEBUGLEVEL(0);
#ifndef NDEBUG  
    Int p = 1;
#endif

    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM [f]; 


    paru_Element **elementList = paruMatInfo->elementList;

    paru_Element *el = elementList[e];
    paru_Element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;
 
    //Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int*)(el+1);

    //Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;

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
        double *dC = curEl_Num + colHash [colInd]*curEl->nrows;
        PRLEVEL (1, ("%% colHash = %ld \n", colHash [colInd] ));
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
            double *dC = curEl_Num + colHash [colInd]*curEl->nrows;
            PRLEVEL (1, ("%% colHash = %ld \n", colHash [colInd] ));

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
