/** =========================================================================  /
 *  ======================  paru_eliminate ==================================  /
 *  ========================================================================= */
/*! @brief  finding the  columns of prior element and eliminate it to the
 * current front 
 *
 *  @author Aznaveh
 */
#include "Parallel_LU.hpp"
#define C 4

void paru_eliminate ( Int e, Int f, 
        std::unordered_map <Int, Int> colHash, 
        paru_matrix *paruMatInfo)

{
    DEBUGLEVEL(1);
#ifndef NDEBUG  
    Int p = 1;
#endif


    //TODO bring the col list
    paru_symbolic *LUsym =  paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM [f]; 


    paru_Element **elementList = paruMatInfo->elementList;
    paru_Element *el = elementList[e];

    paru_Element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;
 
    Int *el_colIndex = (Int*)(el+1);

    //Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int*)(el+1) + 2*nEl +mEl;

    //Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int*) (el+1) + nEl; 

    //double *el_Num = numeric_pointer (el);
    double *el_Num = (double*)((Int*)(el+1) + 2*nEl+ 2*mEl);
    // current elemnt numerical pointer
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

//    if ( el->ncolsleft == 1)
//    {
//        PRLEVEL (p, ("%% 1 col left\n %%"));
//        //TODO linear time go throuh el->lac and assemble it and done
//        double *sC = el_Num + el->lac*mEl; //source column pointer
//        Int colind = el_colIndex [el->lac];
//        ASSERT (colind != LONG_MAX);
//        ASSERT (colind >= 0);
//        double *dC = curEl_Num + colHash [colind]*curEl->nrows;
//        Int nrowsSeen = el->nrowsleft;
//        for (Int i = 0; i < mEl ; i++) 
//        {
//            Int rowIndex = el_rowIndex[i];
//            Int ri = isRowInFront [ri];
//            if (ri >= 0 )
//            {
//                PRLEVEL (1, ("%% sC [%ld] =%2.5lf \n", i, sC [i]));
//                PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", i, dC [ri]));
//                dC [ri ] += sC[i];
//                PRLEVEL (1, ("%% dC [%ld] =%2.5lf \n", i, dC [ri]));
//                if (-- nrowsSeen == 0)
//                    break;
//            }
//        }
//        PRLEVEL (1, ("\n"));
//    }
//    else
    {
        PRLEVEL (p, ("%% more than 1 col left\n %%"));
        for (Int j = el->lac; j < nEl ; j++) 
        {
            PRLEVEL (1, ("%% j =%ld \n", j));
            double *sC = el_Num + mEl*j; //source column pointer
            Int colInd = el_colIndex [j];
            PRLEVEL (1, ("%% colInd =%ld \n", colInd));
            if (colInd < 0) continue;
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
            if (-- el->ncolsleft == 0) break;
            PRLEVEL (1, ("\n"));
        }
    }
}
