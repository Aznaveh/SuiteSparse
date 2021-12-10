////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_full_summed ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  fully sum the pivotal column from a prior contribution block e
 *  assembling starts with el->lac and continues until no assembly is possible
 *  in curent front. If the prior contribution block is empty free it here.
 *
 *
 ***************  assembling the pivotal part of the front 
 *
 *      el           nEl
 *                  6, 7, 11, 12
 *                 _____________
 *              23 | X  Y  .  .     stored in memory like this:
 *          mEl 17 | X  Y  .  .     ..6, 7,11, 12, 23, 17, 2, X, X, X, Y, Y, Y,
 *               2 | X  Y  .  .
 *
 *        It must be assembled in current pivotal fron like this:
 *                                         fp
 *                                     col1, ... , col
 *
 *                                      6, 7, 8, 9, 10
 *                                      ______________
 *                              0   23 | X  Y  .  .  .
 *                   rowCount   1    2 | X  Y  .  .  .
 *                              2    4 | *  *  .  .  .  isRowInFront[4] == 2
 *                              3   17 | X  Y  .  .  .
 *
 *  @author Aznaveh
 */
#include "paru_internal.hpp"

void paru_full_summed(Int e, Int f, paru_matrix *paruMatInfo)

{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    Int p = 1;
#endif
    paru_symbolic *LUsym = paruMatInfo->LUsym;
#ifndef NDEBUG
    Int *snM = LUsym->super2atree;
    Int eli = snM[f];
    PRLEVEL(p, ("%% Fully summing %ld in %ld(%ld)\n", e, f, eli));
#endif

    Int *Super = LUsym->Super;
    Int col1 = Super[f]; /* fornt F has columns col1:col2-1 */
    Int col2 = Super[f + 1];
    PRLEVEL(p, ("%% col1=%ld, col2=%ld\n", col1, col2));

    paru_Element **elementList = paruMatInfo->elementList;

    paru_Element *el = elementList[e];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    // Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

    // Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int *)(el + 1) + nEl;

    paru_fac *LUs = paruMatInfo->partial_LUs;
    Int rowCount = LUs[f].m;
    double *pivotalFront = LUs[f].p;

    // double *el_Num = numeric_pointer (el);
    double *el_Num = (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);

#ifndef NDEBUG  // print the element which is going to be assembled from
    p = 2;
    PRLEVEL(p, ("%% ASSEMBL element= %ld  mEl =%ld ", e, mEl));
    if (p <= 0) paru_print_element(paruMatInfo, e);
#endif

    Int j = el->lac;  // keep record of latest lac
    if (el->ncolsleft == 1)
    // I could add the following condition too:
    //|| next active column is out of current pivotal columns
    // However it requires searching through columns; I will just waste a
    // temp space
    {  // No need for a temp space for active rows; no reuse
        PRLEVEL(p, ("%% 1 col left\n %%"));
        double *sC = el_Num + mEl * el->lac;  // source column pointer
        Int fcolInd = el_colIndex[el->lac] - col1;
#ifndef NDEBUG
        Int colInd = el_colIndex[el->lac];
        PRLEVEL(1, ("%% colInd =%ld \n", fcolInd));
        ASSERT(colInd >= 0);
#endif
        double *dC = pivotalFront + fcolInd * rowCount;
        Int nrows2bSeen = el->nrowsleft;
        for (Int i = 0; i < mEl; i++)
        {
            Int rowInd = el_rowIndex[i];
            PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
            if (rowInd >= 0 && rowRelIndex[i] != -1)
            {  // active and do not contain zero in pivot
                Int ri = rowRelIndex[i];
                PRLEVEL(1, ("%% ri = %ld \n", ri));
                PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
                dC[ri] += sC[i];
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
                el_colIndex[el->lac] = flip(el_colIndex[el->lac]);
                if (--nrows2bSeen == 0) break;
            }
        }
        el->ncolsleft--;
    }
    else
    {
        PRLEVEL(p, ("%% more than 1 col left\n %%"));

        // save the structure of the rows once at first
        Int nrows2assembl = el->nrowsleft - el->nzr_pc;
        // Int tempRow[nrows2assembl];  // C99
        std::vector<Int> tempRow(nrows2assembl);
        Int ii = 0;
        for (Int i = 0; i < mEl; i++)
        {
            Int rowInd = el_rowIndex[i];
            PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
            if (rowInd >= 0 && rowRelIndex[i] != -1)
            {
                tempRow[ii++] = i;
                if (ii == nrows2assembl) break;
            }
        }

#ifndef NDEBUG
        p = 1;
        PRLEVEL(p, ("%% list of the rows to be assembled:\n%%"));
        for (Int i = 0; i < nrows2assembl; i++)
            PRLEVEL(p, ("%ld ", el_rowIndex[tempRow[i]]));
        PRLEVEL(p, ("%% \n"));
#endif
        Int *Depth = LUsym->Depth;
        #pragma omp parallel 
        #pragma omp single
        #pragma omp taskgroup
        for (; j < nEl; j++)
        {  // j already defined out of this scope while it is needed
            PRLEVEL(1, ("%% j =%ld \n", j));
            double *sC = el_Num + mEl * j;  // source column pointer
            Int colInd = el_colIndex[j];
            PRLEVEL(1, ("%% colInd =%ld \n", colInd));
            if (colInd >= col2) break;
            if (colInd < 0) continue;

            Int fcolInd = colInd - col1;

            double *dC = pivotalFront + fcolInd * rowCount;

            #pragma omp task priority(Depth[f]) if(nrows2assembl > TASK_MIN)
            for (Int iii = 0; iii < nrows2assembl; iii++)
            {
                Int i = tempRow[iii];
                Int ri = rowRelIndex[i];

                ASSERT(rowRelIndex[i] != -1);  // I already picked the rows
                // that are not in zero pivots
                ASSERT(el_rowIndex[i] >= 0);  // and also still alive

                PRLEVEL(1, ("%% ri = %ld \n", ri));
                PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
                dC[ri] += sC[i];
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
            }

            el_colIndex[j] = flip(el_colIndex[j]);
            if (--el->ncolsleft == 0) break;
            PRLEVEL(1, ("\n"));
        }
    }

    if (el->ncolsleft == 0)
    {  // free el
        PRLEVEL(p, ("%% element %ld is freed after pivotal assembly\n", e));
        paru_free_el(e, elementList);
    }

    if (elementList[e] != NULL)
    {
        el->lac = j;
        Int *lacList = paruMatInfo->lacList;
        lacList[e] = lac_el(elementList, e);
        PRLEVEL(1, ("%%e = %ld, el->lac= %ld ", e, el->lac));
        PRLEVEL(1, ("el_colIndex[%ld]=%ld :\n", el->lac, el_colIndex[el->lac]));
        ASSERT(j < nEl);
    }
#ifndef NDEBUG  // print the element which has been assembled from
    p = 1;
    PRLEVEL(p, ("%% ASSEMBLED element= %ld  mEl =%ld ", e, mEl));
    if (p <= 0) paru_print_element(paruMatInfo, e);

    // Printing the pivotal front
    p = 2;
    PRLEVEL(p, ("%% After Assemble element %ld\n", e));
    PRLEVEL(p, ("%% x =  \t"));
    for (Int c = col1; c < col2; c++) PRLEVEL(p, ("%ld\t\t", c));
    PRLEVEL(p, (" ;\n"));

    Int *frowList = paruMatInfo->frowList[f];
    for (Int r = 0; r < rowCount; r++)
    {
        PRLEVEL(p, ("%% %ld\t", frowList[r]));
        for (Int c = col1; c < col2; c++)
            PRLEVEL(p, (" %2.5lf\t", pivotalFront[(c - col1) * rowCount + r]));
        PRLEVEL(p, ("\n"));
    }
    p = 1;
#endif
}
