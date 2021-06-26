/** =========================================================================  /
 *  ======================  paru_assemble ===================================  /
 *  ========================================================================= */
/*! @brief  finding the  columns or rows of prior element and fully or partially
 * assemble it  and eliminate it if needed
 *
 *  @author Aznaveh
 */
#include "paru_internal.hpp"

void paru_assemble_all(Int e, Int f, std::vector<Int> &colHash,
                       paru_matrix *paruMatInfo)

{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    Int p = 1;
#endif

    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM[f];
    PRLEVEL(p, ("%% Eliminat all of %ld in %ld\n", e, eli));

    paru_Element **elementList = paruMatInfo->elementList;

    paru_Element *el = elementList[e];
    paru_Element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    // Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

    if (el->cValid != paruMatInfo->time_stamp[f])
        paru_update_rel_ind_col(e, f, colHash, paruMatInfo);

    // Int *colRelIndex = relColInd (paru_Element *el);
    Int *colRelIndex = (Int *)(el + 1) + mEl + nEl;

    // Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int *)(el + 1) + nEl;

    // double *el_Num = numeric_pointer (el);
    double *el_Num = (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);
    // current elemnt numerical pointer
    // double *el_Num = numeric_pointer (curEl);
    double *curEl_Num =
        (double *)((Int *)(curEl + 1) + 2 * curEl->nrows + 2 * curEl->ncols);

    work_struct *Work = paruMatInfo->Work;
    Int *isRowInFront = Work->rowSize;

#ifndef NDEBUG
    paru_fac *Us = paruMatInfo->partial_Us;
    Int *fcolList = paruMatInfo->fcolList[f];
    Int colCount = Us[f].n;
    ASSERT(el_colIndex[el->lac] <= fcolList[colCount - 1]);
    ASSERT(el_colIndex[nEl - 1] <= 0 || fcolList[0] <= el_colIndex[nEl - 1]);
#endif

    PRLEVEL(p, ("%% newColSet.size = %ld\n", colCount));
    PRLEVEL(p, ("%% nEl = %ld\n", nEl));

    if (el->ncolsleft == 1)
    {
        PRLEVEL(p, ("%% 1 col left\n %%"));
        double *sC = el_Num + mEl * el->lac;  // source column pointer
#ifndef NDEBUG
        Int colInd = el_colIndex[el->lac];
        PRLEVEL(1, ("%% colInd =%ld \n", colInd));
        ASSERT(colInd >= 0);
#endif
        Int fcolind = colRelIndex[el->lac];
        double *dC = curEl_Num + fcolind * curEl->nrows;
        Int nrowsSeen = el->nrowsleft;
        for (Int i = 0; i < mEl; i++)
        {
            Int rowInd = el_rowIndex[i];
            PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
            if (rowInd >= 0)
            {
                Int ri = isRowInFront[rowInd];
                PRLEVEL(1, ("%% ri = %ld \n", ri));
                PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
                dC[ri] += sC[i];
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
                if (--nrowsSeen == 0) break;
            }
        }
    }
    else
    {
        PRLEVEL(p, ("%% more than 1 col left\n %%"));

        // save the structure of the rows once at first
        // Int tempRow[el->nrowsleft];  // C99
        std::vector<Int> tempRow(el->nrowsleft);
        Int ii = 0;
        for (Int i = 0; i < mEl; i++)
        {
            Int rowInd = el_rowIndex[i];
            PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
            if (rowInd >= 0)
            {
                tempRow[ii++] = i;
                Int ri = isRowInFront[rowInd];
                rowRelIndex[i] = ri;
                if (ii == el->nrowsleft) break;
            }
        }

        for (Int j = el->lac; j < nEl; j++)
        {
            PRLEVEL(1, ("%% j =%ld \n", j));
            double *sC = el_Num + mEl * j;  // source column pointer
            Int colInd = el_colIndex[j];
            PRLEVEL(1, ("%% colInd =%ld \n", colInd));
            if (colInd < 0) continue;
            Int fcolind = colRelIndex[j];

            double *dC = curEl_Num + fcolind * curEl->nrows;

            for (Int iii = 0; iii < el->nrowsleft; iii++)
            {
                Int i = tempRow[iii];
                Int ri = rowRelIndex[i];

                PRLEVEL(1, ("%% ri = %ld \n", ri));
                PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
                dC[ri] += sC[i];
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
            }

            if (--el->ncolsleft == 0) break;
            PRLEVEL(1, ("\n"));
        }
    }

    paru_free_el(e, elementList);
}

// try to find columns and assemble them to current front. After the first
// column that is not in current front it gets a toll for each column doesn't
// fit

void paru_assemble_cols(Int e, Int f, std::vector<Int> &colHash,
                        paru_matrix *paruMatInfo)

{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    Int p = 1;
    Int c = 0;  // number of columns assembled
#endif
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM[f];

    PRLEVEL(p, ("%% Eliminat some cols of %ld in %ld\n", e, eli));
#ifndef NDEBUG
    p = 1;

    PRLEVEL(p, ("%% %ld :\n", eli));
    if (p <= 0) paru_print_element(paruMatInfo, eli);

    PRLEVEL(p, ("%% %ld :\n", e));
    if (p <= 0) paru_print_element(paruMatInfo, e);
#endif

    paru_Element **elementList = paruMatInfo->elementList;

    paru_Element *el = elementList[e];
    paru_Element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    // Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

    // Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int *)(el + 1) + nEl;

    // double *el_Num = numeric_pointer (el);
    double *el_Num = (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);
    // current elemnt numerical pointer
    // double *el_Num = numeric_pointer (curEl);
    double *curEl_Num =
        (double *)((Int *)(curEl + 1) + 2 * curEl->nrows + 2 * curEl->ncols);

    work_struct *Work = paruMatInfo->Work;
    Int *isRowInFront = Work->rowSize;

    Int *fcolList = paruMatInfo->fcolList[f];

    // Int tempRow[el->nrowsleft];  // C99
    std::vector<Int> tempRow(el->nrowsleft);
    Int tempRow_ready = 0;
    Int toll = 8;  // number of times it continue when do not find anything

    // TOLL FREE zone
    while (paru_find_hash(el_colIndex[el->lac], colHash, fcolList) != -1)
    {
        PRLEVEL(p, ("%% Toll free\n"));
        if (tempRow_ready == 0)
        {
            // save the structure of the rows once at first
            Int ii = 0;
            for (Int i = 0; i < mEl; i++)
            {
                Int rowInd = el_rowIndex[i];
                PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
                if (rowInd >= 0)
                {
                    tempRow[ii++] = i;
                    Int ri = isRowInFront[rowInd];
                    rowRelIndex[i] = ri;
                    if (ii == el->nrowsleft) break;
                }
            }
            tempRow_ready = 1;
        }

        Int colInd = el_colIndex[el->lac];
        Int fcolind = paru_find_hash(colInd, colHash, fcolList);

        PRLEVEL(1, ("%% el->lac =%ld \n", el->lac));
        double *sC = el_Num + mEl * el->lac;  // source column pointer
        PRLEVEL(1, ("%% colInd =%ld \n", colInd));
        ASSERT(colInd >= 0);

        double *dC = curEl_Num + fcolind * curEl->nrows;

        for (Int ii = 0; ii < el->nrowsleft; ii++)
        {
            Int i = tempRow[ii];
            Int ri = rowRelIndex[i];

            PRLEVEL(1, ("%% ri = %ld \n", ri));
            PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
            PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
            dC[ri] += sC[i];
            PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
        }
#ifndef NDEBUG
        c++;
#endif
        el_colIndex[el->lac] = flip(el_colIndex[el->lac]);
        if (--el->ncolsleft == 0) break;
        while (el_colIndex[++el->lac] < 0 && el->lac < el->ncols)
            ;
    }
    // el->lac won't get updated after this
    Int *lacList = paruMatInfo->lacList;
    lacList[e] = el_colIndex[el->lac];

    // TOLL Zone

    for (Int j = el->lac + 1; j < nEl && el->ncolsleft > 0 && toll > 0; j++)
    {
        PRLEVEL(p, ("%% Toll zone\n"));
        toll--;
        if (tempRow_ready == 0)
        {
            // save the structure of the rows once at first
            Int ii = 0;
            for (Int i = 0; i < mEl; i++)
            {
                Int rowInd = el_rowIndex[i];
                PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
                if (rowInd >= 0)
                {
                    tempRow[ii++] = i;
                    Int ri = isRowInFront[rowInd];
                    rowRelIndex[i] = ri;
                    if (ii == el->nrowsleft) break;
                }
            }
            tempRow_ready = 1;
        }

        PRLEVEL(1, ("%% j =%ld \n", j));
        double *sC = el_Num + mEl * j;  // source column pointer
        Int colInd = el_colIndex[j];
        PRLEVEL(1, ("%% colInd =%ld \n", colInd));
        if (colInd < 0) continue;
        Int fcolind = paru_find_hash(colInd, colHash, fcolList);
        if (fcolind == -1) continue;
        toll++;  // if found
        double *dC = curEl_Num + fcolind * curEl->nrows;

        for (Int ii = 0; ii < el->nrowsleft; ii++)
        {
            Int i = tempRow[ii];
            Int ri = rowRelIndex[i];

            PRLEVEL(1, ("%% ri = %ld \n", ri));
            PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
            PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
            dC[ri] += sC[i];
            PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
        }
#ifndef NDEBUG
        c++;
#endif
        el_colIndex[j] = flip(el_colIndex[j]);
        if (--el->ncolsleft == 0) break;
    }

    PRLEVEL(1, ("%%  %ld has found and assembled, ncolsleft %ld\n", c,
                el->ncolsleft));

    if (el->ncolsleft == 0)
    {
        paru_free_el(e, elementList);
    }
}

void paru_assemble_rows(Int e, Int f, std::vector<Int> &colHash,
                        paru_matrix *paruMatInfo)

{
    DEBUGLEVEL(0);
#ifndef NDEBUG
    Int p = 1;
#endif

    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM[f];

    PRLEVEL(p, ("%% Eliminat some rows of %ld in %ld\n", e, eli));

    paru_Element **elementList = paruMatInfo->elementList;

    paru_Element *el = elementList[e];
    paru_Element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    // Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

    // Int *colRelIndex = relColInd (paru_Element *el);
    Int *colRelIndex = (Int *)(el + 1) + mEl + nEl;

    // Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int *)(el + 1) + nEl;

    // Int *el_rowIndex = rowIndex_pointer (curEl);
    Int *curEl_rowIndex = (Int *)(curEl + 1) + curEl->ncols;

    // double *el_Num = numeric_pointer (el);
    double *el_Num = (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);
    // current elemnt numerical pointer
    // double *el_Num = numeric_pointer (curEl);
    double *curEl_Num =
        (double *)((Int *)(curEl + 1) + 2 * curEl->nrows + 2 * curEl->ncols);

    work_struct *Work = paruMatInfo->Work;
    Int *isRowInFront = Work->rowSize;

    std::vector<Int> tempRow;

    // searching for rows
    Int i = 0;
    Int nrowsSeen = el->nrowsleft;
    // Toll free zone
    PRLEVEL(1, ("%% Toll free\n"));
    while (i < mEl && nrowsSeen > 0)
    {
        for (; el_rowIndex[i] < 0; i++)
            ;
        nrowsSeen--;

        Int rowInd = isRowInFront[i];
        if (rowInd > 0 && rowInd < curEl->nrows)
        {
            // coompare their global indices
            if (curEl_rowIndex[rowInd] == el_rowIndex[i])
            {
                PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
                PRLEVEL(1, ("%% curEl_rowIndex[rowInd] =%ld \n",
                            curEl_rowIndex[rowInd]));
                PRLEVEL(1, ("%% i =%ld \n", i));
                PRLEVEL(1, ("%% el_rowIndex[i] =%ld \n", el_rowIndex[i]));
                tempRow.push_back(i);
            }
            else
                break;
        }
        i++;
    }

#ifndef NDEBUG
    if (tempRow.size() > 0)
        PRLEVEL(p, ("%% Toll free zone: %ld rows has been found: \n%%",
                    tempRow.size()));
#endif

    PRLEVEL(1, ("%% TollED \n"));
    Int toll = 8;  // number of times it continue when do not find anything
    // Toll zone
    while (i < mEl && nrowsSeen > 0 && toll > 0)
    // while (i < mEl  && nrowsSeen >0 )
    {
        for (; el_rowIndex[i] < 0; i++)
            ;
        nrowsSeen--;

        Int rowInd = isRowInFront[i];
        if (rowInd > 0 && rowInd < curEl->nrows)
        {
            // coompare their global indices
            if (curEl_rowIndex[rowInd] == el_rowIndex[i])
            {
                PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
                PRLEVEL(1, ("%% curEl_rowIndex[rowInd] =%ld \n",
                            curEl_rowIndex[rowInd]));
                PRLEVEL(1, ("%% i =%ld \n", i));
                PRLEVEL(1, ("%% el_rowIndex[i] =%ld \n", el_rowIndex[i]));

                tempRow.push_back(i);
                toll++;
            }
            else
                toll--;
        }
        i++;
    }

    if (tempRow.empty()) return;

    PRLEVEL(p,
            ("%% %ld rows has been found, toll %ld\n%%", tempRow.size(), toll));
#ifndef NDEBUG
    for (Int ii = 0; ii < (Int)tempRow.size(); ii++)
        PRLEVEL(p, ("%ld ", tempRow[ii]));
    PRLEVEL(p, ("\n "));
#endif
#ifndef NDEBUG
    p = 1;
    PRLEVEL(p, ("%% Before eliminiatine some rows %ld :\n", eli));
    if (p <= 0) paru_print_element(paruMatInfo, eli);

    PRLEVEL(p, ("%% %ld :\n", e));
    if (p <= 0) paru_print_element(paruMatInfo, e);
#endif

    if (el->cValid != paruMatInfo->time_stamp[f])
        paru_update_rel_ind_col(e, f, colHash, paruMatInfo);

    Int ncolsSeen = nEl;

    for (Int j = el->lac; j < nEl; j++)
    {
        PRLEVEL(1, ("%% j =%ld \n", j));
        double *sC = el_Num + mEl * j;  // source column pointer
        Int colInd = el_colIndex[j];
        PRLEVEL(1, ("%% colInd =%ld \n", colInd));
        if (colInd < 0) continue;
        ncolsSeen--;
        Int fcolind = colRelIndex[j];

        PRLEVEL(1, ("%% fcolind=%ld \n", fcolind));
        double *dC = curEl_Num + fcolind * curEl->nrows;

        for (Int ii = 0; ii < (Int)tempRow.size(); ii++)
        {
            Int i1 = tempRow[ii];
            Int rowInd = el_rowIndex[i1];
            Int ri = isRowInFront[rowInd];

            PRLEVEL(1, ("%% ri = %ld \n", ri));
            PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
            PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
            dC[ri] += sC[i1];
            PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
        }

        if (ncolsSeen == 0) break;
        PRLEVEL(1, ("\n"));
    }

    // invalidating assembled rows
    for (Int ii = 0; ii < (Int)tempRow.size(); ii++)
    {
        Int i2 = tempRow[ii];
        el_rowIndex[i2] = -1;
        rowRelIndex[i2] = -1;
    }

    el->nrowsleft -= tempRow.size();
    if (el->nrowsleft == 0)
    {
        paru_free_el(e, elementList);
    }
#ifndef NDEBUG
    p = 1;
    PRLEVEL(p, ("%% After Eliminate some rows %ld :\n", eli));
    if (p <= 0) paru_print_element(paruMatInfo, eli);

    PRLEVEL(p, ("%% %ld :\n", e));
    if (p <= 0) paru_print_element(paruMatInfo, e);
#endif
}

void paru_assemble_el_with0rows(Int e, Int f, std::vector<Int> &colHash,
                                paru_matrix *paruMatInfo)

{
    // This element contributes to both pivotal rows and pivotal columns
    //  However it has zero rows in current pivotal columns therefore
    //  not all rows are there
    // it can be eliminated partially
    //       ________________________________
    //       |      |                         |
    //       |      |                         |
    //       ___xxxxxxxxxxx____________________
    //       |  xxxxxxxxxxx                   |
    //       |  oxxo|oxoxox                   | <- assemble these rows
    //       |  ooxx|oxoxox                   |  and mark them assembled
    //       |  oooo|oxoxox                   |
    //       ---------------------------------
    //          ooooooxxxxx  --> outsidie the front
    //          ooooooxxxxx
    //
    //
    DEBUGLEVEL(0);
#ifndef NDEBUG
    Int p = -1;
#endif
    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int *snM = LUsym->super2atree;
    Int eli = snM[f];
    PRLEVEL(p, ("%% \n+++++++++++++++++++++++++++++++++++++++\n"));
    PRLEVEL(p, ("%% Eliminat elment %ld  with0rows in %ld\n", e, eli));

#ifndef NDEBUG
    p = 1;
    PRLEVEL(p, ("%% %ld :\n", eli));
    if (p <= 0) paru_print_element(paruMatInfo, eli);

    PRLEVEL(p, ("%% %ld :\n", e));
    if (p <= 0) paru_print_element(paruMatInfo, e);

#endif

    paru_Element **elementList = paruMatInfo->elementList;

    paru_Element *el = elementList[e];
    paru_Element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    ASSERT(el->nzr_pc > 0);

    // Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

    if (el->cValid != paruMatInfo->time_stamp[f])
        paru_update_rel_ind_col(e, f, colHash, paruMatInfo);

    // Int *colRelIndex = relColInd (paru_Element *el);
    Int *colRelIndex = (Int *)(el + 1) + mEl + nEl;

    // Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int *)(el + 1) + nEl;

    // double *el_Num = numeric_pointer (el);
    double *el_Num = (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);
    // current elemnt numerical pointer
    // double *el_Num = numeric_pointer (curEl);
    double *curEl_Num =
        (double *)((Int *)(curEl + 1) + 2 * curEl->nrows + 2 * curEl->ncols);

    work_struct *Work = paruMatInfo->Work;
    Int *isRowInFront = Work->rowSize;

#ifndef NDEBUG
    paru_fac *Us = paruMatInfo->partial_Us;
    Int *fcolList = paruMatInfo->fcolList[f];
    Int colCount = Us[f].n;
    ASSERT(el_colIndex[el->lac] <= fcolList[colCount - 1]);
    ASSERT(el_colIndex[nEl - 1] <= 0 || fcolList[0] <= el_colIndex[nEl - 1]);
#endif

    PRLEVEL(p, ("%% newColSet.size = %ld\n", colCount));
    PRLEVEL(p, ("%% nEl = %ld\n", nEl));

    if (el->ncolsleft == 1)
    {
        PRLEVEL(p, ("%% 1 col left\n %%"));
        double *sC = el_Num + mEl * el->lac;  // source column pointer
#ifndef NDEBUG
        Int colInd = el_colIndex[el->lac];
        PRLEVEL(1, ("%% colInd =%ld \n", colInd));
        ASSERT(colInd >= 0);
#endif
        Int fcolind = colRelIndex[el->lac];
        double *dC = curEl_Num + fcolind * curEl->nrows;
        Int nrows2bSeen = el->nrowsleft;
        for (Int i = 0; i < mEl; i++)
        {
            Int rowInd = el_rowIndex[i];
            PRLEVEL(1, ("%% rowInd =%ld \n", rowInd));
            if (rowInd >= 0)
            {
                // FIXME
                if (rowRelIndex[i] != -1)  // row with at least one nz
                {
                    Int ri = isRowInFront[rowInd];
                    PRLEVEL(1, ("%% ri = %ld \n", ri));
                    PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
                    PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
                    dC[ri] += sC[i];
                    PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
                }
                if (--nrows2bSeen == 0) break;
            }
        }
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
            PRLEVEL(1, ("%% rowInd =%ld ", rowInd));
#ifndef NDEBUG
            if (rowRelIndex[i] == -1) PRLEVEL(-1, ("%% row_with0 "));
#endif

            PRLEVEL(1, ("%% \n"));
            if (rowInd >= 0 && rowRelIndex[i] != -1)
            {
                tempRow[ii++] = i;
                Int ri = isRowInFront[rowInd];
                rowRelIndex[i] = ri;
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
        Int ncols2bSeen = el->ncolsleft;
        for (Int j = el->lac; j < nEl; j++)
        {
            PRLEVEL(1, ("%% j =%ld \n", j));
            double *sC = el_Num + mEl * j;  // source column pointer
            Int colInd = el_colIndex[j];
            PRLEVEL(1, ("%% colInd =%ld \n", colInd));
            if (colInd < 0) continue;
            Int fcolind = colRelIndex[j];

            double *dC = curEl_Num + fcolind * curEl->nrows;

            for (Int iii = 0; iii < nrows2assembl; iii++)
            {
                Int i = tempRow[iii];
                Int ri = rowRelIndex[i];
                ASSERT(rowRelIndex[i] != -1);  // I already picked the rows
                // that are not in zero pivots
                ASSERT(rowInd >= 0);  // and also still alive
                PRLEVEL(1, ("%% ri = %ld rowInd=%ld\n", ri, rowInd));
                PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
                dC[ri] += sC[i];
                PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
            }

            if (--ncols2bSeen == 0) break;
            PRLEVEL(1, ("\n"));
        }
    }

    // Mark rows as assembled and updat lac

    Int nrows2bSeen = el->nrowsleft;
    Int new_lac = nEl;
    for (Int ii = 0; ii < mEl; ii++)
    {
        Int rowInd = el_rowIndex[ii];
        if (rowInd < 0) continue;  // already gone

        if (rowRelIndex[ii] == -1)  // row with all zeros in piv
        {                           // update lac
            PRLEVEL(1, ("%%Searching for lac in %ld\n%%", rowInd));
            PRLEVEL(1, ("%%col=%ld\n%%", el->lac));
            for (Int jj = el->lac; jj < new_lac; jj++)
            // searching for the first nz
            {
                if (el_colIndex[jj] < 0) continue;
                // el [rowInd, jj]
                PRLEVEL(-1, ("%% el[%ld,%ld]=%2.5lf\n%%", rowInd, jj,
                             el_Num[mEl * jj + ii]));
                if (el_Num[mEl * jj + ii] != 0)
                {
                    new_lac = jj;
                    PRLEVEL(-1, ("%%Found new-lac in %ld\n%%", jj));
                    break;
                }
            }
        }
        else  // It was assembled here; mark row as assembled
        {
            el_rowIndex[ii] = -1;
        }
        if (--nrows2bSeen == 0) break;
    }
    // updating lac can have effect on number of columns left
    // I should update number of columns left too

    if (new_lac != el->lac)
    {
        Int ncolsleft = 0;
        for (Int j = new_lac; j < nEl; j++)
        {
            if (el_colIndex[j] > 0) ncolsleft++;
        }
        PRLEVEL(-1, ("%%colsleft was %ld and now is %ld\n%%", el->ncolsleft,
                     ncolsleft));
        el->ncolsleft = ncolsleft;
        for (Int j = el->lac; j < new_lac; j++)
        {
            if (el_colIndex[j] >= 0) el_colIndex[j] = flip(el_colIndex[j]);
        }
    }

    el->nrowsleft = el->nzr_pc;
    el->lac = new_lac;
    Int *lacList = paruMatInfo->lacList;
    lacList[e] = el_colIndex[el->lac];
#ifndef NDEBUG
    Int *Super = LUsym->Super;
    Int col1 = Super[f]; /* fornt F has columns col1:col2-1 */
    Int col2 = Super[f + 1];
    p = -1;
    PRLEVEL(p, ("%% %ld(%ld) %ld-%ld :\n", f, eli, col1, col2));
    PRLEVEL(p, ("%%Finally new-lac is %ld ", el->lac));
    PRLEVEL(p, ("nEl=%ld\n lacList[%ld]=%ld nrowsleft=%ld\n", nEl, e,
                lacList[e], el->nrowsleft));

    p = 1;
    if (nEl != new_lac && el_colIndex[new_lac] < col2) p = -2;

    PRLEVEL(p, ("%% %ld :\n", eli));
    if (p <= 0) paru_print_element(paruMatInfo, eli);

    PRLEVEL(p, ("%% %ld :\n", e));
    if (p <= 0) paru_print_element(paruMatInfo, e);
    p = 1;
    ASSERT(nEl == new_lac || col2 <= el_colIndex[new_lac]);

#endif
    if (new_lac == nEl)
    {
#ifndef NDEBUG
        PRLEVEL(p, ("%% %ld is freed inside with0\n", eli));
#endif
        paru_free_el(e, elementList);
    }
}
