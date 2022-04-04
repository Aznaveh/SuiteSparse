////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_assemble //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  finding the  columns or rows of prior element and fully or partially
 * assemble it  and eliminate it if needed
 *
 *  @author Aznaveh
 */
#include "paru_internal.hpp"

void paru_assemble_all(Int e, Int f, std::vector<Int> &colHash,
        paru_work *Work, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;
#ifndef NTIME
    static double tot_assem_time = 0;
    double start_time = PARU_OPENMP_GET_WTIME;
#endif

    ParU_Symbolic *Sym = Work->Sym;
    Int *snM = Sym->super2atree;
    Int eli = snM[f];
    PRLEVEL(PR, ("%% Eliminate all of %ld in %ld(f=%ld) (tid=%d)\n", e, eli, f,
                 PARU_OPENMP_GET_THREAD_ID));

#ifndef NDEBUG
    PR = 0;
    PRLEVEL(PR, ("%% %ld :\n", e));
    if (PR <= 0) paru_print_element(e, Work, Num);

    PRLEVEL(PR, ("%% %ld :\n", eli));
    if (PR <= 0) paru_print_element(eli Work, Num);
    PR = 1;
#endif

    paru_element **elementList = Work->elementList;

    paru_element *el = elementList[e];
    paru_element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    // Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

    if (el->cValid != Work->time_stamp[f])
        paru_update_rel_ind_col(e, f, colHash, Work, Num);

    // Int *colRelIndex = relColInd (paru_element *el);
    Int *colRelIndex = (Int *)(el + 1) + mEl + nEl;

    // Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int *)(el + 1) + nEl;

    // double *el_Num = numeric_pointer (el);
    double *el_Num = (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);
    // current elemnt numerical pointer
    // double *el_Num = numeric_pointer (curEl);
    double *curEl_Num =
        (double *)((Int *)(curEl + 1) + 2 * curEl->nrows + 2 * curEl->ncols);

    Int *isRowInFront = Work->rowSize;

#ifndef NDEBUG
    ParU_Factors *Us = Num->partial_Us;
    Int *fcolList = Num->fcolList[f];
    Int colCount = Us[f].n;
    ASSERT(el_colIndex[el->lac] <= fcolList[colCount - 1]);
    ASSERT(el_colIndex[nEl - 1] <= 0 || fcolList[0] <= el_colIndex[nEl - 1]);
    PRLEVEL(PR, ("%% newColSet.size = %ld\n", colCount));
    PRLEVEL(PR, ("%% nEl = %ld\n", nEl));
#endif

    if (el->ncolsleft == 1)
    {
        PRLEVEL(PR, ("%% 1 col left\n %%"));
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
        PRLEVEL(PR, ("%% more than 1 col left (%ld->%ld)\n %%", e, eli));

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

        Int naft;  // number of active frontal tasks
        #pragma omp atomic read
        naft = Work->naft;
        ParU_Control *Control = Num->Control;
        const Int max_threads = Control->paru_max_threads;

        if (naft > max_threads / 2 || el->nrowsleft * el->ncolsleft < 4096 ||
            el->nrowsleft < 1024)
        {  // not enoght resources or very small assembly
            // sequential
            PRLEVEL(1,
                    ("Seqntial Assembly naft=%ld colsleft=%ld rowsleft=%ld \n",
                     naft, el->ncolsleft, el->nrowsleft));

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
        else
        {
            // enoght threads and big assembly
            // go parallel
            PRLEVEL(1, ("Parallel Assembly naft=%ld colsleft=%ld rowsleft=%ld "
                        "el->lac = %ld nEl=%ld rem =%ld (%ld->%ld)\n",
                        naft, el->ncolsleft, el->nrowsleft, el->lac, nEl,
                        nEl - el->lac, e, eli));

            // // each column a tsk
            //#..pragma omp parallel proc_bind(close)
            // num_threads(max_threads / naft)
            //#..pragma omp single nowait
            //#..pragma omp task untied
            // for (Int j = el->lac; j < nEl; j++)
            //{
            //    PRLEVEL(1, ("%% j =%ld \n", j));
            //    double *sC = el_Num + mEl * j;  // source column pointer
            //    Int colInd = el_colIndex[j];
            //    PRLEVEL(1, ("%% colInd =%ld \n", colInd));
            //    if (colInd < 0) continue;

            //    Int fcolind = colRelIndex[j];

            //    double *dC = curEl_Num + fcolind * curEl->nrows;

            //    #pragma omp task
            //    for (Int iii = 0; iii < el->nrowsleft; iii++)
            //    {
            //        Int i = tempRow[iii];
            //        Int ri = rowRelIndex[i];

            //        PRLEVEL(1, ("%% ri = %ld \n", ri));
            //        PRLEVEL(1, ("%% sC [%ld] =%2.5lf \n", i, sC[i]));
            //        PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", ri, dC[ri]));
            //        dC[ri] += sC[i];
            //        PRLEVEL(1, ("%% dC [%ld] =%2.5lf \n", i, dC[ri]));
            //    }
            //    if (--el->ncolsleft == 0) break;
            //    PRLEVEL(1, ("\n"));
            //}

            ///////////////////////////////////////////////////////////////////
            /////////////// making tasks and such /////////////////////////////
            ///////////////////////////////////////////////////////////////////

            Int ntasks = (max_threads - naft + 1) * 2;
            ntasks = ntasks == 0 ? 1 : ntasks;
            Int task_size = (nEl - el->lac) / ntasks;
            PRLEVEL(1, ("BBB el->lac=%ld nEl=%ld ntasks=%ld task_size=%ld\n",
                        el->lac, nEl, ntasks, task_size));
            if (task_size == 0 || task_size == 1)
            {
                task_size = 1;
                ntasks = nEl - el->lac;
            }
            PRLEVEL(1, ("el->lac=%ld nEl=%ld ntasks=%ld task_size=%ld\n",
                        el->lac, nEl, ntasks, task_size));
            #pragma omp parallel proc_bind(close)
            #pragma omp single
            #pragma omp task
            for (Int t = 0; t < ntasks; t++)
            {
                Int c1 = el->lac + t * task_size;
                Int c2 = el->lac + (t + 1) * task_size;
                c2 = t == ntasks - 1 ? nEl : c2;
                PRLEVEL(1, ("t=%ld c1=%ld c2=%ld\n", t, c1, c2));
                #pragma omp task mergeable
                for (Int j = c1; j < c2; j++)
                {
                    PRLEVEL(1, ("%% j =%ld t=%ld\n", j, t));
                    double *sC = el_Num + mEl * j;  // source column pointer
                    Int colInd = el_colIndex[j];
                    PRLEVEL(1, ("%% colInd =%ld \n", colInd));
                    if (colInd < 0) continue;
                    PRLEVEL(1, ("inside paralle region %d j=%ld (tid=%d)\n",
                                PARU_OPENMP_GET_ACTIVE_LEVEL, j,
                                PARU_OPENMP_GET_THREAD_NUM));
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
                    PRLEVEL(1, ("\n"));
                }
            }
        }
    }
    paru_free_el(e, elementList);
#ifndef NDEBUG
    PR = 0;
    PRLEVEL(PR, ("%% after assembly %ld :\n", eli));
    if (PR <= 0) paru_print_element(eli, Work, Num);
    PR = 1;
#endif

#ifndef NTIME
    double time = PARU_OPENMP_GET_WTIME;
    time -= start_time;
    #pragma omp atomic update
    tot_assem_time += time;
    if (f > Sym->nf - 5)
        PRLEVEL(-1, ("%% assemble all %ld\t->%ld\t took %lf seconds tot=%lf\n",
                     e, eli, time, tot_assem_time));
#endif
}

// try to find columns and assemble them to current front. After the first
// column that is not in current front it gets a toll for each column doesn't
// fit

void paru_assemble_cols(Int e, Int f, std::vector<Int> &colHash,
        paru_work *Work, ParU_Numeric *Num)

{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;
#ifndef NDEBUG
    Int c = 0;  // number of columns assembled
#endif
    ParU_Symbolic *Sym = Work->Sym;
    Int *snM = Sym->super2atree;
    Int eli = snM[f];

    PRLEVEL(PR, ("%% Eliminat some cols of %ld in %ld\n", e, eli));
#ifndef NDEBUG
    PR = 1;

    PRLEVEL(PR, ("%% %ld :\n", eli));
    if (PR <= 0) paru_print_element(eli, Work, Num);

    PRLEVEL(PR, ("%% %ld :\n", e));
    if (PR <= 0) paru_print_element(e, Work, Num);
#endif

    paru_element **elementList = Work->elementList;

    paru_element *el = elementList[e];
    paru_element *curEl = elementList[eli];

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

    Int *isRowInFront = Work->rowSize;

    Int *fcolList = Num->fcolList[f];

    // Int tempRow[el->nrowsleft];  // C99
    std::vector<Int> tempRow(el->nrowsleft);
    Int tempRow_ready = 0;
    Int toll = 8;  // number of times it continue when do not find anything

    // Int naft; //number of active frontal tasks
    // pragma omp atomic read
    // naft = Num->naft;
    // const Int max_threads = Num->paru_max_threads;
    ////Int *Depth = Sym->Depth;
    // pragma omp parallel proc_bind(close) num_threads(max_threads/naft)
    // if (naft < max_threads/2 &&
    //        el->nrowsleft*el->ncolsleft < 4096 && el->nrowsleft < 1024 )
    // pragma omp single nowait
    // pragma omp task untied

    // TOLL FREE zone
    while (paru_find_hash(el_colIndex[el->lac], colHash, fcolList) != -1)
    {
        PRLEVEL(PR, ("%% Toll free\n"));
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

        // pragma omp task
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
    Int *lacList = Num->lacList;
    lacList[e] = el_colIndex[el->lac];

    // TOLL Zone
    //**//pragma omp parallel
    //**//pragma omp single nowait
    //**//pragma omp taskgroup
    for (Int j = el->lac + 1; j < nEl && el->ncolsleft > 0 && toll > 0; j++)
    {
        PRLEVEL(1, ("%% Toll zone\n"));
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

        //**//pragma omp task priority(Depth[f]) if(el->nrowsleft > 1024)
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

#ifndef NDEBUG
    PRLEVEL(1, ("%%  %ld has found and assembled, ncolsleft %ld\n", c,
                el->ncolsleft));
#endif

    if (el->ncolsleft == 0)
    {
        paru_free_el(e, elementList);
    }
}

void paru_assemble_rows(Int e, Int f, std::vector<Int> &colHash,
        paru_work *Work, ParU_Numeric *Num)

{
    DEBUGLEVEL(0);
    PARU_DEFINE_PRLEVEL;

    ParU_Symbolic *Sym = Work->Sym;
    Int *snM = Sym->super2atree;
    Int eli = snM[f];

    PRLEVEL(PR, ("%% Eliminat some rows of %ld in %ld\n", e, eli));

    paru_element **elementList = Work->elementList;

    paru_element *el = elementList[e];
    paru_element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    // Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

    // Int *colRelIndex = relColInd (paru_element *el);
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
        PRLEVEL(PR, ("%% Toll free zone: %ld rows has been found: \n%%",
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

    PRLEVEL(PR,
            ("%% %ld rows has been found, toll %ld\n%%", tempRow.size(), toll));
#ifndef NDEBUG
    for (Int ii = 0; ii < (Int)tempRow.size(); ii++)
        PRLEVEL(PR, ("%ld ", tempRow[ii]));
    PRLEVEL(PR, ("\n "));
#endif
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%% Before eliminiatine some rows %ld :\n", eli));
    if (PR <= 0) paru_print_element(eli, Work, Num);

    PRLEVEL(PR, ("%% %ld :\n", e));
    if (PR <= 0) paru_print_element(e, Work, Num);
#endif

    if (el->cValid != Work->time_stamp[f])
        paru_update_rel_ind_col(e, f, colHash, Work, Num);

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
    PR = 1;
    PRLEVEL(PR, ("%% After Eliminate some rows %ld :\n", eli));
    if (PR <= 0) paru_print_element(eli, Work, Num);

    PRLEVEL(PR, ("%% %ld :\n", e));
    if (PR <= 0) paru_print_element(e, Work, Num);
#endif
}

void paru_assemble_el_with0rows(Int e, Int f, std::vector<Int> &colHash,
        paru_work *Work, ParU_Numeric *Num)

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
    PARU_DEFINE_PRLEVEL;

    ParU_Symbolic *Sym = Work->Sym;
    Int *snM = Sym->super2atree;
    Int eli = snM[f];
    PRLEVEL(PR, ("%% \n+++++++++++++++++++++++++++++++++++++++\n"));
    PRLEVEL(PR, ("%% Eliminat elment %ld  with0rows in %ld\n", e, eli));

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%% %ld :\n", eli));
    if (PR <= 0) paru_print_element(eli, Work, Num);

    PRLEVEL(PR, ("%% %ld :\n", e));
    if (PR <= 0) paru_print_element(e, Work, Num);

#endif

    paru_element **elementList = Work->elementList;

    paru_element *el = elementList[e];
    paru_element *curEl = elementList[eli];

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    ASSERT(el->nzr_pc > 0);

    // Int *el_colIndex = colIndex_pointer (el);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *rowRelIndex = relRowInd (el);
    Int *rowRelIndex = (Int *)(el + 1) + 2 * nEl + mEl;

    if (el->cValid != Work->time_stamp[f])
        paru_update_rel_ind_col(e, f, colHash, Work, Num);

    // Int *colRelIndex = relColInd (paru_element *el);
    Int *colRelIndex = (Int *)(el + 1) + mEl + nEl;

    // Int *el_rowIndex = rowIndex_pointer (el);
    Int *el_rowIndex = (Int *)(el + 1) + nEl;

    // double *el_Num = numeric_pointer (el);
    double *el_Num = (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);
    // current elemnt numerical pointer
    // double *el_Num = numeric_pointer (curEl);
    double *curEl_Num =
        (double *)((Int *)(curEl + 1) + 2 * curEl->nrows + 2 * curEl->ncols);

    Int *isRowInFront = Work->rowSize;

#ifndef NDEBUG
    ParU_Factors *Us = Num->partial_Us;
    Int *fcolList = Num->fcolList[f];
    Int colCount = Us[f].n;
    ASSERT(el_colIndex[el->lac] <= fcolList[colCount - 1]);
    ASSERT(el_colIndex[nEl - 1] <= 0 || fcolList[0] <= el_colIndex[nEl - 1]);
    PRLEVEL(PR, ("%% newColSet.size = %ld\n", colCount));
    PRLEVEL(PR, ("%% nEl = %ld\n", nEl));
#endif

    if (el->ncolsleft == 1)
    {
        PRLEVEL(PR, ("%% 1 col left\n %%"));
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
        PRLEVEL(PR, ("%% more than 1 col left\n %%"));

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
            if (rowRelIndex[i] == -1) PRLEVEL(1, ("%% row_with0 "));
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
        PR = 1;
        PRLEVEL(PR, ("%% list of the rows to be assembled:\n%%"));
        for (Int i = 0; i < nrows2assembl; i++)
            PRLEVEL(PR, ("%ld ", el_rowIndex[tempRow[i]]));
        PRLEVEL(PR, ("%% \n"));
#endif
        Int ncols2bSeen = el->ncolsleft;
        // Int *Depth = Sym->Depth;
        //**//pragma omp parallel
        //**//pragma omp single nowait
        //**//pragma omp taskgroup
        for (Int j = el->lac; j < nEl; j++)
        {
            PRLEVEL(1, ("%% j =%ld \n", j));
            double *sC = el_Num + mEl * j;  // source column pointer
            Int colInd = el_colIndex[j];
            PRLEVEL(1, ("%% colInd =%ld \n", colInd));
            if (colInd < 0) continue;
            Int fcolind = colRelIndex[j];

            double *dC = curEl_Num + fcolind * curEl->nrows;

            //**//pragma omp task priority(Depth[f]) if(nrows2assembl > 1024)
            for (Int iii = 0; iii < nrows2assembl; iii++)
            {
                Int i = tempRow[iii];
                Int ri = rowRelIndex[i];
                ASSERT(rowRelIndex[i] != -1);  // I already picked the rows
                // that are not in zero pivots
                PRLEVEL(1, ("%% ri = %ld \n", ri));
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
                PRLEVEL(1, ("%% el[%ld,%ld]=%2.5lf\n%%", rowInd, jj,
                            el_Num[mEl * jj + ii]));
                if (el_Num[mEl * jj + ii] != 0)
                {
                    new_lac = jj;
                    PRLEVEL(1, ("%%Found new-lac in %ld\n%%", jj));
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
        PRLEVEL(1, ("%%colsleft was %ld and now is %ld\n%%", el->ncolsleft,
                    ncolsleft));
        el->ncolsleft = ncolsleft;
        for (Int j = el->lac; j < new_lac; j++)
        {
            if (el_colIndex[j] >= 0) el_colIndex[j] = flip(el_colIndex[j]);
        }
    }

    el->nrowsleft = el->nzr_pc;
    el->lac = new_lac;
    Int *lacList = Num->lacList;
    lacList[e] = el_colIndex[el->lac];
#ifndef NDEBUG
    Int *Super = Sym->Super;
    Int col1 = Super[f]; /* fornt F has columns col1:col2-1 */
    Int col2 = Super[f + 1];
    PR = 1;
    PRLEVEL(PR, ("%% %ld(%ld) %ld-%ld :\n", f, eli, col1, col2));
    PRLEVEL(PR, ("%%Finally new-lac is %ld ", el->lac));
    PRLEVEL(PR, ("nEl=%ld\n lacList[%ld]=%ld nrowsleft=%ld\n", nEl, e,
                 lacList[e], el->nrowsleft));

    PR = 1;
    if (nEl != new_lac && el_colIndex[new_lac] < col2) PR = -2;

    PRLEVEL(PR, ("%% %ld :\n", eli));
    if (PR <= 0) paru_print_element(eli, Work, Num);

    PRLEVEL(PR, ("%% %ld :\n", e));
    if (PR <= 0) paru_print_element(e, Work, Num);
    PR = 1;
    ASSERT(nEl == new_lac || col2 <= el_colIndex[new_lac]);

#endif
    if (new_lac == nEl)
    {
#ifndef NDEBUG
        PRLEVEL(PR, ("%% %ld is freed inside with0\n", eli));
#endif
        paru_free_el(e, elementList);
    }
}
