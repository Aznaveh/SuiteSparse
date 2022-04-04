////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_print   ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief printing datas are implemented in this file
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
void paru_print_element(Int e, paru_work *Work, ParU_Numeric *Num)
{
    // print out contribution blocks
    paru_element **elementList;
    elementList = Work->elementList;
    paru_element *curEl = elementList[e];

    Int morign = Num->m;
    Int nf = Num->Sym->nf;

    if (e > morign + nf + 1)
    {
        printf("%% paru_element %ld is out of range; just %ld elements \n", e,
               morign + nf + 1);
        return;
    }

    if (curEl == NULL)
    {
        printf("%% paru_element %ld is empty\n", e);
        return;
    }

    Int m, n;
    m = curEl->nrows;
    n = curEl->ncols;

    Int *el_colIndex = colIndex_pointer(curEl);
    Int *el_rowIndex = rowIndex_pointer(curEl);

    // Int *rel_col = relColInd (curEl);
    // Int *rel_row = relRowInd (curEl);

    double *el_colrowNum = numeric_pointer(curEl);

    printf("\n");
    printf("%% paru_element %ld is %ld x %ld:\n", e, m, n);

    printf("\t");
    //    for (int j = 0; j < n; j++)
    //        printf("%% %ld\t", rel_col[j] );
    //    printf("\n\t");
    for (int j = 0; j < n; j++) printf("%% %ld\t", el_colIndex[j]);

    printf("\n");
    for (int i = 0; i < m; i++)
    {
        //     printf("%% %ld\t %ld\t",rel_row[i], el_rowIndex [i] );
        printf("%% %ld\t", el_rowIndex[i]);
        for (int j = 0; j < n; j++)
        {
            double value = el_colrowNum[j * m + i];
            printf("%2.4lf\t", value);
        }
        printf("\n");
    }
}

void paru_print_paru_tupleList(paru_tupleList *listSet, Int index)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% listSet =%p\n", listSet));

    if (listSet == NULL)
    {
        printf("%% Empty tuple\n");
        return;
    }

    paru_tupleList cur = listSet[index];
    Int numTuple = cur.numTuple;
    paru_tuple *l = cur.list;

    printf("%% There are %ld tuples in this list:\n %%", numTuple);
    for (Int i = 0; i < numTuple; i++)
    {
        paru_tuple curTpl = l[i];
        printf(" (%ld,%ld)", curTpl.e, curTpl.f);
    }
    printf("\n");
}
