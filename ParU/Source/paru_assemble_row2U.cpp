////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_assemble_row2U.cpp ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief  assemble numbers in U part of the matrix.
 *          It is per row, and the matrices are stored in column,
 *          therefore it can reduce the performance in some cases.
 *
 * @author Aznaveh
 *
 */

#include "paru_internal.hpp"

void paru_assemble_row_2U(Int e, Int f, Int sR, Int dR,
                          std::vector<Int> &colHash, 
                          paru_work *Work, ParU_Numeric *Num)
{
    DEBUGLEVEL(0);

    paru_element **elementList = Work->elementList;
    paru_element *el = elementList[e];

    if (el->cValid != Work->time_stamp[f])
        // if not updatated
        paru_update_rel_ind_col(e, f, colHash, Work, Num);

    ParU_Factors *Us = Num->partial_Us;
    double *uPart = Us[f].p;  // uPart

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    ParU_Factors *LUs = Num->partial_LUs;
    Int fp = LUs[f].n;

    // Int *el_colIndex = colIndex_pointer (curEl);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *colRelIndex = relColInd (paru_element *el);
    Int *colRelIndex = (Int *)(el + 1) + mEl + nEl;

    // double *el_Num = numeric_pointer (el);
    double *sM = (double *)((Int *)(el + 1) + 2 * nEl + 2 * mEl);

    Int ncolsSeen = el->ncolsleft;
    for (Int j = el->lac; j < nEl; j++)
    {
        Int rj = colRelIndex[j];
        if (el_colIndex[j] >= 0)
        {  // If still valid
            ncolsSeen--;
            PRLEVEL(1,
                    ("%% sM [%ld] =%2.5lf \n", mEl * j + sR, sM[mEl * j + sR]));
            PRLEVEL(1, ("%% uPart [%ld] =%2.5lf \n", rj * fp + dR,
                        uPart[rj * fp + dR]));
            //**//#pragma omp atomic
            uPart[rj * fp + dR] += sM[mEl * j + sR];
            PRLEVEL(1, ("%% uPart [%ld] =%2.5lf \n", rj * fp + dR,
                        uPart[rj * fp + dR]));
            if (ncolsSeen == 0) break;
        }
    }
}
