/** =========================================================================  /
 * =======================  paru_numeric_assemble ===========================  /
 * ========================================================================== */
/*! @brief  assemble numbers in U part of the matrix.
 *          It is per row, and the matrices are stored in column,
 *          therefore it can reduce the performance in some cases
 *
 * @author Aznaveh
 *
 */

#include "paru_internal.hpp"

void assemble_row_toU(Int e, Int f, Int sR, Int dR, std::vector<Int> &colHash,
                      paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);

    paru_Element **elementList = paruMatInfo->elementList;
    paru_Element *el = elementList[e];

    if (el->cValid != paruMatInfo->time_stamp[f])
        // if not updatated
        paru_update_rel_ind_col(e, f, colHash, paruMatInfo);

    paru_fac *Us = paruMatInfo->partial_Us;
    double *uPart = Us[f].p;  // uPart

    Int nEl = el->ncols;
    Int mEl = el->nrows;

    paru_fac *LUs = paruMatInfo->partial_LUs;
    Int fp = LUs[f].n;

    // Int *el_colIndex = colIndex_pointer (curEl);
    Int *el_colIndex = (Int *)(el + 1);

    // Int *colRelIndex = relColInd (paru_Element *el);
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
            uPart[rj * fp + dR] += sM[mEl * j + sR];
            PRLEVEL(1, ("%% uPart [%ld] =%2.5lf \n", rj * fp + dR,
                        uPart[rj * fp + dR]));
            if (ncolsSeen == 0) break;
        }
    }
}
