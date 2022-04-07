////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_init_rel /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0

/*! @brief Initiazing element f's  time_stamp
 *          chekc all f's children and find the maximum time_stamp
 *          this time is checking for validation or invalidation of elements
 *
 *
 * @param Num: pointer to matrix info
 *        f: front that is going to be initialized;
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
void paru_init_rel(Int f, paru_work *Work)
{
    DEBUGLEVEL(0);
    ParU_Symbolic *Sym = Work->Sym;
    Int *time_stamp = Work->time_stamp;

    Int *Child = Sym->Child;
    Int *Childp = Sym->Childp;
    Int max_time = 0;

    PRLEVEL(1, ("%% begining=%ld end=%ld \n", Childp[f], Childp[f + 1]));
    PRLEVEL(1, ("%% children of %ld  are:\n", f));
    for (Int p = Childp[f]; p <= Childp[f + 1] - 1; p++)
    {
        Int child_rel;
        ASSERT(Child[p] >= 0);
        child_rel = time_stamp[Child[p]];
        PRLEVEL(1, ("%% Child[%ld]= %ld  ", p, Child[p]));
        max_time = max_time > child_rel ? max_time : child_rel;
    }
    time_stamp[f] = ++max_time;
    PRLEVEL(1, ("%% max_time=%ld \n", max_time));
}
