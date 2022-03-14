////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_init_rel /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief Initiazing element f's  time_stamp
 *          chekc all f's children and find the maximum time_stamp
 *          this time is checking for validation or invalidation of elements
 *
 *
 * @param paruMatInfo: pointer to matrix info
 *        f: front that is going to be initialized;
 *
 * @author Aznaveh
 */
#include "paru_internal.hpp"
void paru_init_rel(Int f, paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(0);
    paru_symbolic *Sym = paruMatInfo->Sym;
    Int *time_stamp = paruMatInfo->time_stamp;

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
