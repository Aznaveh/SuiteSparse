/** =========================================================================  /
 * =======================  paru_init_rel ===================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"
/*! @brief Initiazing element f's rowRelIndValid and colRelIndValid
 *          chekc all f's children and find the maximum relIndex and ensure that
 *          f's relative index is bigger than all its children
 *          if any chidlren is NULL it means that it has been vaporized along 
 *          the way and it's children must be checked. It is done recursively
 *
 * @param LUsym: pointer to symboicl analysis
 *        elementList: list of element pointers that are intiailized by NULL 
 *        f: front that is going to be initialized; its element number is
 *        computed in augmented tree
 *        maxrValid, maxcValid: local varibales in paru_assemble that will be
 *        updated
 *
 */

void paru_init_rel  (
    paru_symbolic *LUsym,
    Element **elementList,
    Int f,
    Int *maxrValid, 
    Int *maxcValid)
{
    DEBUGLEVEL(0);
   Int *Child = LUsym->Child;
    Int *Childp = LUsym->Childp;
    Int *snM = LUsym->super2atree;
    PRLEVEL (1, ("%% children of %ld  are:\n", f));
    Int maxrV, maxcV;
    maxrV= maxcV= -1;
    for(Int p = Childp[f]; p <= Childp[f+1]-1; p++){
        Int rV, cV;
        ASSERT(Child[p] >= 0);
        Element *childEl = elementList[snM[Child[p]]];
        if(childEl == NULL){
            PRLEVEL (1, ("%% Vaporized element=%ld  ",snM[Child[p]] ));
            PRLEVEL (1, ("%% Child[%ld]= %ld  ",p,Child[p] ));
            paru_init_rel (LUsym, elementList, Child[p], &rV, &cV);
        }
        else{
            rV = childEl->rValid;
            cV = childEl->cValid;
        }
        PRLEVEL (0, ("%% Child[%ld]= %ld  ",p,Child[p] ));
        PRLEVEL (0, ("%% rV=%ld cV= %ld\n  ",rV,cV ));
        maxrV = maxrV > rV? maxrV : rV;
        maxcV = maxcV > cV? maxcV : cV;
    }
    *maxrValid = ++maxrV; 
    *maxcValid = ++maxcV; 
}


