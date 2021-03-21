/** =========================================================================  /
 * =======================  paru_create_element  ============================  /
 * ========================================================================== */
/*! @brief Initializing an empty element
 *    if (Init) use calloc .i.e x = 0;
 *                     
 *                    V RRRRRRRRRRRRRR
 *                      IIIIIIIIIIIIII            I = global index
 *                  V                             r = relative index 
 *                  R I  xxxxxxxxxxxxxx           V = row/col valid bit
 *                  R I  xxxxxxxxxxxxxx            if V == f it is valid for
 *                  R I  xxxxxxxxxxxxxx             current front
 *                  R I  xxxxxxxxxxxxxx
 *                  R I  xxxxxxxxxxxxxx
 *
 * @author Aznaveh
 *  */
#include "Parallel_LU.hpp"
paru_Element *paru_create_element (Int nrows, Int ncols, 
        Int init, paru_matrix *paruMatInfo, cholmod_common *cc)
{
    DEBUGLEVEL(0);

    PRLEVEL (1, ("%% creating %ldx%ld element ", nrows, ncols));
    paru_Element *curEl;
    Int tot_size = sizeof(paru_Element)+ sizeof(Int)*(2*(nrows+ncols))+
                sizeof(double)*nrows*ncols;

    std::pmr::synchronized_pool_resource &pool = *paruMatInfo->pool_p;

    curEl = (paru_Element*) pool.allocate (tot_size, 8);
    /*
    if (init)
        curEl = (paru_Element*) paru_calloc(1, tot_size , cc);
    else
        curEl = (paru_Element*) paru_alloc(1, tot_size , cc);
    */

    if (curEl == NULL) return NULL; // do not do error checking

    PRLEVEL (1, (" with size of %ld in %p\n", tot_size, curEl));
    
    // Initializing current element
    curEl->nrowsleft = curEl->nrows = nrows;
    curEl->ncolsleft = curEl->ncols = ncols;
    curEl->rValid = -1;
    curEl->cValid = -1;

    curEl->lac = 0;

    return curEl;
}


