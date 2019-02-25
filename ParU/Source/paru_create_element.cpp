/** =========================================================================  /
 * =======================  paru_create_element  ============================  /
 * ==========================================================================  /
 * @brief Initializing an empty element
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
Element *paru_create_element (Int nrows, Int ncols, 
        Int init, cholmod_common *cc)
{
    DEBUGLEVEL(0);

    PRLEVEL (1, ("%% creating %ldx%ld element ", nrows, ncols));
    Element *curEl;
    Int tot_size = sizeof(Element)+ sizeof(Int)*(2*(nrows+ncols))+
                sizeof(double)*nrows*ncols;
    if (init)
        curEl = (Element*) paru_calloc(1, tot_size , cc);
    else
        curEl = (Element*) paru_alloc(1, tot_size , cc);
    if (curEl == NULL) return NULL; // do not do error checking

    PRLEVEL (1, ("%% with size of %ld in %p\n", tot_size, curEl));
    
    // Initializing current element
    curEl->nrowsleft = curEl->nrows = nrows;
    curEl->ncolsleft = curEl->ncols = ncols;
    curEl->rValid = -1;
    curEl->cValid = -1;

    return curEl;
}


