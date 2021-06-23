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
#include "paru_internal.hpp"
paru_Element *paru_create_element(Int nrows, Int ncols, Int init)
{
    DEBUGLEVEL(0);

    PRLEVEL(1, ("%% creating %ldx%ld element ", nrows, ncols));
    paru_Element *curEl;
    size_t tot_size = sizeof(paru_Element) +
                      sizeof(Int) * (2 * (nrows + ncols)) +
                      sizeof(double) * nrows * ncols;
    if (init)
    {
        curEl = (paru_Element *)paru_calloc(1, tot_size);
    }
    else
    {
        curEl = (paru_Element *)paru_alloc(1, tot_size);
    }
    if (curEl == NULL) return NULL;  // do not do error checking

    PRLEVEL(1, (" with size of %ld in %p\n", tot_size, curEl));

    // Initializing current element
    curEl->nrowsleft = curEl->nrows = nrows;
    curEl->ncolsleft = curEl->ncols = ncols;
    curEl->rValid = -1;
    curEl->cValid = -1;

    curEl->lac = 0;
    return curEl;
}
