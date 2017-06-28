#include "Parallel_LU.hpp"

/** =========================================================================  /
 * =======================  paru_add_tuples =================================  /
 * ==========================================================================  /
/*! \brief assembling a front and updating correspoing elelment
 *
 *
 * \param  a list of tuples and the the tuple we want to add
 * \return 0 on sucess 
 */
void paru_assemble(
    paru_matrix *paruMatInfo,
    //Row/Col list/tuples and LUsym handle
    Int f)

{
    DEBUGLEVEL(1);
    paru_symbolic *LUsym =  paruMatInfo->LUsym;

    Int m,n,nf;
    m = paruMatInfo-> m;
    n = paruMatInfo-> n;
    nf = paruMatInfo->LUsym->nf;
    Int *Super = LUsym->Super;
    /* ---------------------------------------------------------------------- */
    /* get the front F  */
    /* ---------------------------------------------------------------------- */

    PRLEVEL (1, ("Working on Front %ld\n", f));
    /* pivotal columns Super [f] ... Super [f+1]-1 */
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];
    Int *Rp = LUsym->Rp;
    Int *Rj = LUsym->Rj;
    Int p1 = Rp [f];        /* Rj [p1:p2-1] = columns in F */
    Int p2 = Rp [f+1];
    Int fp = col2 - col1;   /* first fp columns are pivotal */ 
    Int fn = p2 - p1;           /* upper bound in number of columns of F */

    PRLEVEL (1, ("%ld pivotal columns:%ld...%ld\n", fp, col1, col2));
    PRLEVEL (1, ("upper bound on columns: %ld ... %ld\n", Rj [p1], Rj [p2-1]));
    PRLEVEL (1, ("p1=%ld p2=%ld\n",p1 ,p2));
    ASSERT (fp > 0 );
    ASSERT (fp <= fn );
}
