#include "Parallel_LU.hpp"

/** =========================================================================  /
 * =======================  paru_assemble   =================================  /
 * ==========================================================================  /
/*! \brief assembling a front and updating correspoing elelment
 *
 *
 * \param  a list of tuples and the the tuple we want to add
 * \return 0 on sucess 
 */
void paru_assemble (
    paru_matrix *paruMatInfo,
    /* RowCol list/tuples and LUsym handle */
    Int f,
    /* front need to be assembled */
    cholmod_common *cc)

{
    DEBUGLEVEL(0);
    paru_symbolic *LUsym =  paruMatInfo->LUsym;

    Int m,n,nf;
    m = paruMatInfo-> m;
    n = paruMatInfo-> n;
    nf = paruMatInfo->LUsym->nf;
    Int *Super = LUsym->Super;
    /* ---------------------------------------------------------------------- */
    /* get the front F  */
    /* ---------------------------------------------------------------------- */

    PRLEVEL (1, ("Assemble Front %ld\n", f));
    /* pivotal columns Super [f] ... Super [f+1]-1 */
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];
    Int *Rp = LUsym->Rp;
    Int *Rj = LUsym->Rj;
    Int p1 = Rp [f];        /* Rj [p1:p2-1] = columns in F */
    Int p2 = Rp [f+1];
    Int fp = col2 - col1;   /* first fp columns are pivotal */ 
    Int fn = p2 - p1;          /* Upper bound number of columns of F */ 

    PRLEVEL (1, ("fp=%ld pivotal columns:clo1=%ld...col2=%ld\n", 
                fp, col1, col2-1));
    PRLEVEL (1, ("Upper bound number of columns: Rj[%ld]=%ld ... Rj[%ld]=%ld\n", 
                p1, Rj [p1], p2, Rj [p2-1]));
    ASSERT (fp > 0 );
    ASSERT (fp <= fn );

    tupleList *ColList = paruMatInfo->ColList;
    Int rowsP = 0;
    Int setSize = fn; /*! TODO: find a good initialize     */
    Int *rowSet = (Int*) paru_alloc (setSize, sizeof (Int), cc);
    /*! TODO: Check for memory in all places!!! DO IT    */
    for (Int c = col1; c < col2; c++){ /* computing number of rows, set union */
        tupleList *cur = &ColList[c];
        Int numTuple = cur->numTuple;
        ASSERT (numTuple >= 0);
        Tuple *l = cur->list;
        PRLEVEL (1, ("numTuple =%ld\n", numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = l [i];
            Int e = curTpl.e;
            Element **elementList = paruMatInfo->elementList;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int *el_colrowIndex = (Int*)(curEl+1);     // pointers to element index 
            PRLEVEL (1, ("Tuple %ld\n",i));
            for (Int rEl = 0; rEl < mEl; rEl++){
                Int rS;
                for (rS = 0; rS < rowsP; rS++){
                    PRLEVEL (1, (" %ld %ld\n", rowSet[rS], el_colrowIndex [rEl]));
                    if (el_colrowIndex [rEl] == rowSet [rS])
                        break; // row is already counted
                }
                if ( rS == rowsP){ // count the new row
                    PRLEVEL (1, (" rS=%ld rEl=%ld\n", rS, rEl));
                    if (setSize < rowsP)
                        rowSet [rowsP++] = el_colrowIndex [rEl] ;
                    else{//realloc
//                        Int *rowSettem = 
                    }

                }


            }
        }
    }

#ifndef NDEBUG
    PRLEVEL (0, ("There are %ld rows in this front: ", rowsP));
    for (Int i = 0; i < rowsP; i++)
        PRLEVEL (0, (" %ld", rowSet[i]));
    PRLEVEL (0, ("\n"));
#endif 

    paru_free (setSize, sizeof (Int), rowSet, cc);
}
