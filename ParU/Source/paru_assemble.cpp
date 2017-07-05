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

    PRLEVEL (0, ("Assemble Front %ld\n", f));
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
        PRLEVEL (1, ("c =%ld numTuple = %ld\n", c, numTuple));
        for (Int i = 0; i < numTuple; i++){
            Tuple curTpl = l [i];
            Int e = curTpl.e;
            Element **elementList = paruMatInfo->elementList;
            Element *curEl = elementList[e];
            Int mEl = curEl->nrows;
            Int *el_colrowIndex = (Int*)(curEl+1);     // pointers to element index 
            PRLEVEL (1, ("mEl =%ld rowsP=%ld\n", mEl, rowsP));
            Int rS;
            for (Int rEl = 0; rEl < mEl; rEl++){
                for (rS = 0; rS < rowsP; rS++){
                    PRLEVEL (1, ("rS =%ld rEl=%ld\n", rS, rEl));
                    if (el_colrowIndex [rEl] == rowSet [rS])
                        break; 
                }
                if ( rS == rowsP){ // count the new row
                    if (rowsP >= setSize){
                        PRLEVEL (0, ("rowsP =%ld setSize=%ld\n", 
                                    rowsP, setSize));
                        PRLEVEL (0, ("setSize*2 =%ld\n", setSize*2));
                        Int *newSet= (Int*) paru_realloc (setSize*2, 
                                sizeof (Int), rowSet, &setSize, cc);

                        if(newSet== NULL){
                            printf("Error in allocating memory for rows\n");
                            return;
                        }
                        rowSet = newSet;
                    }
                    rowSet [rowsP++] = el_colrowIndex [rEl] ;
                }
            }
        }
    }



    double *pF = (double*) paru_calloc (rowsP*fp, sizeof (double), cc);
    /*! TODO: Here:     */
    /* assembling the pivotal part of the front */

#ifndef NDEBUG
    PRLEVEL (0, ("There are %ld rows in this front: ", rowsP));
    for (Int i = 0; i < rowsP; i++)
        PRLEVEL (0, (" %ld", rowSet[i]));
    PRLEVEL (0, ("\n"));
#endif 

    paru_free (setSize, sizeof (Int), rowSet, cc);
    paru_free (rowsP*fp, sizeof (Int), pF, cc);
}
