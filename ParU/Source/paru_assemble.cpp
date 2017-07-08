#include "Parallel_LU.hpp"

#ifndef NDEBUG  // using STL for debugging
#include <iostream>
#include <algorithm>
#include <set>
#endif

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

    /* computing number of rows, set union */
    tupleList *ColList = paruMatInfo->ColList;

#if 0
    Int rowsP = 0;
    Int setSize = fn;  // it can be initialized with fm
    Int *rowSet = (Int*) paru_alloc (setSize, sizeof (Int), cc);
#endif

    Int listP= 0;
    work_struct *Work =  paruMatInfo->Work;
    Int *isInSet = Work->all_Zero; /*  I can reduce the size with bitwise 
                                       operation. it might cause race in GPU */
    Int *rowList = Work->scratch;
    PRLEVEL (1, ("rowList(scratch)=%p isInSet(all_Zero)=%p\n", 
                rowList, isInSet));

#ifndef NDEBUG /* chekcing first part of Work to be zero */
    Int s = 0;
    for (Int i = 0; i < m; i++) s += isInSet [i];
    ASSERT (s == 0);
#endif 


#ifndef NDEBUG
    std::set<Int> stl_rowSet;
    std::set<Int>::iterator it;
#endif 
    for (Int c = col1; c < col2; c++){
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
            Int nEl = curEl->ncols;
            Int *el_colrowIndex = (Int*)(curEl+1);     // pointers to element index 
            Int *el_rowIndex = el_colrowIndex + nEl;   // pointers to row indices
            PRLEVEL (1, ("element= %ld  mEl =%ld \n",e, mEl));
            Int rS;
            for (Int rEl = 0; rEl < mEl; rEl++){
                Int curRow = el_rowIndex [rEl]; 
                PRLEVEL (1, ("curRow =%ld\n", curRow));
                ASSERT (curRow < m ) ;
#ifndef NDEBUG
                stl_rowSet.insert (curRow );
#endif
                PRLEVEL (1, ("%p ---> isInSet [%ld]=%ld\n", 
                            isInSet+curRow, curRow, isInSet[curRow]));

                if (!isInSet[curRow]){
                    PRLEVEL (1, ("curRow =%ld listP=%ld\n", curRow, listP));
                    rowList [listP] = curRow;
                    isInSet [curRow] = FLIP (listP); /* set to some nonzero
                                                        I want to use listP to
                                                        know reverse perumtation
                                                        too and listP can also 
                                                        be zero */
                    PRLEVEL (0, ("listP=%ld FLIP(listP)=%ld\n", 
                                listP, FLIP (listP) ));
                    listP++;
                }
                ASSERT (listP <= m); 

#if 0
                for (rS = 0; rS < rowsP; rS++){
                    PRLEVEL (1, ("rS =%ld rEl=%ld\n", rS, rEl));
                    if (curRow == rowSet [rS])
                        break; 
                }
                if ( rS == rowsP){ // count the new row
                    if (rowsP >= setSize){
                        PRLEVEL (1, ("rowsP =%ld setSize=%ld\n", 
                                    rowsP, setSize));
                        PRLEVEL (1, ("setSize*2 =%ld\n", setSize*2));
                        Int *newSet= (Int*) paru_realloc (setSize*2, 
                                sizeof (Int), rowSet, &setSize, cc);

                        if(newSet== NULL){
                            printf("Error in allocating memory for rows\n");
                            return;
                        }
                        rowSet = newSet;
                    }
                    rowSet [rowsP++] = curRow;
                }
#endif                
            }
        }
    }

#if 0    
    //shrinking allocated space
    rowSet= (Int*) paru_realloc (rowsP, sizeof (Int), rowSet, &setSize, cc);
#endif    

#ifndef NDEBUG /* Checking if pivotal rows are correct */
    PRLEVEL (0, ("There are %ld rows in this front: \n", listP));
    for (Int i = 0; i < listP; i++)
        PRLEVEL (0, (" %ld", rowList [i]));
    PRLEVEL (0, ("\n"));
    Int stl_size = stl_rowSet.size();
    if (listP != stl_size){
        PRLEVEL (0, ("#######################\n"));
        PRLEVEL (0, ("STL %ld:\n",stl_size));
        for (it = stl_rowSet.begin(); it != stl_rowSet.end(); it++)
            PRLEVEL (0, (" %ld", *it));
        PRLEVEL (0, ("\nMy Set %ld:\n",listP));
        for (Int i = 0; i < listP; i++)
            PRLEVEL (1, (" %ld", rowList [i]));
        PRLEVEL (0, ("\n"));
    }
    ASSERT (listP == stl_size );

#endif 

    /* clearing W for next iteration */
    for (Int i = 0; i < listP; i++){
        Int curRow = rowList [i];
        ASSERT (curRow < m );
        ASSERT (isInSet [curRow] != 0);
        isInSet  [curRow] = 0;
    }

#ifndef NDEBUG /* chekcing isInSet to be zero */
    s = 0;
    for (Int i = 0; i < m; i++) s+=isInSet [i];
    ASSERT (s == 0);
#endif 



    double *pF = (double*) paru_calloc (listP*fp, sizeof (double), cc);
    /*! TODO: Here:     */
    /* assembling the pivotal part of the front */




    paru_free (listP*fp, sizeof (Int), pF, cc);

#if 0
    paru_free (setSize, sizeof (Int), rowSet, cc);
#endif
}
