#include "Parallel_LU.hpp"
/** =========================================================================  /
 * =======================  paru_print   ====================================  /
 * ========================================================================== */
void paru_print_element (paru_matrix *paruMatInfo, Int e){
    DEBUGLEVEL(0);
    // print out contribution blocks
    Element **elementList; 
    elementList = paruMatInfo->elementList;
    Element *curEl = elementList[e];

    Int morign = paruMatInfo->m;
    Int nf = paruMatInfo->LUsym->nf;

    if ( e > morign + nf +1){
        printf("Element %ld is out of range; just %ld elements \n", 
                e,  morign + nf +1);
        return;
    }

    if (curEl == NULL){
        printf("Element %ld is empty\n",e );
        return;
    }

    Int m,n;
    m = curEl->nrows;
    n = curEl->ncols;
   
    Int *el_colrowIndex = colIndex_pointer (curEl);
    double *el_colrowNum = numeric_pointer (curEl);

    PRLEVEL (1, ("el_colrowIndex =%p, el_colrowNum = %p \n", 
                el_colrowIndex, el_colrowNum));

    printf("\n"); 
    printf("Element %ld is %ld x %ld:\n", e, m, n);


    printf("\t"); 
    for (int j = 0; j < n; j++) 
        printf("%ld\t", el_colrowIndex [j] );
    printf("\n"); 
    for (int i = 0; i < m; i++) {
        printf("%ld\t",el_colrowIndex [n+i] );
        for (int j = 0; j < n; j++) {
            double value =  el_colrowNum [i*m + j];
            printf("%2.4lf\t",value );
        }
        printf("\n"); 
    }

}

void paru_print_tupleList (tupleList *listSet, Int index){
    DEBUGLEVEL(0);
    PRLEVEL (1, ("listSet =%p\n", listSet));

    if (listSet == NULL) {
       printf("Empty tuple\n"); 
       return;
    }

    tupleList cur= listSet [index];
    Int numTuple = cur.numTuple;
    Tuple *l = cur.list;

    printf(" There are %ld tuples in this list:\n", numTuple);
    for (Int i = 0; i < numTuple; i++) {
       Tuple curTpl = l [i];
        printf(" (%ld,%ld)", curTpl.e, curTpl.f);
    }
    printf("\n"); 
}
