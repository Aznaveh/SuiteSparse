#include "Parallel_LU.hpp"
/** =========================================================================  /
 * =======================  paru_print   ====================================  /
 * ========================================================================== */
void paru_print_element (paru_matrix *paruMatInfo, Int e){
    // print out contribution blocks
    Element **elementList; 
    elementList = paruMatInfo->elementList;
    Element *curEl = elementList[e];

    Int morign = paruMatInfo->m;
    Int nf = paruMatInfo->LUsym->nf;

    if ( e > morign + nf +1){
        printf("%% Element %ld is out of range; just %ld elements \n", 
                e,  morign + nf +1);
        return;
    }

    if (curEl == NULL){
        printf("%% Element %ld is empty\n",e );
        return;
    }

    Int m,n;
    m = curEl->nrows;
    n = curEl->ncols;
   
    Int *el_colIndex = colIndex_pointer (curEl);
    Int *el_rowIndex = rowIndex_pointer (curEl);

    Int *rel_col = relColInd (curEl);
    Int *rel_row = relRowInd (curEl);
    
    double *el_colrowNum = numeric_pointer (curEl);

    printf("\n"); 
    printf("%% Element %ld is %ld x %ld:\n", e, m, n);


    printf("\t"); 
//    for (int j = 0; j < n; j++) 
//        printf("%% %ld\t", rel_col[j] );
//    printf("\n\t"); 
    for (int j = 0; j < n; j++) 
        printf("%% %ld\t", el_colIndex [j] );

    printf("\n"); 
    for (int i = 0; i < m; i++) {
   //     printf("%% %ld\t %ld\t",rel_row[i], el_rowIndex [i] );
        printf("%% %ld\t", el_rowIndex [i] );
        for (int j = 0; j < n; j++) { 
            double value =  el_colrowNum [j*m + i];
            printf("%2.4lf\t",value );
        }
        printf("\n"); 
    }

}

void paru_print_tupleList (tupleList *listSet, Int index){
    DEBUGLEVEL(0);
    PRLEVEL (1, ("%% listSet =%p\n", listSet));

    if (listSet == NULL) {
        printf("%% Empty tuple\n"); 
        return;
    }

    tupleList cur= listSet [index];
    Int numTuple = cur.numTuple;
    Tuple *l = cur.list;

    printf("%% There are %ld tuples in this list:\n %%", numTuple);
    for (Int i = 0; i < numTuple; i++) {
        Tuple curTpl = l [i];
        printf(" (%ld,%ld)", curTpl.e, curTpl.f);
    }
    printf("\n"); 
}
