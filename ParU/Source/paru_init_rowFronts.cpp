/* ========================================================================== */
/* =======================  paru_init_rowFronts  ============================ */
/* ========================================================================== */
/* Assembling row fronts; other fronts will be assembled from their childs.
 *      Initializing Row and column list: 
 *          allocating memory and updating lists and initializing matrix 
 *          structre
 *      Assemble each front:
 *          Adding numerical values, allocating data
 *          updating the list
 * 
 * Updating row and column list and allocating memory for fronts will
 *  be implemnted in different files and will be reused for other fronts      */

#include "Parallel_LU.hpp"
void paru_init_rowFronts(
        // inputs, not modified
        cholmod_sparse *A,
        //symbolic analysis
        paru_symbolic *LUsym,
        // workspace and parameters
        cholmod_common *cc
){

    paru_matrix *Amat;
    Amat = (paru_matrix*) paralloc (1,sizeof(paru_matrix),cc);
    if (Amat == NULL){   //out of memory
        return;
    }

    Int m,n;     Int slackRow=2,slackCol=2;
    m= Amat->m= LUsym->m;   n= Amat->n= LUsym->n; 
    tupleList *RowList,*ColList;
    RowList=Amat->RowList =
        (tupleList*) paralloc (slackRow, m*sizeof(tupleList), cc);
    ColList=Amat->ColList =
        (tupleList*) paralloc (slackCol, n*sizeof(tupleList), cc);
    Element **elementList; Int nf = LUsym->nf;
    elementList = Amat->elementList = // Initialize with NULL
        (Element**) parcalloc (1, (m+nf+1)*sizeof(Element), cc);

    //Double check with Dr Davis
    cholmod_sparse *C = cholmod_transpose (A, 2, cc);

    double *Cx;
    Int *Cp, *Ci, *Cnz;
    Int p, pend, ncolC, xtype;
    ncolC = C->ncol;        Cp = (Int*) C->p;
    Ci = (Int*) C->i;       Cx = (double*) C->x;
    Cnz = (Int*) C->nz;     xtype = C->xtype;
    ASSERT(ncolC == m);
    if(xtype != CHOLMOD_REAL){
        printf("Error: Input matrix is not double\n");
        return; // Error: Just working with real for now
    }

    for(Int i = 0; i < m ; i++){
        Int e = LUsym->row2atree[i]; //element number in augmented tree
        Int nrows = 1,
            ncols = Cnz[i]; 

        Element *curEl = elementList[e] =
            (Element*) paralloc(1,
                    sizeof(Element)+sizeof(Int)*(nrows+ncols)+ 
                    sizeof(double)*nrows*ncols, cc);
        curEl->nrowsleft = curEl->nrows = nrows;
        curEl->ncolsleft = curEl->ncols = ncols;

        Int *colrowIndex = (Int*)(curEl+1);
        double *colrowNum = (double*)(colrowIndex+nrows+ncols);
        Int j = 0;
        /* TODO: Initializing tuple list */
        p = Cp [i];
        pend = p + Cnz[i];
        for ( ; p < pend ; p++){
            colrowIndex[j] = Ci[p];
            colrowNum[j++] =   Cx[p];
            /* TODO Add col tuples:  <06-06-17, Me> */
        }
        colrowIndex[j++] = i;  //initializing row one item  
        /*! TODO: add Row tuple	 */

    }
}
