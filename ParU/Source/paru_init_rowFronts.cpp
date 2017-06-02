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
 *  be implemnted in different files and will be reused for other fronts    */

#include "Parallel_LU.hpp"
void paru_init_rowFronts(
        // inputs, not modified
        cholmod_sparse *A,
        //symbolic analysis
        paru_symbolic *LUsym,
        // workspace and parameters
        cholmod_common *cc
        ){

    paru_matrix* Amat;
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
    Element** elementList; Int nf=LUsym->nf;
    elementList=Amat->elementList=
        (Element**) paralloc (1, (m+nf+1)*sizeof(Element), cc);


    for(Int i=0; i<m ; i++){
        Int e= LUsym->row2atree[i]; //element number in augmented tree
        Int nrows=1,ncols; // Initializing number of columns of current element
        /*TODO*/
        //INITIALIZE NROWS and numberCols
        Element* curEl= elementList[e]=
           (Element*) paralloc (1,
                   sizeof(Element)+nrows+ncols+ nrows*ncols, cc);
        curEl->nrowsleft=curEl->nrows=nrows;
        curEl->ncolsleft=curEl->ncols=ncols;

        double *colrowIndex = (double *)(curEl+1);
        for (Int j = 0;  j<ncols ; j++) {
            /*TODO Initalizing indices and Update tuple list*/
            //initializae colrowIndex[j]
        }
        colrowIndex[j++]=i;  //initializing row /TODO update rowtuplelist
        double *numericIndex= (double *)(curEl+1)+ncols+1;
         for (Int j = 0;  j<ncols ; j++) {
            /*TODO Numerics*/
            //numericIndex[j]
        }
        //TODO assemble nth row
        //

        //Add nth row to Row list and update Col list
    }

}

