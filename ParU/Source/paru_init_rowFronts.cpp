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
   Amat->Row_list =(listEl*) paralloc (slackRow, m*sizeof(listEl), cc);
   Amat->Col_list =(listEl*) paralloc (slackCol, n*sizeof(listEl), cc);

   for(Int i=0; i<m ; i++){
       Int e= LUsym->row2atree[i]; //element number in augmented tree
       Int nrows,ncols;
    //   (Element*) paralloc (1,sizeof(Element)+nrows*ncols, cc);
       //assemble nth row
       //Add nth row to Row list and update Col list
   }

}
