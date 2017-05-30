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

   paru_work* Amat;
   Amat = (paru_work *) paralloc (1,sizeof(paru_work),cc);
   if (Amat == NULL){   //out of memory
       return;
   }
   Int m,n;
   m= Amat->nrows= LUsym->m;   n= Amat->ncols= LUsym->n;
   Amat->Row_list =(listEl*) paralloc (1,m*sizeof(listEl),cc);
   Amat->Col_list =(listEl*) paralloc (1,n*sizeof(listEl),cc);

}
