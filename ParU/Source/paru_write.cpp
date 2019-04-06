/** =========================================================================  /
 * =======================  paru_write  =====================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

/*! @brief Writing the results into a file
 *    it must be called after ther result are computed
 *  @author Aznaveh
 */
void paru_write( paru_matrix *paruMatInfo, cholmod_common *cc){
    paru_symbolic *LUsym = paruMatInfo-> LUsym;
    Int nf = LUsym->nf;

    FILE *fptr;
    fptr = ( fopen("./out.txt","w"));
    if (fptr == NULL ){
        printf ("Error in opening a file");
        return;
    }

    fprintf (fptr, "%%MatrixMarket matrix coordinate real general\n");

    for(Int f = 0; f < nf ; f++){   
        paru_fac *LUs =  paruMatInfo->partial_LUs;
        paru_fac *Us =  paruMatInfo->partial_Us;

        Int colCount =  paruMatInfo->fcolCount[f];
        Int *fcolList =  paruMatInfo->fcolList[f];

        Int rowCount =  paruMatInfo->frowCount[f];
        Int *frowList =  paruMatInfo->frowList[f];

        Int *Super = LUsym->Super;
        Int col1 = Super [f];     
        Int col2 = Super [f+1];
        Int fp = col2-col1;

        Int nnz = 0;
        //Printing LU part
        double *pivotalFront= LUs[f].p  ;
        for (Int j = col1 ; j < col2; j++)
            for (Int i = 0; i < rowCount ; i++){
                fprintf (fptr, "%ld  %ld %.16g\n", j, frowList[i], 
                        pivotalFront[rowCount*(j-col1)+i]);
                nnz++;
            }


        //Printing U part
        double *uPart = Us[f].p  ;
        for (Int j = 0; j < colCount; j++)
            for (Int i = 0; i < fp; i++){
                fprintf (fptr, "%ld  %ld %.16g\n", fcolList[j], frowList[i], 
                        pivotalFront[fp*j+i]);
                nnz++;
            }



    }
    fclose (fptr);

}
