/** =========================================================================  /
 * =======================  paru_write  =====================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

/*! @brief Writing the results into a file
 *    it must be called after ther result are computed
 *  @author Aznaveh
 */
void paru_write( paru_matrix *paruMatInfo, cholmod_common *cc){

    DEBUGLEVEL(0);
    paru_symbolic *LUsym = paruMatInfo-> LUsym;
    Int nf = LUsym->nf;
    Int m = LUsym->m;
    Int n = LUsym->n;
    Int *PLinv =  LUsym->PLinv;  
    Int *Qfill =  LUsym->Qfill;

    paru_fac *LUs =  paruMatInfo->partial_LUs;
    paru_fac *Us =  paruMatInfo->partial_Us;
    Int *Super = LUsym->Super;

    //-------------------- writing column permutation to a file
    FILE *colfptr;
    colfptr = ( fopen("./col.txt","w"));
    if (colfptr == NULL ){
        printf ("Error in opening a file");
        return;
    }
    for (Int col = 0 ; col < n ; col++){    // for each column of A(:,Qfill)
        Int j = Qfill ? Qfill [col] : col ; // col of S is column j of A
        fprintf (colfptr, "%ld\n", j );
    }
 
    fprintf (colfptr, "%%cols\n");
    fclose(colfptr);
    //--------------------
 
    
    //-------------------- writing row permutation to a file
    FILE *rowfptr;
    rowfptr = ( fopen("./row.txt","w"));
    if (rowfptr == NULL ){
        printf ("Error in opening a file");
        return;
    }
    fprintf (rowfptr, "%%rows\n");
    Int *oldRofS= (Int*) paru_alloc ( m, sizeof (Int), cc); // S -> LU
    Int *PofA= (Int*) paru_alloc ( m, sizeof (Int), cc);    // inv(PLinv)

    for(Int k = 0; k < m ; k++){
        PofA[PLinv[k]] = k;
    }
 
    Int ip = 0;
    for(Int f = 0; f < nf ; f++){   
        Int col1 = Super [f];     
        Int col2 = Super [f+1];
        Int fp = col2-col1;
        Int *frowList =  paruMatInfo->frowList[f];
        for (Int k = 0; k<fp ; k++){
            oldRofS[ip++] = frowList[k];
            fprintf (rowfptr, "%ld\n", PofA[frowList[k]] );
        }
    }
    paru_free ( m, sizeof (Int), PofA, cc);
    fclose(rowfptr);
    //--------------------


    
    Int *newRofS= (Int*) paru_alloc ( m, sizeof (Int), cc);

    for(Int k = 0; k < m ; k++){
        newRofS[oldRofS[k]] = k;
    }
    

    FILE *LUfptr;
    LUfptr = ( fopen("./out.txt","w"));
    if (LUfptr == NULL ){
        printf ("Error in opening a file");
        return;
    }

    fprintf (LUfptr, "%%MatrixMarket matrix coordinate real general\n");


    for(Int f = 0; f < nf ; f++){   
        Int colCount =  paruMatInfo->fcolCount[f];
        Int *fcolList =  paruMatInfo->fcolList[f];
        Int rowCount =  paruMatInfo->frowCount[f];
        Int *frowList =  paruMatInfo->frowList[f];
        Int col1 = Super [f];     
        Int col2 = Super [f+1];
        Int fp = col2-col1;

        Int nnz = 0;
        //Printing LU part
        double *pivotalFront= LUs[f].p  ;
        for (Int j = col1 ; j < col2; j++)
            for (Int i = 0; i < rowCount ; i++){
                fprintf (LUfptr, "%ld  %ld %.16g\n", newRofS[frowList[i]], j, 
                        pivotalFront[(j-col1)*rowCount+i]);
                nnz++;
            }

#ifndef NDEBUG  // Printing the pivotal front
        Int p = 1;
        PRLEVEL (p, ("\n%% Luf{%ld}= [",f+1));
        for (Int r = 0; r < rowCount; r++){
            PRLEVEL (p, (" "));
            for (Int c = col1; c < col2; c++){
                PRLEVEL (p, (" %.16g ", pivotalFront [(c-col1)*rowCount + r]));
            }
            PRLEVEL (p, (";\n  %% "));
        }
            PRLEVEL (p, (";\n"));
#endif

        //Printing U part
        double *uPart = Us[f].p  ;
        for (Int j = 0; j < colCount; j++)
            for (Int i = 0; i < fp; i++){
                fprintf (LUfptr, "%ld  %ld %.16g\n", newRofS[frowList[i]],
                        fcolList[j], uPart[fp*j+i]);
                nnz++;
            }
#ifndef NDEBUG  // Printing the  U part
        p = 1;
        PRLEVEL (p, ("\n"));
        for (Int i = 0; i < fp; i++){
            for (Int j = 0; j < colCount; j++){
                PRLEVEL (p, (" %2.5lf\t", uPart[j*fp+i]));
            }
            PRLEVEL (p, (";\n  %% "));
        }

        PRLEVEL (p, ("\n"));
#endif




    }
    fclose (LUfptr);

    paru_free ( m, sizeof (Int), oldRofS, cc);
    paru_free ( m, sizeof (Int), newRofS, cc);
}
