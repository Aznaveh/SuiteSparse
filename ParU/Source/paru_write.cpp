/** =========================================================================  /
 * =======================  paru_write  =====================================  /
 * ========================================================================== */

#include "Parallel_LU.hpp"

/*! @brief Writing the results into a file
 *    it must be called after ther result are computed
 *  @author Aznaveh
 */
void paru_write( paru_matrix *paruMatInfo, int scale,
        char* id,  cholmod_common *cc){

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

    char default_name[] = "0";
    char *name;
    if (id)
        name = id;
    else 
        name = default_name;
    


//    char dpath[] = "/export/scratch/multifrontal/aznaveh/myResults/";
    char dpath[] = "/users/aznaveh/SuiteSparse/ParU/Demo/Res/";


    //-------------------- writing column permutation to a file
    {
        FILE *colfptr;

        char fname [100] = "";
        strcat (fname,dpath);
        strcat (fname,name);
        strcat (fname,"_col.txt");
        colfptr = ( fopen(fname,"w"));

    //colfptr = ( fopen("/users/aznaveh/SuiteSparse/ParU/Demo/%s_col.txt","w"));
        //    colfptr = ( fopen("./col.txt","w"));
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
    }
    //--------------------


    Int *oldRofS= (Int*) paru_alloc ( m, sizeof (Int), cc); // S -> LU
    Int *PofA= (Int*) paru_alloc ( m, sizeof (Int), cc);    // inv(PLinv)
    //-------------------- writing row permutation to a file
    {
        FILE *rowfptr;
        char fname [100] = "";
        strcat (fname,dpath);
        strcat (fname,name);
        strcat (fname,"_row.txt");
        rowfptr = ( fopen(fname,"w"));


        // rowfptr = ( fopen("/users/aznaveh/SuiteSparse/ParU/Demo/row.txt","w"));
        if (rowfptr == NULL ){
            printf ("Error in opening a file");
            return;
        }
        fprintf (rowfptr, "%%rows\n");
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
    }
    //--------------------


    //-------- computing the direct permutation
    Int *newRofS= (Int*) paru_alloc ( m, sizeof (Int), cc);
    for(Int k = 0; k < m ; k++){
        newRofS[oldRofS[k]] = k;
    }

    //--------------------


    //-------------------- writing row scales to a file
    if (scale) {
        double *scale_row = paruMatInfo->scale_row;
        FILE *scalefptr;
        char fname [100] = "";
        strcat (fname,dpath);
        strcat (fname,name);
        strcat (fname,"_scale.txt");
        scalefptr = ( fopen(fname,"w"));


        if (scalefptr== NULL ){
            printf ("Error in opening a file");
            return;
        }
        for(Int row = 0; row < m ; row++)  
            fprintf (scalefptr, "%16g\n",1/scale_row[row]);
        fclose(scalefptr);
    }
    //--------------------

    //-------------------- writing info to a file
    {
        FILE *infofptr;
        char fname [100] = "";
        strcat (fname,dpath);
        strcat (fname,name);
        strcat (fname,"_info.txt");
        infofptr = ( fopen(fname,"w"));


        if (infofptr == NULL ){
            printf ("Error in opening a file");
            return;
        }
        fprintf (infofptr, "%16g\n",paruMatInfo->time);
        fclose(infofptr);
    }
    //--------------------

    //-------------------- writing results to a file
    FILE *LUfptr;
    char fname [100] = "";
    strcat (fname,dpath);
    strcat (fname,name);
    strcat (fname,"_LU.txt");
    LUfptr = ( fopen(fname,"w"));


    //LUfptr = ( fopen("/users/aznaveh/SuiteSparse/ParU/Demo/out.txt","w"));
    if (LUfptr == NULL ){
        printf ("Error in opening a file");
        return;
    }

    //computing nnnz
    Int nnz = 0;
    for(Int f = 0; f < nf ; f++){   
        Int colCount =  paruMatInfo->fcolCount[f];
        Int rowCount =  paruMatInfo->frowCount[f];
        Int col1 = Super [f];     
        Int col2 = Super [f+1];
        Int fp = col2-col1;
        nnz += fp*(rowCount + colCount);
    }

    fprintf (LUfptr, "%%%MatrixMarket matrix coordinate real general\n");
    fprintf (LUfptr, "%%-----------produced by ParU ---------------\n");
    fprintf (LUfptr, "%ld  %ld %ld\n",m ,n ,nnz );

    for(Int f = 0; f < nf ; f++){   
        Int colCount =  paruMatInfo->fcolCount[f];
        Int *fcolList =  paruMatInfo->fcolList[f];
        Int rowCount =  paruMatInfo->frowCount[f];
        Int *frowList =  paruMatInfo->frowList[f];
        Int col1 = Super [f];     
        Int col2 = Super [f+1];
        Int fp = col2-col1;

        //Printing LU part
        double *pivotalFront= LUs[f].p  ;
        PRLEVEL (0, ("%% pivotalFront =%p \n", pivotalFront));
        for (Int j = col1 ; j < col2; j++)
            for (Int i = 0; i < rowCount ; i++){
                fprintf (LUfptr, "%ld  %ld %.16g\n",
                        newRofS[frowList[i]]+1, j+1, 
                        pivotalFront[(j-col1)*rowCount+i]);
            }

#ifndef NDEBUG  // Printing the pivotal front
        Int p = 0;
        PRLEVEL (p, ("\n%%Inside paru_write Luf{%ld}= [",f+1));
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
                fprintf (LUfptr, "%ld  %ld %.16g\n", newRofS[frowList[i]]+1,
                        fcolList[j]+1, uPart[fp*j+i]);
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
