/** =========================================================================  /
 * =======================  paru_write  =====================================  /
 * ========================================================================== */
/*! @brief Writing the results into a file
 *    it must be called after ther result are computed
 *  @author Aznaveh
 */

#include "Parallel_LU.hpp"
void paru_write( paru_matrix *paruMatInfo, int scale, 
        char* id,  cholmod_common *cc)
{

    DEBUGLEVEL(0);
    PRLEVEL (1, ("%% Start Writing\n"));
    paru_symbolic *LUsym = paruMatInfo-> LUsym;
    Int nf = LUsym->nf;
    
    Int m = LUsym->m;
    Int n = LUsym->n;
    Int n1 = LUsym->n1; // row+col singletons

//    Int m = paruMatInfo->m;
//    Int n = paruMatInfo->n;

    Int *Pinv =  LUsym->Pinv;  
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

    //char dpath[] = "./";
    //it works in most cases!!
    //char dpath[] = "./Res/";
    //I have to use this one when I am calling from MATLAB
    char dpath[] = "/home/grads/a/aznaveh/SuiteSparse/ParU/Demo/Res/";
    // used to work on Parasol systems
    //char dpath[] = "/users/aznaveh/SuiteSparse/ParU/Demo/Res/";


    //-------------------- writing column permutation to a file
    {
        FILE *colfptr;

        char fname [100] = "";
        strcat (fname,dpath);
        strcat (fname,name);
        strcat (fname,"_col.txt");
        colfptr = ( fopen(fname,"w"));

        if (colfptr == NULL )
        {
            printf ("Error in making %s to write the results!\n", fname);
            return;
        }
        for (Int col = 0 ; col < n ; col++)
        {    // for each column of A(:,Qfill)
            Int j = Qfill ? Qfill [col] : col ; // col of S is column j of A
            fprintf (colfptr, "%ld\n", j );
        }

        fprintf (colfptr, "%%cols\n");
        fclose(colfptr);
        PRLEVEL (1, ("%% column permutaion DONE\n"));
    }
    //--------------------

    //--------------------computing and  writing row permutation to a file

    //some working memory that is freed in this function 
    Int *oldRofS = NULL; Int *PofA = NULL; Int *newRofS = NULL;    

    oldRofS = (Int*) paru_alloc ( m, sizeof (Int), cc); // S -> LU P
    PofA = (Int*) paru_alloc ( m, sizeof (Int), cc);    // P direct of A 
    newRofS = (Int*) paru_alloc ( m, sizeof (Int), cc); //Pinv of S

    if (oldRofS == NULL || PofA == NULL || newRofS == NULL)
    {
        printf ("memory problem for writing into files\n");
        paru_free  ( m, sizeof (Int), oldRofS, cc);
        paru_free ( m, sizeof (Int), PofA, cc);
        paru_free ( m, sizeof (Int), newRofS, cc);
        return;
    }

    {
        FILE *rowfptr;
        char fname [100] = "";
        strcat (fname,dpath);
        strcat (fname,name);
        strcat (fname,"_row.txt");
        rowfptr = ( fopen(fname,"w"));


        if (rowfptr == NULL )
        {
            printf ("Error in opening a file");
            return;
        }
        fprintf (rowfptr, "%%rows\n"); 

        for(Int k = 0; k < m ; k++)
        { //direct permutation from Pinv
            PofA[Pinv[k]] = k;       // actually I have it now with UMFPACK
                                     // but didn't save it
        }

        Int ip = 0; //number of rows seen so far
        for (Int k = 0; k < n1 ; k++)
         //first singletons
            fprintf (rowfptr, "%ld\n", PofA[k] );
        
        for(Int f = 0; f < nf ; f++)
        {  // rows for each front 
            Int col1 = Super [f];     
            Int col2 = Super [f+1];
            Int fp = col2-col1;
            Int *frowList =  paruMatInfo->frowList[f];


            for (Int k = 0; k < fp ; k++)
            {
                oldRofS[ip++] = frowList[k];   // computing permutation for S
                                               // P[k] = i 
                fprintf (rowfptr, "%ld\n", PofA[frowList[k]] );
            }
        }
        paru_free ( m, sizeof (Int), PofA, cc);
        fclose(rowfptr);
        PRLEVEL (1, ("%% row permutaion DONE\n"));
    }
    //--------------------


    //-------- computing the direct permutation of S
    for(Int k = 0; k < m-n1 ; k++)
     // Inv permutation for S Pinv[i] = k;
        newRofS[oldRofS[k]] = k;

    //--------------------


    //-------------------- writing row scales to a file
    if (scale) 
    {
        double *scale_row = paruMatInfo->scale_row;
        FILE *scalefptr;
        char fname [100] = "";
        strcat (fname,dpath);
        strcat (fname,name);
        strcat (fname,"_scale.txt");
        scalefptr = ( fopen(fname,"w"));


        if (scalefptr== NULL )
        {
            printf ("Error in opening a file");
            return;
        }
        for(Int row = 0; row < m ; row++)  
            fprintf (scalefptr, "%.17g\n",scale_row[row]);
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


        if (infofptr == NULL )
        {
            printf ("Error in opening a file");
            return;
        }
        fprintf (infofptr, "%.17g\n",paruMatInfo->my_time);
        fprintf (infofptr, "%.17g\n",paruMatInfo->umf_time);

#ifdef COUNT_FLOPS
        fprintf (infofptr, "%.17g\n",paruMatInfo->flp_cnt_dgemm);
        fprintf (infofptr, "%.17g\n",paruMatInfo->flp_cnt_trsm);
        fprintf (infofptr, "%.17g\n",paruMatInfo->flp_cnt_dger);
#endif
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


    if (LUfptr == NULL )
    {
        printf ("Error in opening a file");
        return;
    }

    //computing nnz of factorized S
    Int nnz = 0;
    for(Int f = 0; f < nf ; f++)
    {   
        Int colCount =  paruMatInfo->fcolCount[f];
        Int rowCount =  paruMatInfo->frowCount[f];
        Int col1 = Super [f];     
        Int col2 = Super [f+1];
        Int fp = col2-col1;
        nnz += fp*(rowCount + colCount);
    }
    nnz += LUsym->anz - LUsym->snz; //adding singletons

    fprintf (LUfptr, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf (LUfptr, "%%-----------produced by ParU ---------------\n");
    fprintf (LUfptr, "%ld  %ld %ld\n",m ,n ,nnz );

    // writing the singletons
    // TODO

    // writing the L and U factors
    for(Int f = 0; f < nf ; f++)
    {   
        Int colCount =  paruMatInfo->fcolCount[f];
        Int *fcolList =  paruMatInfo->fcolList[f];
        Int rowCount =  paruMatInfo->frowCount[f];
        Int *frowList =  paruMatInfo->frowList[f];
        Int col1 = Super [f];     
        Int col2 = Super [f+1];
        Int fp = col2-col1;

        //Printing LU part
        double *pivotalFront= LUs[f].p  ;
        PRLEVEL (1, ("%% pivotalFront =%p \n", pivotalFront));
        for (Int j = col1 ; j < col2; j++)
            for (Int i = 0; i < rowCount ; i++)
            {
                fprintf (LUfptr, "%ld  %ld %.17g\n",
                        newRofS[frowList[i]]+n1+1, j+n1+1, 
                        pivotalFront[(j-col1)*rowCount+i]);
            }

#ifndef NDEBUG  // Printing the pivotal front
        Int p = 1;
        PRLEVEL (p, ("\n%%Inside paru_write Luf{%ld}= [",f+1));
        for (Int r = 0; r < rowCount; r++)
        {
            PRLEVEL (p, (" "));
            for (Int c = col1; c < col2; c++)
                PRLEVEL (p, (" %.17g ", pivotalFront [(c-col1)*rowCount + r]));
            PRLEVEL (p, (";\n%% "));
        }
        PRLEVEL (p, (";]\n"));
#endif

        //Printing U part
        double *uPart = Us[f].p  ;
        for (Int j = 0; j < colCount; j++)
            for (Int i = 0; i < fp; i++)
            {
                fprintf (LUfptr, "%ld  %ld %.17g\n", newRofS[frowList[i]]+n1+1,
                        fcolList[j]+n1+1, uPart[fp*j+i]);
            }
#ifndef NDEBUG  // Printing the  U part
        p = 1;
        PRLEVEL (p, ("\n"));
        PRLEVEL (p, ("%% fp = %ld, colCount = %ld\n", fp, colCount));
        if (colCount != 0)
        {
            for (Int i = 0; i < fp; i++)
            {
                for (Int j = 0; j < colCount; j++)
                    PRLEVEL (p, (" %2.5lf\t", uPart[j*fp+i]));
                PRLEVEL (p, (";\n  %% "));
            }

            PRLEVEL (p, ("\n"));
        }
#endif
    }

    fclose (LUfptr);

    paru_free ( m, sizeof (Int), oldRofS, cc);
    paru_free ( m, sizeof (Int), newRofS, cc);
}
