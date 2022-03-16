////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_write /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*! @brief Writing the results into a file
 *    it must be called after ther result are computed
 *  @author Aznaveh
 */
#include "paru_internal.hpp"
void paru_write(ParU_Numeric *Num, int scale, char *id)
{
    DEBUGLEVEL(0);
    PRLEVEL(1, ("%% Start Writing\n"));
    ParU_Symbolic *Sym = Num->Sym;
    Int nf = Sym->nf;

    Int m = Sym->m;
    Int n = Sym->n;
    Int n1 = Sym->n1;  // row+col singletons

    Int *Qfill = Sym->Qfill;

    ParU_Factors *LUs = Num->partial_LUs;
    ParU_Factors *Us = Num->partial_Us;
    Int *Super = Sym->Super;

    char default_name[] = "0";
    char *name;
    if (id)
        name = id;
    else
        name = default_name;

    char dpath[] = "../Demo/Res/";

    //-------------------- writing column permutation to a file
    {
        FILE *colfptr;

        char fname[100] = "";
        strcat(fname, dpath);
        strcat(fname, name);
        strcat(fname, "_col.txt");
        colfptr = (fopen(fname, "w"));

        if (colfptr == NULL)
        {
            printf("Error in making %s to write the results!\n", fname);
            return;
        }
        fprintf(colfptr, "%%cols\n");

        for (Int col = 0; col < n; col++)
        {                                      // for each column of A(:,Qfill)
            Int j = Qfill ? Qfill[col] : col;  // col of S is column j of A
            fprintf(colfptr, "%ld\n", j);
        }

        fclose(colfptr);
        PRLEVEL(1, ("%% column permutaion DONE\n"));
    }
    //--------------------

    //--------------------computing and  writing row permutation to a file

    // some working memory that is freed in this function
    Int *oldRofS = NULL;
    Int *newRofS = NULL;
    Int *Pinit = Sym->Pinit;

    oldRofS = (Int *)paru_alloc(m, sizeof(Int));  // S -> LU P
    newRofS = (Int *)paru_alloc(m, sizeof(Int));  // Pinv of S

    if (oldRofS == NULL || newRofS == NULL)
    {
        printf("memory problem for writing into files\n");
        paru_free(m, sizeof(Int), oldRofS);
        paru_free(m, sizeof(Int), newRofS);
        return;
    }

    //   I have this transition   A ---> S ---> LU
    //   Mostly in the code I have rows of S.
    //   It is why I compute oldRofS and newRofS
    //              oldRofS
    //             -------->
    //           S           LU
    //             <-------
    //              newRofS
    //

    {
        FILE *rowfptr;
        char fname[100] = "";
        strcat(fname, dpath);
        strcat(fname, name);
        strcat(fname, "_row.txt");
        rowfptr = (fopen(fname, "w"));

        if (rowfptr == NULL)
        {
            printf("Error in opening a file");
            return;
        }
        fprintf(rowfptr, "%%rows\n");

        Int ip = 0;  // number of rows seen so far
        for (Int k = 0; k < n1; k++)
            // first singletons
            fprintf(rowfptr, "%ld\n", Pinit[k]);

        for (Int f = 0; f < nf; f++)
        {  // rows for each front
            Int col1 = Super[f];
            Int col2 = Super[f + 1];
            Int fp = col2 - col1;
            Int *frowList = Num->frowList[f];

            for (Int k = 0; k < fp; k++)
            {
                oldRofS[ip++] = frowList[k];  // computing permutation for S
                // P[k] = i
                fprintf(rowfptr, "%ld\n", Pinit[frowList[k]]);
            }
        }
        fclose(rowfptr);
        PRLEVEL(1, ("%% row permutaion DONE\n"));
    }
    //--------------------

    //-------- computing the direct permutation of S
    for (Int k = 0; k < m - n1; k++)
        // Inv permutation for S Pinv[i] = k;
        newRofS[oldRofS[k]] = k;

    //--------------------

    //-------------------- writing row scales to a file
    if (scale)
    {
        double *scale_row = Sym->scale_row;
        FILE *scalefptr;
        char fname[100] = "";
        strcat(fname, dpath);
        strcat(fname, name);
        strcat(fname, "_scale.txt");
        scalefptr = (fopen(fname, "w"));

        if (scalefptr == NULL)
        {
            printf("Error in opening a file");
            return;
        }
        for (Int row = 0; row < m; row++)
            fprintf(scalefptr, "%.17g\n", scale_row[row]);
        fclose(scalefptr);
    }
    //--------------------

    //-------------------- writing info to a file
    {
        FILE *infofptr;
        char fname[100] = "";
        strcat(fname, dpath);
        strcat(fname, name);
        strcat(fname, "_info.txt");
        infofptr = (fopen(fname, "w"));

        if (infofptr == NULL)
        {
            printf("Error in opening a file");
            return;
        }
        fprintf(infofptr, "%.17g\n", Num->my_time + Sym->my_time);
        fprintf(infofptr, "%.17g\n", Num->umf_time);

#ifdef COUNT_FLOPS
        fprintf(infofptr, "%.17g\n", Num->flp_cnt_dgemm);
        fprintf(infofptr, "%.17g\n", Num->flp_cnt_trsm);
        fprintf(infofptr, "%.17g\n", Num->flp_cnt_dger);
        fprintf(infofptr, "%.17g\n", Num->flp_cnt_real_dgemm);
#endif
        fclose(infofptr);
    }
    //--------------------

    //-------------------- writing results to a file
    FILE *LUfptr;
    char fname[100] = "";
    strcat(fname, dpath);
    strcat(fname, name);
    strcat(fname, "_LU.txt");
    LUfptr = (fopen(fname, "w"));

    if (LUfptr == NULL)
    {
        printf("Error in opening a file");
        return;
    }

    // computing nnz of factorized S
    Int nnz = 0;
    for (Int f = 0; f < nf; f++)
    {
        Int colCount = Num->fcolCount[f];
        Int rowCount = Num->frowCount[f];
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;
        // nnz += fp * (rowCount + colCount);

        double *pivotalFront = LUs[f].p;
        double *uPart = Us[f].p;
        for (Int j = col1; j < col2; j++)
        {
            for (Int i = 0; i < rowCount; i++)
            {
                if (pivotalFront[(j - col1) * rowCount + i] != 0.0) nnz++;
            }
        }

        for (Int j = 0; j < colCount; j++)
            for (Int i = 0; i < fp; i++)
            {
                {
                    if (uPart[fp * j + i] != 0.0) nnz++;
                }
            }
    }
    nnz += Sym->anz - Sym->snz;  // adding singletons

    fprintf(LUfptr, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(LUfptr, "%%-----------produced by ParU ---------------\n");
    fprintf(LUfptr, "%ld  %ld %ld\n", m, n, nnz);

    // writing the singletons
    // TODO (I haven't add that part yet)

    // writing the L and U factors
    for (Int f = 0; f < nf; f++)
    {
        Int colCount = Num->fcolCount[f];
        Int *fcolList = Num->fcolList[f];
        Int rowCount = Num->frowCount[f];
        Int *frowList = Num->frowList[f];
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;

        // Printing LU part
        double *pivotalFront = LUs[f].p;
        PRLEVEL(1, ("%% pivotalFront =%p \n", pivotalFront));
        for (Int j = col1; j < col2; j++)
            for (Int i = 0; i < rowCount; i++)
            {
                if (pivotalFront[(j - col1) * rowCount + i] != 0.0)
                    fprintf(LUfptr, "%ld  %ld %.17g\n",
                            newRofS[frowList[i]] + n1 + 1, j + n1 + 1,
                            pivotalFront[(j - col1) * rowCount + i]);
            }

#ifndef NDEBUG  // Printing the pivotal front
        Int p = 1;
        PRLEVEL(p, ("\n%%Inside paru_write Luf{%ld}= [", f + 1));
        for (Int r = 0; r < rowCount; r++)
        {
            PRLEVEL(p, (" "));
            for (Int c = col1; c < col2; c++)
                PRLEVEL(p,
                        (" %.17g ", pivotalFront[(c - col1) * rowCount + r]));
            PRLEVEL(p, (";\n%% "));
        }
        PRLEVEL(p, (";]\n"));
#endif

        // Printing U part
        double *uPart = Us[f].p;
        for (Int j = 0; j < colCount; j++)
            for (Int i = 0; i < fp; i++)
            {
                if (uPart[fp * j + i] != 0.0)
                    fprintf(LUfptr, "%ld  %ld %.17g\n",
                            newRofS[frowList[i]] + n1 + 1, fcolList[j] + n1 + 1,
                            uPart[fp * j + i]);
            }
#ifndef NDEBUG  // Printing the  U part
        p = 1;
        PRLEVEL(p, ("\n"));
        PRLEVEL(p, ("%% fp = %ld, colCount = %ld\n", fp, colCount));
        if (colCount != 0)
        {
            for (Int i = 0; i < fp; i++)
            {
                for (Int j = 0; j < colCount; j++)
                    PRLEVEL(p, (" %2.5lf\t", uPart[j * fp + i]));
                PRLEVEL(p, (";\n  %% "));
            }

            PRLEVEL(p, ("\n"));
        }
#endif
    }

    fclose(LUfptr);

    paru_free(m, sizeof(Int), oldRofS);
    paru_free(m, sizeof(Int), newRofS);
}
