/* =========================================================================   /
 * ============================== paru_perm  ===============================   /
 * =========================================================================   /
 * @brief Computing and saving row permutation. This must be doen after
 * factorization.
 *   I have this transition   A ---> S ---> LU
 *   There are both col and row permutation form A to S.
 *   However there is no column permuation from S to LU. Therefore the overall
 *   column permutaion is the same with S. (Qfill)
 *   Row permutation happens from S to LU.
 *   Row permutation and inverse permutation is computed here
 *
 *                         ------p--->
 *                         A         LU
 *                         <----q----
 *
 *                          Pfin (COMPUTED HERE)
 *                    ------------------>   
 *                      Pinit      oldRofS (paru_write)
 *                    --------> -------->   
 *                   A         S           LU
 *                    <-------   <-------  
 *         (paru_analyze)Pinv     newRofS (paru_write)
 *                
 * 
 *  We need these permuataions for compuing Ax = b
 *        x = b (p)
 *        x = L\x
 *        x = U\x
 *        b(q) = x
 *
 *
 *  This must be done after factorizing.
 *

 * @author Aznaveh
 * */
#include "paru_internal.hpp"
void paru_perm(paru_matrix *paruMatInfo)
{

    paru_symbolic *LUsym = paruMatInfo->LUsym;
    Int nf = LUsym->nf;

    Int m = LUsym->m;
    Int n = LUsym->n;
    Int n1 = LUsym->n1;  // row+col singletons

    paru_fac *LUs = paruMatInfo->partial_LUs;
    paru_fac *Us = paruMatInfo->partial_Us;
    Int *Super = LUsym->Super;

    // some working memory that is freed in this function
    Int *oldRofS = NULL;
    Int *newRofS = NULL;
    Int *Pinit = LUsym->Pinit;

    oldRofS = (Int *)paru_alloc(m, sizeof(Int));  // S -> LU P
    newRofS = (Int *)paru_alloc(m, sizeof(Int));  // Pinv of S

    if (oldRofS == NULL || newRofS == NULL)
    {
        printf("memory problem for writing into files\n");
        paru_free(m, sizeof(Int), oldRofS);
        paru_free(m, sizeof(Int), newRofS);
        return;
    }

    Int ip = 0;  // number of rows seen so far
    for (Int k = 0; k < n1; k++)
        // first singletons
        //fprintf(rowfptr, "%ld\n", Pinit[k]);

    for (Int f = 0; f < nf; f++)
    {  // rows for each front
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;
        Int *frowList = paruMatInfo->frowList[f];

        for (Int k = 0; k < fp; k++)
        {
            oldRofS[ip++] = frowList[k];  // computing permutation for S
            // P[k] = i
            //fprintf(rowfptr, "%ld\n", Pinit[frowList[k]]);
        }
    }
    //fclose(rowfptr);
    PRLEVEL(1, ("%% row permutaion DONE\n"));

    //-------- computing the direct permutation of S
    for (Int k = 0; k < m - n1; k++)
        // Inv permutation for S Pinv[i] = k;
        newRofS[oldRofS[k]] = k;


}
