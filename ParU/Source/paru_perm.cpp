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
 *                         ------P--->
 *                         A         LU
 *                     **********    
 *                     **********     The rest is identity
 *                     ***#######    #######    
 *                     ***#######    #######    
 *                         <----q----
 *
 *                          Pfin (COMPUTED HERE)
 *                    ------------------>
 *                      Pinit     Ps = (compute here) newRofS (paru_write)
 *                    --------> -------->
 *                   A         S           LU
 *                    <-------   <-------
 *         (paru_analyze)Pinv     oldRofS (paru_write)
 *
 *
 *  We need these permuataions for compuing Ax = b
 *        x = b (P)
 *        x = L\x
 *        x = U\x
 *        b(q) = x
 *

 * @author Aznaveh
 * */
#include "paru_internal.hpp"
void paru_perm(paru_matrix *paruMatInfo)
{
    DEBUGLEVEL(1);
    paru_symbolic *LUsym = paruMatInfo->LUsym;

    if (LUsym->Pfin != NULL)  // it must have been computed
        return;
    Int nf = LUsym->nf;

    Int m = LUsym->m;

    Int *Super = LUsym->Super;

    // some working memory that is freed in this function
    Int *Pfin = NULL;
    Int *Ps = NULL;
    Int *Pinit = LUsym->Pinit;

    LUsym->Pfin = Pfin = (Int *)paru_alloc(m, sizeof(Int));
    LUsym->Ps = Ps = (Int *)paru_alloc(m, sizeof(Int));

    PRLEVEL(1, ("%% Inside Perm\n"));
    if (Pfin == NULL || Ps == NULL)
    {
        printf("memory problem inside perm\n");
        return;
    }

#ifndef NDEBUG
    //Int *oldRofS = (Int *)paru_alloc(m, sizeof(Int));
    //Int *newRofS = (Int *)paru_alloc(m, sizeof(Int));
    //if (oldRofS == NULL || newRofS == NULL)
    //{
    //    printf("memory problem inside perm in debug mode\n");
    //    return;
    //}
    Int PR = 1;

    PRLEVEL(PR, ("%% Initial row permutaion is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(PR, (" %ld, ", Pinit[k]));
    }
    PRLEVEL(PR, (" \n"));
#endif

    Int n1 = LUsym->n1;  // row+col singletons 
    Int ip = 0;          // number of rows seen so far
    PRLEVEL(PR, ("%% singlton part"));
    for (Int k = 0; k < n1; k++)
    {  // first singletons
        Pfin[ip++] = Pinit[k];
        PRLEVEL(PR, ("(%ld)%ld ", ip-1, Pfin[ip-1]));
    }
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("%% the rest\n"));
    for (Int f = 0; f < nf; f++)
    {  // rows for each front
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = col2 - col1;
        Int *frowList = paruMatInfo->frowList[f];

        for (Int k = 0; k < fp; k++)
        {
            // P[k] = i
#ifndef NDEBUG
     //       oldRofS[ip - n1] = frowList[k];  // computing permutation for S
     //       PRLEVEL(1, ("%% frowList[%ld]= %ld-", k, frowList[k]));
#endif
            Ps[frowList[k]] = ip - n1;
            Pfin[ip++] = Pinit[frowList[k]+n1];
            PRLEVEL(PR, ("(%ld)%ld\n ", ip-1, Pfin[ip-1]));
        }
    }
    PRLEVEL(PR, ("\n"));

#ifndef NDEBUG
    //-------- computing the direct permutation of S only for debug
//    for (Int k = 0; k < m - n1; k++)
//    {
//        // Inv permutation for S Pinv[i] = k;
//        newRofS[oldRofS[k]] = k;
//        ASSERT(Ps[oldRofS[k]] == newRofS[oldRofS[k]]);
//    }
//
//    paru_free(m, sizeof(Int), oldRofS);
//    paru_free(m, sizeof(Int), newRofS);
//


    PRLEVEL(PR, ("%% Final Ps:\n%%"));
    for (Int k = 0; k < m-n1; k++)
    {
        PRLEVEL(PR, (" %ld, ", Ps[k]));
    }
    PRLEVEL(PR, (" \n"));
    PR = 1;
    PRLEVEL(PR, ("%% n1=%ld Final row permutaion is:\n%%",n1));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(PR, (" %ld, ", Pfin[k]));
    }
    PRLEVEL(PR, (" \n"));
#endif
}
///////////////apply perm x = b(P) ///////////////////////////////////
Int paru_apply_perm(const Int *P, const double *b, double *x, Int m)
{
    DEBUGLEVEL(0);
    if (!x || !b) return (0);

#ifndef NDEBUG
    PRLEVEL(1, ("%% Inside apply permutaion P is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %ld, ", P[k]));
    }
    PRLEVEL(1, (" \n"));

    PRLEVEL(1, ("%% and b is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", b[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    for (Int k = 0; k < m; k++)
    {
        Int j = P[k];  // k-new and j-old; P(new) = old
        x[k] = b[j];
    }

#ifndef NDEBUG
    PRLEVEL(1, ("%% after applying permutaion x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
///////////////apply inverse perm x = b(pinv) //////////////////////////////////
Int paru_apply_inv_perm(const Int *P, const double *b, double *x, Int m)
{
    DEBUGLEVEL(0);
    if (!x || !b) return (0);
#ifndef NDEBUG
    PRLEVEL(1, ("%% Inside apply inv permutaion P is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %ld, ", P[k]));
    }
    PRLEVEL(1, (" \n"));

    PRLEVEL(1, ("%% before applying inverse permutaion b is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", b[k]));
    }
    PRLEVEL(1, (" \n"));
#endif

    for (Int k = 0; k < m; k++)
    {
        Int j = P[k];  // k-new and j-old; P(new) = old
        x[j] = b[k];   // Pinv(old) = new
    }

#ifndef NDEBUG
    PRLEVEL(1, ("%% after applying inverse permutaion x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.8lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}

///////////////apply scale x= s.b  /////////////////////////////////////////////
Int paru_apply_scale(const double *s, const Int *Ps, double *x, Int m, Int n1)
{
    DEBUGLEVEL(1);
#ifndef NDEBUG
    PRLEVEL(1, ("%% before applying scale x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.2lf, ", x[k]));
    }
    PRLEVEL(1, ("\n%% and s is\n"));

    for (Int k = n1; k < m; k++)
    {
        PRLEVEL(1, (" %lf, ", s[k-n1]));
    }
    PRLEVEL(1, (" \n"));
#endif

    if (!x || !s) return (0);

    for (Int k = n1; k < m; k++)
    {
        x[k] = x[k] / s[Ps[k - n1]];
    }

#ifndef NDEBUG
    PRLEVEL(1, ("%% after applying scale x is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %.8lf, ", x[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    return (1);
}
