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
    DEBUGLEVEL(0);
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
        PRLEVEL(PR, ("(%ld)%ld ", ip - 1, Pfin[ip - 1]));
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
            Ps[frowList[k]] = ip - n1;
            Pfin[ip++] = Pinit[frowList[k] + n1];
            PRLEVEL(PR, ("(%ld)%ld\n ", ip - 1, Pfin[ip - 1]));
        }
    }
    PRLEVEL(PR, ("\n"));

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%% Final Ps:\n%%"));
    for (Int k = 0; k < m - n1; k++)
    {
        PRLEVEL(PR, (" %ld, ", Ps[k]));
    }
    PRLEVEL(PR, (" \n"));
    PRLEVEL(PR, ("%% n1=%ld Final row permutaion is:\n%%", n1));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(PR, (" %ld, ", Pfin[k]));
    }
    PRLEVEL(PR, (" \n"));
#endif
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

///////////////apply perm and scale x = sb(P) //////////////////////////////////
Int paru_apply_perm_scale(const Int *P, const double *s, const double *b,
                          double *x, Int m)
{
    DEBUGLEVEL(0);
    if (!x || !b) return (0);

#ifndef NDEBUG
    PRLEVEL(1, ("%% Inside apply permutaion and scale P is:\n%%"));
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

    PRLEVEL(1, ("%% and s is\n%%"));

    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %lf, ", s[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
    for (Int k = 0; k < m; k++)
    {
        Int j = P[k];  // k-new and j-old; P(new) = old

#ifndef NDEBUG
        PRLEVEL(1, ("x[%ld]= %lf ", k, x[k]));
        if (s != NULL) PRLEVEL(1, ("s[%ld]=%lf, ", j, s[j]));
#endif
        x[k] = (s == NULL) ? b[j] : b[j] / s[j];
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
