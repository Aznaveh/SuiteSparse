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
 *                      Pinit     Ps = (compute here) newRofS (paru_write)
 *                    --------> -------->
 *                   A         S           LU
 *                    <-------   <-------
 *         (paru_analyze)Pinv     oldRofS (paru_write)
 *
 *
 *  We need these permuataions for compuing Ax = b
 *        x = b (p)
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
    // Int n = LUsym->n;
    Int n1 = LUsym->n1;  // row+col singletons

    // paru_fac *LUs = paruMatInfo->partial_LUs;
    // paru_fac *Us = paruMatInfo->partial_Us;
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
    Int *oldRofS = (Int *)paru_alloc(m, sizeof(Int));
    Int *newRofS = (Int *)paru_alloc(m, sizeof(Int));
    if (oldRofS == NULL || newRofS == NULL)
    {
        printf("memory problem inside perm in debug mode\n");
        return;
    }
#endif

    Int ip = 0;  // number of rows seen so far
    // TODO: singletons shouldn't affect this.
    // for (Int k = 0; k < n1; k++)
    //    // first singletons
    //    Pfin[ip++] = Pinit[k];

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
            oldRofS[ip] = frowList[k];  // computing permutation for S
#endif
            Ps[frowList[k]] = ip;
            Pfin[ip++] = Pinit[frowList[k]];
        }
    }

#ifndef NDEBUG
    //-------- computing the direct permutation of S only for debug
    for (Int k = 0; k < m - n1; k++)
    {
        // Inv permutation for S Pinv[i] = k;
        newRofS[oldRofS[k]] = k;
        ASSERT(Ps[oldRofS[k]] == newRofS[oldRofS[k]]);
    }

    paru_free(m, sizeof(Int), oldRofS);
    paru_free(m, sizeof(Int), newRofS);
    PRLEVEL(1, ("%% Final row permutaion is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %ld, ", Pfin[k]));
    }
    PRLEVEL(1, (" \n"));
#endif
}
///////////////apply perm x = b(p) /////////////////////////////////////////////
Int paru_apply_perm(const Int *p, const double *b, double *x, Int m)
{
    DEBUGLEVEL(1);
    if (!x || !b) return (0);

#ifndef NDEBUG
    PRLEVEL(1, ("%% Inside apply permutaion p is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %ld, ", p[k]));
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
        Int j = p[k];  // k-new and j-old; P(new) = old
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
Int paru_apply_inv_perm(const Int *p, const double *b, double *x, Int m)
{
    DEBUGLEVEL(1);
    if (!x || !b) return (0);
#ifndef NDEBUG
    PRLEVEL(1, ("%% Inside apply inv permutaion p is:\n%%"));
    for (Int k = 0; k < m; k++)
    {
        PRLEVEL(1, (" %ld, ", p[k]));
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
        Int j = p[k];  // k-new and j-old; P(new) = old
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
