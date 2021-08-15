/* =========================================================================   /
 * ============================== paru_analyze =============================   /
 * =========================================================================   /
 * @brief Computing etree and do the symbolic analysis. In this file I am going
 * to use umfpack symbolic analysis instaed of spqrsym. However, I will keep the
 * style of spqr mostly.
 *
 *  Relaxed amalgamation:
 *
 *                      Example: ./Matrix/b1_ss.mtx
 * original post ordered etree:
 *
 *   5(2) <-- front(number of pivotal cols)
 *   |   \
 *   0(1) 4(1)__
 *        |     \
 *        1(1)   3(1)
 *                \
 *                2(1)
 *
 * Relaxed tree:  threshold=3  front(number of pivotal cols)##oldfront
 *   3(2)##5
 *   |       \
 *   0(1)##0 2(3)##2,3,4
 *            |
 *            1(1)##1
 *
 *            0 1 2 3 4 5 -1
 *      fmap: 0 1 2 2 2 3 -1  last one is necessary for my check
 *
 *      fmap[oldf] == fmap [oldf+1]  amalgamated
 *      fmap[oldf-1] == fmap [oldf]  amalgamated  and root of subtree
 *
 *  Augmented tree creation:
 *
 *                      Example: ./Matrix/problem.mtx
 * original post ordered etree:          augmented tree: (Row elements)
 *
 *   4___                                     16_____
 *   |   \                                    |   \  (14-15)
 *   0    3__                                 3    13_____
 *        |  \                              (0-2)  |  \   \
 *        1   2                                    7   10 (11-12)
 *                                                 |    \
 *                                               (4-6) (8-9)
 *
 *      Note: the original etree use a dummy parent for all the tree(forest)
 *              the augmented tree does not
 * @author Aznaveh
 * */
#include "paru_internal.hpp"
paru_symbolic *paru_analyze(
    // inputs, not modified
    cholmod_sparse *A)
{
    DEBUGLEVEL(0);
    double my_start_time = omp_get_wtime();
    paru_symbolic *LUsym;

    LUsym = (paru_symbolic *)paru_alloc(1, sizeof(paru_symbolic));
    // ... check for LUsym NULL ...
    if (LUsym == NULL)
        // out of memory
        return NULL;

    Int anz;

    Int *Ap = (Int *)A->p;
    Int *Ai = (Int *)A->i;
    double *Ax = (double *)A->x;
    Int m = A->nrow;
    Int n = A->ncol;

    if (m != n)
    {
        printf("Input matrix is not square!\n");
        return NULL;
    }

    // Initializaing pointers with NULL; just in case for an early exit
    // not to free an uninitialized space
    LUsym->Chain_start = LUsym->Chain_maxrows = LUsym->Chain_maxcols = NULL;
    LUsym->Parent = LUsym->Super = LUsym->Child = LUsym->Childp = NULL;
    LUsym->Qfill = LUsym->Pfin = LUsym->Pinit = LUsym->Ps = NULL;
    LUsym->Sp = LUsym->Sj = LUsym->Sleft = NULL;
    LUsym->Sx = NULL;
    LUsym->Fm = LUsym->Cm = LUsym->Rj = LUsym->Rp = NULL;
    LUsym->aParent = LUsym->aChildp = LUsym->aChild = LUsym->row2atree = NULL;
    LUsym->super2atree = NULL;
    LUsym->first = NULL;
    LUsym->stree_flop_bound = LUsym->front_flop_bound = NULL;

    LUsym->ustons.Sup = NULL;
    LUsym->ustons.Suj = NULL;
    LUsym->ustons.Sux = NULL;

    LUsym->lstons.Slp = NULL;
    LUsym->lstons.Sli = NULL;
    LUsym->lstons.Slx = NULL;

    //############  Calling UMFPACK and retrieving data structure ##############

    /* ---------------------------------------------------------------------- */
    /*    The varialbes are needed for the UMFPACK symbolic analysis phase    */
    /* ---------------------------------------------------------------------- */

    Int nr, nc,  // A is nrxnc, I will use mxn; they should be the same anyway

        n1,  // The number of pivots with zero Markowitz cost.
        // Info[UMFPACK_COL_SINGLETONS]+Info[UMFPACK_ROW_SINGLETONS]
        // They apper first in the output permutations P and Q
        //
        //
        //  --- Copied from UMFPACK
        // 	---   x x x x x x x x x
        // 	---   . x x x x x x x x
        // 	---   . . x x x x x x x
        // 	---   . . . x . . . . .
        // 	---   . . . x x . . . .
        // 	---   . . . x x s s s s
        // 	---   . . . x x s s s s
        // 	---   . . . x x s s s s
        // 	---   . . . x x s s s s
        //
        //  ---   The above example has 3 column singletons (the first three
        //  ---   columns and their corresponding pivot rows) and 2 row
        //  ---   singletons.  The singletons are ordered first, because they
        //  ---   have zero Markowitz cost. The LU factorization for these first
        //  ---   five rows and columns is free - there is no work to do (except
        //  ---   to scale the pivot columns for the 2 row singletons), and no
        //  ---   fill-in occurs.  The remaining submatrix (4-by-4 in the above
        //  ---   example) has no rows or columns with degree one.  It may have
        //  ---   empty rows or columns.
        //
        //
        //        _______________
        //       |\**************r
        //       |  \************r -> UMFPACK_COL_SINGLETONS
        //       |   \***********r
        //       |    *\         |
        //       |    ***\       |
        //       |    ***xx\xxxxx|            |
        //       |    ***xxxx\xxx|            +   = n1
        //       -----ccc--------            /
        //             |
        //            UMFPACK_ROW_SINGLETONS

        nfr,  // The number of frontam matrices; nf in SPQR analysis

        nchains,  // The frontal matrices are related to one another by the
        // supernodal column elimination tree. Each nod in this tree
        // is one frontal matrix. The tree is partitioned into a set
        // of disjoint paths, and a frontal matrix chaing is one path
        // in this tree.  UMFPACK uses unifrontal technique to
        // factroize chains, with a single working array that holds
        // each frontal matrix in the chain, one at a time. nchains is
        // in the range 0 to nfr

        *Pinit,  // The inital row permutation. If P [k] = i, then this means
        // that row i is the kth row in the pre-ordered matrix.
        // For the unsymmetric strategy, P defines the row-merge
        // order. Let j be the column index of the leftmost nonzero
        // entry in row i of A*Q. The P defines a sort of the rows
        // according to this value. A row can appear earlier in this
        // ordering if it is aggressively absorbed before it can
        // become a pivot row. If P [k]= i, row i typically will not
        // be the kth pivot row.
        // For the symmetric strategy, P = Q. If no pivoting occurs
        // during numerical factorization, P[k] = i also defines the
        // fianl permutation of umfpack_*_numeric, for the symmetric
        // strategy.
        // NOTE: SPQR uses Pinv for stairecase structure and that is
        // an invert permutation that I have to compute the direct
        // permutation in paru_write.

        *Qinit,  // The inital column permutation. If Q [k] = j, then this
        // means that column j is the kth pivot column in pre-ordered
        // matrix. Q is not necessearily the same as final column
        // permutation in UMFPACK. In UMFPACK if the matrix is
        // structurally singular, and if the symmetric strategy is
        // used (or if Control [UMFPACK_FIXQ] > 0), then this Q will
        // be the same as the final column permutaion.
        // NOTE: there is no column permutation in paru. So the
        // initial permutaion would stay. SPQR uses Qfill for
        // staircase structure and that is the column permutation for
        // paru also.

        *Front_npivcol,  // size = n_col +1;  actual size = nfr+1
        // NOTE: This is not the case for SPQR
        // I think SPQR is easier:
        // Front_npivcol [f] = Super [f+1] - Super [f]
        // In the case of amalgamation of a child to the parent
        // then the number of pivotal column of the amalgamted node is:
        //  Super [Parent[f]+1] - Super [f]
        // Front_parent [nfr+1] is a place holder for columns
        // with no entries

        *Front_parent,  // size = n_col +1;  actual size = nfr+1
        // NOTE: This is not the case for SPQR
        // Parent is the one I should use instead.

        *Front_1strow,  // size = n_col +1;  actual size = nfr+1
        // Front_1strow [k] is the row index of the first row in
        // A (P,Q) whose leftmost entry is in pivot column for
        // kth front.
        // This is necessary only to properly factorize singular
        // matrices. Rows in the range Front_1strow [k] to
        // Front_1strow [k+1]-1 first become pivot row candidate
        // at the kth front. Any rows not eliminated in the kth
        // front maybe selected as pivot rows in the parent of k
        // (Front_1strow [k]) and so on up the tree.
        // Aznaveh: I am now using it at least for the rowMarks.

        *Front_leftmostdesc,  // size = n_col +1;  actual size = nfr+1
        // Aznaveh: I have a module computing leftmostdesc
        // for my augmented tree; so maybe do not need it

        *Chain_start,  // size = n_col +1;  actual size = nfr+1
        // The kth frontal matrix chain consists of frontal
        // matrices Chain_start [k] through Chain_start [k+1]-1.
        // Thus, Chain_start [0] is always 0 and
        // Chain_start[nchains] is the total number of frontal
        // matrices, nfr. For two adjacent fornts f and f+1
        // within a single chian, f+1 is always the parent of f
        // (that is, Front_parent [f] = f+1).
        //
        *Chain_maxrows,  // size = n_col +1;  actual size = nfr+1
        *Chain_maxcols;  // The kth frontal matrix chain requires a single
    // working array of dimension Chain_maxrows [k] by
    // Chain_maxcols [k], for the unifrontal technique that
    // factorizes the frontal matrix chain. Since the
    // symbolic factorization only provides

    void *Symbolic;  // Output argument in umf_dl_symbolc;
    // holds a pointer to the Symbolic object  if succesful
    // and NULL otherwise

    double status,           // Info [UMFPACK_STATUS]
        Info[UMFPACK_INFO],  // Contains statistics about the symbolic analysis

        Control[UMFPACK_CONTROL];  // it is set in umfpack_dl_defaults and
    // is used in umfpack_dl_symbolic; if
    // passed NULL it will use the defaults

    /* ---------------------------------------------------------------------- */
    /*    Setting up umfpakc symbolic analysis and do the  analysis phase     */
    /* ---------------------------------------------------------------------- */

    /* get the default control parameters */

    // Here is where the pre-ordering strategy is being chosen. I need a nested
    // dissection method here like metis:
    //      Control [UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    //      However I am using the default for now; Page 40 UMFPACK_UserGuide
    //      Page 22 UserGuide
    umfpack_dl_defaults(Control);
    Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    Control[UMFPACK_FIXQ] = -1;

#ifndef NDEBUG
    /* print the control parameters */
    Int PR = 1;
    if (PR <= 0) umfpack_dl_report_control(Control);
#endif

    /* performing the symbolic analysis */
    //  status = umfpack_dl_symbolic (m, n, Ap, Ai, Ax, &Symbolic, Control,
    //  Info);

    void *SW;
    status = umfpack_dl_azn_symbolic(m, n, Ap, Ai, Ax,
                                     NULL,   // user provided ordering
                                     FALSE,  // No user ordering
                                     NULL,   // user params
                                     &Symbolic,
                                     &SW,  // new in/out
                                     Control, Info);

    if (status < 0)
    {
        umfpack_dl_report_info(Control, Info);
        umfpack_dl_report_status(Control, status);
        printf("umfpack_dl_symbolic failed");
        umfpack_dl_azn_free_sw(&SW);
        umfpack_dl_free_symbolic(&Symbolic);
        paru_freesym(&LUsym);
        return NULL;
    }

    Int cs1 = Info[UMFPACK_COL_SINGLETONS];
    Int rs1 = Info[UMFPACK_ROW_SINGLETONS];

    /* ---------------------------------------------------------------------- */
    /* startegy UMFPACK used*/
    /* ---------------------------------------------------------------------- */

    Int strategy = Info[UMFPACK_STRATEGY_USED];
    LUsym->strategy = strategy;

#ifndef NDEBUG
    PR = 0;
    if (strategy == UMFPACK_STRATEGY_SYMMETRIC)
    {
        PRLEVEL(PR, ("\n%% strategy used:  symmetric\n"));
        if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_AMD)
        {
            PRLEVEL(PR, ("%% ordering used:  amd on A+A'\n"));
        }
        else if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_GIVEN)
        {
            PRLEVEL(PR, ("%% ordering used: user perm.\n"));
        }
        else if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_USER)
        {
            PRLEVEL(PR, ("%% ordering used:  user function\n"));
        }
        else if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_NONE)
        {
            PRLEVEL(PR, ("%% ordering used: none\n"));
        }
        else if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_METIS)
        {
            PRLEVEL(PR, ("%% ordering used: metis on A+A'\n"));
        }
        else
        {
            PRLEVEL(PR, ("%% ordering used: not computed\n"));
        }
    }
    else
    {
        PRLEVEL(PR, ("\n%% strategy used:unsymmetric\n"));
        if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_AMD)
        {
            PRLEVEL(PR, ("%% ordering used: colamd on A\n"));
        }
        else if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_GIVEN)
        {
            PRLEVEL(PR, ("%% ordering used: user perm.\n"));
        }
        else if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_USER)
        {
            PRLEVEL(PR, ("%% ordering used: user function\n"));
        }
        else if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_NONE)
        {
            PRLEVEL(PR, ("%% ordering used: none\n"));
        }
        else if (Info[UMFPACK_ORDERING_USED] == UMFPACK_ORDERING_METIS)
        {
            PRLEVEL(PR, ("%% ordering used: metis on A'A\n"));
        }
        else
        {
            PRLEVEL(PR, ("%% ordering used: not computed\n"));
        }
    }
#endif

#ifndef NDEBUG
    PR = -1;
    PRLEVEL(PR, ("\n%% Symbolic factorization of A: "));
    if (PR <= 0) (void)umfpack_dl_report_symbolic(Symbolic, Control);
    PRLEVEL(PR, ("\n%%\tcolsingleton = %ld, rowsingleton=%ld\n", cs1, rs1));
    PR = 1;
#endif

    /* ---------------------------------------------------------------------- */
    /*    Copy the contents of Symbolic in my data structure                  */
    /* ---------------------------------------------------------------------- */
    Pinit = (Int *)paru_alloc((m + 1), sizeof(Int));
    Qinit = (Int *)paru_alloc((n + 1), sizeof(Int));
    Front_npivcol = (Int *)paru_alloc((n + 1), sizeof(Int));
    Front_1strow = (Int *)paru_alloc((n + 1), sizeof(Int));
    Front_leftmostdesc = (Int *)paru_alloc((n + 1), sizeof(Int));
    Front_parent = (Int *)paru_alloc((n + 1), sizeof(Int));
    Chain_start = (Int *)paru_alloc((n + 1), sizeof(Int));
    Chain_maxrows = (Int *)paru_alloc((n + 1), sizeof(Int));
    Chain_maxcols = (Int *)paru_alloc((n + 1), sizeof(Int));

    // temp amalgamation data structure
    Int *fmap = (Int *)paru_alloc((n + 1), sizeof(Int));
    Int *newParent = (Int *)paru_alloc((n + 1), sizeof(Int));

    if (!Pinit || !Qinit || !Front_npivcol || !Front_parent || !Chain_start ||
        !Chain_maxrows || !Chain_maxcols || !newParent || !Front_1strow ||
        !Front_leftmostdesc || !fmap)
    {
        paru_free((m + 1), sizeof(Int), Pinit);
        paru_free((n + 1), sizeof(Int), Qinit);
        paru_free((n + 1), sizeof(Int), Front_npivcol);
        paru_free((n + 1), sizeof(Int), Front_1strow);
        paru_free((n + 1), sizeof(Int), Front_leftmostdesc);
        paru_free((n + 1), sizeof(Int), Front_parent);
        paru_free((n + 1), sizeof(Int), Chain_start);
        paru_free((n + 1), sizeof(Int), Chain_maxrows);
        paru_free((n + 1), sizeof(Int), Chain_maxcols);

        paru_free((n + 1), sizeof(Int), fmap);
        paru_free((n + 1), sizeof(Int), newParent);

        printf("out of memory");

        paru_freesym(&LUsym);
        umfpack_dl_free_symbolic(&Symbolic);
        umfpack_dl_azn_free_sw(&SW);

        return NULL;
    }

    status = umfpack_dl_get_symbolic(
        &nr, &nc, &n1, &anz, &nfr, &nchains, Pinit, Qinit, Front_npivcol,
        Front_parent, Front_1strow, Front_leftmostdesc, Chain_start,
        Chain_maxrows, Chain_maxcols, Symbolic);
    if (status < 0)
    {
        printf("symbolic factorization invalid");

        // free memory
        paru_free((m + 1), sizeof(Int), Pinit);
        paru_free((n + 1), sizeof(Int), Qinit);
        paru_free((n + 1), sizeof(Int), Front_npivcol);
        paru_free((n + 1), sizeof(Int), Front_1strow);
        paru_free((n + 1), sizeof(Int), Front_leftmostdesc);
        paru_free((n + 1), sizeof(Int), Front_parent);
        paru_free((n + 1), sizeof(Int), Chain_start);
        paru_free((n + 1), sizeof(Int), Chain_maxrows);
        paru_free((n + 1), sizeof(Int), Chain_maxcols);

        paru_freesym(&LUsym);
        umfpack_dl_free_symbolic(&Symbolic);
        umfpack_dl_azn_free_sw(&SW);

        return NULL;
    }
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%%%% n1 is %ld\n", n1));
    PRLEVEL(PR, ("From the Symbolic object,\
                C is of dimension %ld-by-%ld\n",
                 nr, nc));
    PRLEVEL(PR, ("   with nz = %ld, number of fronts = %ld,\n", anz, nfr));
    PR = 1;
    PRLEVEL(PR, ("   number of frontal matrix chains = %ld\n", nchains));

    PRLEVEL(1, ("\nPivot columns in each front, and parent of each front:\n"));
    Int k = 0;

    PR = 1;
    for (Int i = 0; i < nfr; i++)
    {
        Int fnpiv = Front_npivcol[i];
        PRLEVEL(PR, ("Front %ld: parent front: %ld number of pivot cols: %ld\n",
                     i, Front_parent[i], fnpiv));
        PRLEVEL(PR, ("%% first row is %ld\n", Front_1strow[i]));

        for (Int j = 0; j < fnpiv; j++)
        {
            Int col = Qinit[k];
            PRLEVEL(PR, ("%ld-th pivot column is column %ld"
                         " in original matrix\n",
                         k, col));
            k++;
        }
    }
    PR = 1;

    PRLEVEL(PR, ("\nTotal number of pivot columns "
                 "in frontal matrices: %ld\n",
                 k));

    PRLEVEL(PR, ("\nFrontal matrix chains:\n"));
    for (Int j = 0; j < nchains; j++)
    {
        PRLEVEL(PR, ("Frontal matrices %ld to %ld in chain\n", Chain_start[j],
                     Chain_start[j + 1] - 1));
        PRLEVEL(PR, ("\tworking array of size %ld-by-%ld\n", Chain_maxrows[j],
                     Chain_maxcols[j]));
    }
#endif

    umfpack_dl_free_symbolic(&Symbolic);

    paru_free((n + 1), sizeof(Int), Front_leftmostdesc);
    paru_free((n + 1), sizeof(Int), Front_1strow);

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    ASSERT(m == nr);
    ASSERT(n == nc);

    LUsym->m = m;
    LUsym->n = n;
    LUsym->n1 = n1;
    LUsym->rs1 = rs1;
    LUsym->cs1 = cs1;
    LUsym->anz = anz;
    Int nf = LUsym->nf = nfr;
    LUsym->Chain_start = Chain_start;
    LUsym->Chain_maxrows = Chain_maxrows;
    LUsym->Chain_maxcols = Chain_maxcols;
    LUsym->Qfill = Qinit;

    PRLEVEL(0, ("%% A  is  %ld x %ld \n", m, n));
    PRLEVEL(-1, ("LU = zeros(%ld,%ld);\n", m, n));
    PRLEVEL(1, ("npivots =[]; \n"));
    PRLEVEL(1, ("S = zeros(%ld,%ld); %% n1 = %ld\n", m, n, n1));
    PRLEVEL(1, ("%% nf=%ld\n", nf));
    //
    /* ---------------------------------------------------------------------- */
    /*           Fixing Parent and computing Children datat structure         */
    /* ---------------------------------------------------------------------- */

    // Parent size is nf+1 potentially smaller than what UMFPACK allocate
    size_t size = n + 1;
    Int *Parent = (Int *)paru_realloc(nf + 1, sizeof(Int), Front_parent, &size);
    ASSERT(size <= (size_t)n + 1);
    // TODO: nf == 0 is a weird condintion I have to check
    if (Parent == NULL || nf == 0)
    {  // should not happen anyway it is always shrinking
        printf("memory problem");
        // free memory
        paru_free((n + 1), sizeof(Int), Front_npivcol);
        paru_free((n + 1), sizeof(Int), Front_parent);
        paru_free((m + 1), sizeof(Int), Pinit);
        paru_freesym(&LUsym);
        umfpack_dl_azn_free_sw(&SW);
        return NULL;
    }
    LUsym->Parent = Parent;

    // Making Super data structure
    // like SPQR: Super[f]<= pivotal columns of (f) < Super[f+1]
    Int *Super = LUsym->Super = (Int *)paru_alloc((nf + 1), sizeof(Int));
    if (Super == NULL)
    {
        printf("memory problem");
        paru_free((m + 1), sizeof(Int), Pinit);
        paru_freesym(&LUsym);
        umfpack_dl_azn_free_sw(&SW);
        return NULL;
    }
    Super[0] = 0;
    for (Int k = 1; k <= nf; k++)
    {
        Super[k] = Front_npivcol[k - 1];
    }
    paru_cumsum(nf + 1, Super);

    /* ---------------------------------------------------------------------- */
    /*                          Relaxed amalgamation                          */
    /* ---------------------------------------------------------------------- */
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%%%% Before relaxed amalgmation\n"));

    PRLEVEL(PR, ("%%%% Super:\n"));
    for (Int k = 0; k <= nf; k++) PRLEVEL(PR, ("  %ld", Super[k]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("%%%% Parent:\n"));
    for (Int k = 0; k <= nf; k++) PRLEVEL(PR, ("  %ld", Parent[k]));
    PRLEVEL(PR, ("\n"));
    PR = 1;
#endif

    // Relax amalgamation before making the list of children
    // The gist is to make a cumsum over Pivotal columns
    // then start from children and see if I merge it to the father how many
    // pivotal columns will it have;
    // if it has less than threshold merge the child to the father
    //
    // Super and Parent  and upperbounds
    // Parent needs an extra work space
    // Super can be changed in space upperbound can be changed in space
    // Upperbound how to do: maximum of pervious upperbounds
    // Number of the columns of the root of each subtree
    //
    Int threshold = 32;
    Int newF = 0;

    for (Int f = 0; f < nf; f++)
    {  // finding representative for each front
        Int repr = f;
        // amalgamate till number of pivot columns is small
        PRLEVEL(PR, ("%% repr = %ld Parent =%ld\n", repr, Parent[repr]));
        PRLEVEL(PR, ("%%size of Potential pivot= %ld\n",
                     Super[Parent[repr] + 1] - Super[f]));
        while (Super[Parent[repr] + 1] - Super[f] < threshold &&
               Parent[repr] != -1)
        {
            repr = Parent[repr];
            PRLEVEL(PR, ("%%Middle stage f= %ld repr = %ld\n", f, repr));
            PRLEVEL(PR, ("%%number of pivot cols= %ld\n",
                         Super[repr + 1] - Super[f]));
            PRLEVEL(PR, ("%%number of pivot cols if Parent collapsed= %ld\n",
                         Super[Parent[repr] + 1] - Super[f]));
        }

        PRLEVEL(PR, ("%% newF = %ld for:\n", newF));
        for (Int k = f; k <= repr; k++)
        {
            PRLEVEL(PR, ("%%  %ld ", k));
            fmap[k] = newF;
        }
        PRLEVEL(PR, ("%%repr = %ld\n", repr));
        newF++;
        f = repr;
    }

#ifndef NDEBUG
    PR = 1;
#endif

    Int newNf = newF;  // new size of number of fronts
    fmap[nf] = -1;
    // nf =  LUsym->nf = newF;
    // newParent size is newF+1 potentially smaller than nf
    newParent = (Int *)paru_realloc(newF + 1, sizeof(Int), newParent, &size);
    ASSERT(newF <= nf);
    // TODO: add memory guard?
    // Int newSuper[newNf+2];

    for (Int oldf = 0; oldf < nf; oldf++)
    {  // maping old to new
        Int newf = fmap[oldf];
        Int oldParent = Parent[oldf];
        newParent[newf] = oldParent >= 0 ? fmap[oldParent] : -1;
        // newSuper[newf] = Super[oldf] ;
        // Super[newf] = Super[oldf] ;
    }
    PRLEVEL(-1, ("%% newF = %ld and nf=%ld\n", newNf, nf));
    // newSuper[newNf] = Super[nf] ;

    /* ---------------------------------------------------------------------- */
    /*         Finding the Upper bound of rows and cols                       */
    /* ---------------------------------------------------------------------- */
    SWType *mySW = (SWType *)SW;
    Int *Front_nrows = (Int *)mySW->Front_nrows;
    Int *Front_ncols = (Int *)mySW->Front_ncols;

    LUsym->Fm = NULL;  // Upper bound on number of rows including pivots
    LUsym->Cm = NULL;  // Upper bound on number of columns excluding pivots
    //    LUsym->Fm = QRsym->Fm;
    //    QRsym->Fm = NULL;
    Int *Fm = (Int *)paru_calloc((newNf + 1), sizeof(Int));
    Int *Cm = (Int *)paru_alloc((newNf + 1), sizeof(Int));
    LUsym->Fm = Fm;
    LUsym->Cm = Cm;

    // TODO: I have not checked memory problems after changin the code
    if (Fm == NULL || Cm == NULL)
    {
        printf("memory problem");
        paru_freesym(&LUsym);
        umfpack_dl_azn_free_sw(&SW);
        return NULL;
    }

    // after relaxed amalgamation
    // Copying first nf+1 elements of Front_nrows(UMFPACK) into Fm(SPQR like)
    for (Int oldf = 0; oldf < nf; oldf++)
    {
        PRLEVEL(PR, ("oldf=%ld\n", oldf));
        Int newf = fmap[oldf];
        PRLEVEL(PR, ("newf=%ld\n", newf));
        PRLEVEL(PR, ("next=%ld\n", fmap[oldf + 1]));

        if (newf != fmap[oldf + 1])
        {                                   // either root or not amalgamated
            Fm[newf] += Front_nrows[oldf];  // + Front_npivcol[oldf];
            Cm[newf] = Front_ncols[oldf] - Front_npivcol[oldf];
            // newSuper[newf+1] = Super[oldf+1] ;
            Super[newf + 1] = Super[oldf + 1];
            PRLEVEL(PR, ("Fm[newf]=%ld\n", Fm[newf]));
            PRLEVEL(PR, ("Cm[newf]=%ld\n", Cm[newf]));
        }
        else
        {
            Fm[newf] += Front_npivcol[oldf];
            // Cm[newf] += Front_npivcol[oldf];
            PRLEVEL(PR, ("Fm[newf]=%ld\n", Fm[newf]));
            PRLEVEL(PR, ("Cm[newf]=%ld\n", Cm[newf]));
        }
    }
    Super[newNf] = Super[nf];

    LUsym->Rj = NULL;
    //    LUsym->Rj = QRsym->Rj;
    //    QRsym->Rj = NULL;

    LUsym->Rp = NULL;
    //    LUsym->Rp = QRsym->Rp;
    //    QRsym->Rp = NULL;

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%%%% After relaxed amalgmation\n"));

    PRLEVEL(PR, ("Cm =\n"));
    for (Int i = 0; i < nf + 1; i++) PRLEVEL(PR, ("%ld ", Cm[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Fm =\n"));
    for (Int i = 0; i < nf + 1; i++) PRLEVEL(PR, ("%ld ", Fm[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Pivot cols=\n"));
    for (Int i = 0; i < nf + 1; i++) PRLEVEL(PR, ("%ld ", Front_npivcol[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Upper bound on Rows =\n"));
    for (Int i = 0; i < nf + 1; i++) PRLEVEL(PR, ("%ld ", Front_nrows[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Upper bound on Cols=\n"));
    for (Int i = 0; i < nf + 1; i++) PRLEVEL(PR, ("%ld ", Front_ncols[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("%%%% Super:\n"));
    for (Int k = 0; k <= nf; k++) PRLEVEL(PR, ("  %ld", Super[k]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("%%%% fmap:\n"));
    for (Int k = 0; k <= nf; k++) PRLEVEL(PR, ("  %ld", fmap[k]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("%%%% newParent:\n"));
    for (Int k = 0; k <= newNf; k++) PRLEVEL(PR, ("  %ld", newParent[k]));
    PRLEVEL(PR, ("\n"));

#endif

    paru_free(nf + 1, sizeof(Int), LUsym->Parent);
    LUsym->Parent = Parent = newParent;  // TODO:free stuff not necessary
    nf = LUsym->nf = newNf;

    umfpack_dl_azn_free_sw(&SW);
    paru_free((n + 1), sizeof(Int), Front_npivcol);
    paru_free((n + 1), sizeof(Int), fmap);

    //////////////////////end of relaxed amalgamation/////////////////////////

    // Making Children list and computing the bound sizes
    Int *Childp = (Int *)paru_calloc((nf + 2), sizeof(Int));
    LUsym->Childp = Childp;
    if (Childp == NULL)
    {
        printf("memory problem");
        paru_free((m + 1), sizeof(Int), Pinit);
        paru_freesym(&LUsym);
        // umfpack_dl_azn_free_sw (&SW);
        return NULL;
    }

#ifndef NDEBUG
    Int Us_bound_size = 0;
    Int LUs_bound_size = 0;
    Int row_Int_bound = 0;
    Int col_Int_bound = 0;
#endif
    for (Int f = 0; f < nf; f++)
    {
        // TODO
        // work[f]= fp*fm*fn;
        // sum(work[first[f]] ...work[f]) task size
#ifndef NDEBUG
        Int fp = Super[f + 1] - Super[f];
        Int fm = LUsym->Fm[f];
        Int fn = LUsym->Cm[f]; /* Upper bound number of cols of F */
        Us_bound_size += fp * fn;
        LUs_bound_size += fp * fm;
        row_Int_bound += fm;
        col_Int_bound += fn;
#endif
        if (Parent[f] > 0) Childp[Parent[f] + 1]++;
    }
    // see GraphBLAS/Source/GB_cumsum.c
    paru_cumsum(nf + 2, Childp);
#ifndef NDEBUG
    LUsym->Us_bound_size = Us_bound_size;
    LUsym->LUs_bound_size = LUs_bound_size;
    LUsym->row_Int_bound = row_Int_bound;
    LUsym->col_Int_bound = col_Int_bound;
    PR = 1;
    PRLEVEL(PR, ("%%row_Int_bound=%ld, col_Int_bound=%ld", row_Int_bound,
                 col_Int_bound));
    PRLEVEL(PR,
            ("%%-Us_bound_size = %ld LUs_bound_size = %ld sum = %ld\n",
             Us_bound_size, LUs_bound_size,
             row_Int_bound + col_Int_bound + Us_bound_size + LUs_bound_size));
    PR = 1;
    PRLEVEL(PR, ("%%%%-Chidlp-----"));
    for (Int f = 0; f < nf + 2; f++) PRLEVEL(PR, ("%ld ", Childp[f]));
    PRLEVEL(PR, ("\n"));
#endif
    Int *Child = (Int *)paru_calloc((nf + 1), sizeof(Int));
    LUsym->Child = Child;
    if (Child == NULL)
    {
        printf("memory problem");
        paru_free((m + 1), sizeof(Int), Pinit);
        paru_freesym(&LUsym);
        // umfpack_dl_azn_free_sw (&SW);
        return NULL;
    }

    // copy of Childp using Work for other places also
    Int *Work = (Int *)paru_alloc((MAX(m, n) + 2), sizeof(Int));
    Int *cChildp = Work;
    if (cChildp == NULL)
    {
        printf("memory problem");
        paru_free((m + 1), sizeof(Int), Pinit);
        paru_free((MAX(m, n) + 2), sizeof(Int), Work);
        paru_freesym(&LUsym);
        // umfpack_dl_azn_free_sw (&SW);
        return NULL;
    }

    // see GraphBLAS/Source/GB_memcpy
    paru_memcpy(cChildp, Childp, (nf + 2) * sizeof(Int));

    for (Int f = 0; f < nf; f++)
    {
        if (Parent[f] > 0) Child[cChildp[Parent[f]]++] = f;
    }
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%%%%_cChidlp_____"));
    for (Int f = 0; f < nf + 2; f++) PRLEVEL(PR, ("%ld ", cChildp[f]));
    PRLEVEL(PR, ("\n"));
#endif

    /* ---------------------------------------------------------------------- */
    /*                   computing the Staircase structures                   */
    /* ---------------------------------------------------------------------- */

    Int *Sp = LUsym->Sp = (Int *)paru_calloc(m + 1 - n1, sizeof(Int));
    Int *Sleft = LUsym->Sleft = (Int *)paru_alloc(n + 2 - n1, sizeof(Int));
    Int *Pinv = (Int *)paru_alloc(m + 1, sizeof(Int));

    if (Sp == NULL || Sleft == NULL || Pinv == NULL)
    {
        printf("memory problem");

        paru_free((m + 1), sizeof(Int), Pinit);
        paru_free((MAX(m, n) + 2), sizeof(Int), Work);

        paru_freesym(&LUsym);
        // umfpack_dl_azn_free_sw (&SW);

        paru_free(m, sizeof(Int), Pinv);
        return NULL;
    }

    //-------- computing the inverse permutation for P
    for (Int i = 0; i < m; i++)
    {
        Pinv[Pinit[i]] = i;
    }

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("Qinit =\n"));
    for (Int j = 0; j < m; j++) PRLEVEL(PR, ("%ld ", Qinit[j]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Pinit =\n"));
    for (Int i = 0; i < m; i++) PRLEVEL(PR, ("%ld ", Pinit[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Pinv =\n"));
    for (Int i = 0; i < m; i++) PRLEVEL(PR, ("%ld ", Pinv[i]));
    PRLEVEL(PR, ("\n"));

#endif

    Int *Ps;  // new row permutation for just the Submatrix part

    Int *Sup = NULL;   // Singlton u p
    Int *cSup = NULL;  // copy of Singlton u p
    Int *Slp = NULL;   // Singlton l p
    Int *cSlp = NULL;  // copy Singlton l p
    Ps = (Int *)paru_calloc(m - n1, sizeof(Int));
    if (cs1 != 0)
    {
        Sup = LUsym->ustons.Sup = (Int *)paru_calloc(cs1 + 1, sizeof(Int));
        cSup = (Int *)paru_calloc(cs1 + 1, sizeof(Int));
    }
    if (rs1 != 0)
    {
        Slp = LUsym->lstons.Slp = (Int *)paru_calloc(rs1 + 1, sizeof(Int));
        cSlp = (Int *)paru_calloc(rs1 + 1, sizeof(Int));
    }

    if (((Slp == NULL || cSlp == NULL) && rs1 != 0) ||
        ((Sup == NULL || cSup == NULL) && cs1 != 0) || Ps == NULL)
    {
        printf("rs1=%ld cs1=%ld memory problem\n", rs1, cs1);
        paru_free((cs1 + 1), sizeof(Int), Sup);
        paru_free((cs1 + 1), sizeof(Int), cSup);

        paru_free((rs1 + 1), sizeof(Int), Slp);
        paru_free((rs1 + 1), sizeof(Int), cSlp);

        paru_free((m + 1), sizeof(Int), Pinit);
        paru_free((MAX(m, n) + 2), sizeof(Int), Work);
        paru_free(m, sizeof(Int), Pinv);
        // umfpack_dl_azn_free_sw (&SW);
        paru_freesym(&LUsym);
        return NULL;
    }
    Int sunz = 0;  // U nnz: singlteton nnzero of s
    Int slnz = 0;  // L nnz: singlteton nnzero of s
    Int snz = 0;   // s nonzero: nnz in submatrix excluding singletons
    Int rowcount = 0;
    Sleft[0] = 0;
    // counting number of entries in each row of submatrix Sp and also making
    // the singelton matrices
#ifndef NDEBUG
    PR = -1;
    PRLEVEL(PR, ("Computing Staircase Structure and singleton structure\n"));
    PRLEVEL(PR, ("rs1= %ld cs1=%ld\n", rs1, cs1));
    PR = 1;
#endif
    for (Int newcol = 0; newcol < n1; newcol++)
    {  // The columns that are just in singleton
        Int oldcol = Qinit[newcol];
        PRLEVEL(PR, ("newcol = %ld oldcol=%ld\n", newcol, oldcol));
        for (Int p = Ap[oldcol]; p < Ap[oldcol + 1]; p++)
        {
            Int oldrow = Ai[p];
            Int newrow = Pinv[oldrow];
            PRLEVEL(PR, ("newrow=%ld oldrow=%ld\n", newrow, oldrow));
            if (newrow < cs1)
            {  // inside U singletons
                PRLEVEL(PR, ("Inside U singletons\n"));
                sunz++;
                Sup[newrow + 1]++;
            }
            else
            {  // inside L singletons
#ifndef NDEBUG
                PR = 1;
#endif
                PRLEVEL(PR, ("Inside L singletons\n"));
                slnz++;
                Slp[newcol - cs1 + 1]++;
            }
        }
    }
    LUsym->ustons.nnz = sunz;
    LUsym->lstons.nnz = slnz;
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("Sup and Slp in the middle\n"));
    if (cs1 > 0)
    {
        PRLEVEL(PR, ("(%ld) Sup =", sunz));
        for (Int k = 0; k <= cs1; k++) PRLEVEL(PR, ("%ld ", Sup[k]));
        PRLEVEL(PR, ("\n"));
    }
    if (rs1 > 0)
    {
        PRLEVEL(PR, ("(%ld) Slp =", slnz));
        for (Int k = cs1; k <= n1; k++) PRLEVEL(PR, ("%ld ", Slp[k - cs1]));
        PRLEVEL(PR, ("\n"));
    }
#endif
    for (Int newcol = n1; newcol < n; newcol++)
    {
        Int oldcol = Qinit[newcol];
        PRLEVEL(1, ("newcol = %ld oldcol=%ld\n", newcol, oldcol));
        for (Int p = Ap[oldcol]; p < Ap[oldcol + 1]; p++)
        {
            Int oldrow = Ai[p];
            Int newrow = Pinv[oldrow];
            Int srow = newrow - n1;
            PRLEVEL(1, ("\tnewrow=%ld oldrow=%ld srow=%ld\n", newrow, oldrow,
                        srow));
            if (srow >= 0)
            {  // it is insdie S otherwise it is part of singleton
                if (Sp[srow + 1] == 0)
                {  // first time seen
                    PRLEVEL(1, ("\tPs[%ld]= %ld\n", rowcount, srow));
                    Ps[rowcount] = srow;
                    Pinit[n1 + rowcount] = oldrow;
                    rowcount++;
                }
                snz++;
                Sp[srow + 1]++;
            }
            else
            {  // inside the U singletons
                Sup[newrow + 1]++;
                sunz++;
            }
        }
        Sleft[newcol - n1 + 1] = rowcount;
    }
    Sleft[n - n1 + 1] = m - n1 - rowcount;  // empty rows of S if any
    LUsym->snz = snz;
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("Sup and Slp finished (before cumsum)U-sing =%ld L-sing=%ld\n", 
                sunz, slnz));
    if (cs1 > 0)
    {
        PRLEVEL(PR, ("(%ld) Sup =", sunz));
        for (Int k = 0; k <= cs1; k++) PRLEVEL(PR, ("%ld ", Sup[k]));
        PRLEVEL(PR, ("\n"));
    }
    if (rs1 > 0)
    {
        PRLEVEL(PR, ("(%ld) Slp =", slnz));
        for (Int k = cs1; k <= n1; k++) PRLEVEL(PR, ("%ld ", Slp[k - cs1]));
        PRLEVEL(PR, ("\n"));
    }
    PR = 1;
    PRLEVEL(PR, ("Ps =\n"));
    for (Int k = 0; k < rowcount; k++) PRLEVEL(PR, ("%ld ", Ps[k]));
    PRLEVEL(PR, ("\n"));
#endif

    if (rowcount < m - n1)
    {
        // I think that must not happen anyway while umfpack finds it
        printf("Empty rows in submatrix\n");

#ifndef NDEBUG
        PRLEVEL(1, ("m = %ld, n1 = %ld, rowcount = %ld, snz = %ld\n", m, n1,
                    rowcount, snz));
        for (Int srow = 0; srow < m - n1; srow++)
        {
            if (Sp[srow] == 0)
            {
                PRLEVEL(1, ("Row %ld is empty\n", srow));
            }
        }
#endif
        paru_free((m - n1), sizeof(Int), Ps);
        paru_free((m + 1), sizeof(Int), Pinit);
        paru_free((MAX(m, n) + 2), sizeof(Int), Work);

        paru_free(m, sizeof(Int), Pinv);
        paru_freesym(&LUsym);
        // umfpack_dl_azn_free_sw (&SW);
        return NULL;  // Free memory
    }
    ASSERT(rowcount == m - n1);

    // Update Sp based on new permutation for Sp[1...m+1]
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("Before permutation Sp =\n"));
    for (Int i = n1; i <= m; i++) PRLEVEL(PR, ("%ld ", Sp[i - n1]));
    PRLEVEL(PR, ("\n"));
#endif

    // update Pinv
    for (Int i = n1; i < m; i++) Pinv[Pinit[i]] = i;

    ///////////////////////////////////////////////////////////////
    Int *cSp = Work;
    for (Int i = n1; i < m; i++)
    {
        Int row = i - n1;
        PRLEVEL(1, ("Permutation row = %ld, Ps[row]= %ld, Sp[row]=%ld\n", row,
                    Ps[row], Sp[row]));
        cSp[row + 1] = Sp[Ps[row] + 1];
    }

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("After permutation Sp =\n"));
    for (Int i = n1; i <= m; i++) PRLEVEL(PR, ("%ld ", cSp[i - n1]));
    PRLEVEL(PR, ("\n"));
#endif
    // TODO is it too much memcpy?
    // Prior falsy version:
    // paru_memcpy(Sp, cSp, (m + 1 - n1) * sizeof(Int));
    // paru_cumsum(m + 1 - n1, Sp);
    // paru_memcpy(cSp, Sp, (m + 1 - n1) * sizeof(Int));

    paru_cumsum(m + 1 - n1, cSp);
    paru_memcpy(Sp, cSp, (m + 1 - n1) * sizeof(Int));

    if (cs1 > 0)
    {
        paru_cumsum(cs1+1, Sup);
        paru_memcpy(cSup, Sup, cs1);
    }
    if (rs1 > 0)
    {
        paru_cumsum(rs1+1, Slp);
        paru_memcpy(cSlp, Slp, rs1);
    }

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("After cumsum  Sp =\n"));
    for (Int i = n1; i <= m; i++) PRLEVEL(PR, ("%ld ", Sp[i - n1]));
    PRLEVEL(PR, ("\n"));
#endif

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("After Stair case Pinv =\n"));
    for (Int i = 0; i < m; i++) PRLEVEL(PR, ("%ld ", Pinv[i]));
    PRLEVEL(PR, ("\n"));
#endif

    // PofA
    LUsym->Pinit = Pinit;
    ///////////////////////////////////////////////////////////////
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("Before Sj Sp =\n"));
    for (Int i = n1; i <= m; i++) PRLEVEL(PR, ("%ld ", cSp[i - n1]));
    PRLEVEL(PR, ("\n"));
#endif

    Int *Suj = NULL;
    double *Sux = NULL;
    if (cs1 > 0)
    {
        Suj = (Int *)paru_alloc(sunz, sizeof(Int));
        Sux = (double *)paru_alloc(sunz, sizeof(double));
    }
    LUsym->ustons.Sux = Sux;
    LUsym->ustons.Suj = Suj;

    Int *Sli = NULL;
    double *Slx = NULL;
    if (rs1 > 0)
    {
        Sli = (Int *)paru_alloc(slnz, sizeof(Int));
        Slx = (double *)paru_alloc(slnz, sizeof(double));
    }
    LUsym->lstons.Slx = Slx;
    LUsym->lstons.Sli = Sli;

    // Updating Sj and Sx using copy of Sp
    Int *Sj = (Int *)paru_alloc(snz, sizeof(Int));
    double *Sx = (double *)paru_alloc(snz, sizeof(double));

    LUsym->Sj = Sj;
    LUsym->Sx = Sx;

    if (Sj == NULL || Sx == NULL || (cs1 > 0 && (Suj == NULL || Sux == NULL)) ||
        (rs1 > 0 && (Sli == NULL || Slx == NULL)))
    {
        printf("memory problem");
        paru_free(m, sizeof(Int), Pinv);
        paru_freesym(&LUsym);
        // umfpack_dl_azn_free_sw (&SW);
        return NULL;
    }

    // construct Sj, Sx and singltons
#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("Constructing Sj and singletons\n"));
#endif
    for (Int newcol = 0; newcol < n1; newcol++)
    {  // The columns that are just in singleton
        Int oldcol = Qinit[newcol];
        PRLEVEL(PR, ("newcol = %ld oldcol=%ld\n", newcol, oldcol));
        for (Int p = Ap[oldcol]; p < Ap[oldcol + 1]; p++)
        {
            Int oldrow = Ai[p];
            Int newrow = Pinv[oldrow];
            PRLEVEL(PR, ("newrow=%ld oldrow=%ld\n", newrow, oldrow));
            if (newrow < cs1)
            {  // inside U singletons CSR
                PRLEVEL(PR, ("Inside U singletons\n"));
                if (newcol == newrow)
                {  // diagonal entry
                    Suj[Sup[newrow]] = newcol;
                    Sux[Sup[newrow]] = Ax[p];
                }
                else
                {
                    Suj[++cSup[newrow]] = newcol;
                    Sux[cSup[newrow]] = Ax[p];
                }
            }
            else
            {  // inside L singletons CSC
                PRLEVEL(PR, ("Inside L singletons\n"));
                if (newcol == newrow)
                {  // diagonal entry
                    Sli[Slp[newcol - cs1]] = newrow;
                    Slx[Slp[newcol - cs1]] = Ax[p];
                }
                else
                {
                    Sli[++cSlp[newcol - cs1]] = newrow;
                    Slx[cSlp[newcol - cs1]] = Ax[p];
                }
            }
        }
    }

    PRLEVEL(PR, ("The rest of matrix after singletons\n"));
    for (Int newcol = n1; newcol < n; newcol++)
    {
        Int oldcol = Qinit[newcol];
        for (Int p = Ap[oldcol]; p < Ap[oldcol + 1]; p++)
        {
            Int oldrow = Ai[p];
            Int newrow = Pinv[oldrow];
            Int srow = newrow - n1;
            Int scol = newcol - n1;
            if (srow >= 0)
            {  // it is insdie S otherwise it is part of singleton
                Sj[cSp[srow]] = scol;
                Sx[cSp[srow]++] = Ax[p];
            }
            else
            {  // inside the U singletons
                PRLEVEL(PR, ("Usingleton rest newcol = %ld newrow=%ld\n",
                             newcol, newrow));
                ASSERT(newrow != newcol);  // not a diagonal entry
                Suj[++cSup[newrow]] = newcol;
                Sux[cSup[newrow]] = Ax[p];
            }
        }
    }
    if (cs1 > 0) paru_free((cs1 + 1), sizeof(Int), cSup);
    if (rs1 > 0) paru_free((rs1 + 1), sizeof(Int), cSlp);

    paru_free((m - n1), sizeof(Int), Ps);
    PRLEVEL(PR, ("Constructing Sj and singletons finished here\n"));
#ifndef NDEBUG
    PR = -1;
    PRLEVEL(PR, ("Sup and Slp after mading Sux Slx\n"));
    if (cs1 > 0)
    {
        PRLEVEL(PR, ("(%ld) Sup =", sunz));
        for (Int k = 0; k <= cs1; k++) PRLEVEL(PR, ("%ld ", Sup[k]));
        PRLEVEL(PR, ("\n"));

        for(Int newrow = 0; newrow < cs1; newrow ++)
        {
            PRLEVEL(PR, ("row = %ld\n", newrow));
            for(Int p = Sup[newrow]; p < Sup[newrow+1]; p++)
            {
                PRLEVEL(PR, (" (%ld)%.2lf", Suj[p],Sux[p]));
            }
            PRLEVEL(PR, ("\n"));
        }
    }
    if (rs1 > 0)
    {
        PRLEVEL(PR, ("(%ld) Slp =", slnz));
        for (Int k = cs1; k <= n1; k++) PRLEVEL(PR, ("%ld ", Slp[k - cs1]));
        PRLEVEL(PR, ("\n"));

        for(Int newcol = cs1; newcol < n1; newcol++)
        {
            PRLEVEL(PR, ("col = %ld\n", newcol+cs1));
            for(Int p = Slp[newcol-cs1]; p < Slp[newcol-cs1+1]; p++)
            {
                PRLEVEL(PR, (" (%ld)%.2lf", Sli[p],Slx[p]));
            }
            PRLEVEL(PR, ("\n"));
        }
    }
    PR = 1;
    PRLEVEL(PR, ("Sp =\n"));
    for (Int i = n1; i <= m; i++) PRLEVEL(PR, ("%ld ", Sp[i - n1]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("Sj =\n"));
    for (Int k = 0; k < snz; k++) PRLEVEL(PR, ("%ld ", Sj[k]));
    PRLEVEL(PR, ("\n"));
    PR = 1;
    for (Int i = 0; i < m - n1; i++)
    {
        PRLEVEL(PR, (" i=%ld cSp=%ld Sp=%ld\n", i, cSp[i], Sp[i + 1]));
        ASSERT(cSp[i] == Sp[i + 1]);
    }
#endif
    paru_free((MAX(m, n) + 2), sizeof(Int), Work);

    /* ---------------------------------------------------------------------- */
    /*    computing the augmented tree and the data structures related        */
    /* ---------------------------------------------------------------------- */

    /*Computing augmented tree */
    Int *aParent = LUsym->aParent = NULL;  // augmented tree size m+nf
    Int *aChildp = LUsym->aChildp = NULL;  // size m+nf+2
    Int *aChild = LUsym->aChild = NULL;    // size m+nf+1
    Int *rM = LUsym->row2atree = NULL;     // row map
    Int *snM = LUsym->super2atree = NULL;  // and supernode map
    Int *first = LUsym->first = NULL;      // first descendent in the tree
    // augmented tree size nf+1

#ifndef NDEBUG
    PR = 1;
    // HERE
    /* print fronts*/
    for (Int f = 0; f < nf; f++)
    {
        Int col1 = Super[f];
        Int col2 = Super[f + 1];
        Int fp = Super[f + 1] - Super[f];
        PRLEVEL(PR, ("%% Front=%ld Parent=%ld\n", f, Parent[f]));
        PRLEVEL(PR, ("%% %ld...%ld npivotal=%ld\n", col1, col2, fp));
        PRLEVEL(PR, ("%% list of %ld children: ", Childp[f + 1] - Childp[f]));
        for (Int i = Childp[f]; i <= Childp[f + 1] - 1; i++)
            PRLEVEL(PR, ("%ld ", Child[i]));
        PRLEVEL(PR, ("\n"));
    }
    PR = 1;
#endif
    // submatrix is msxns
    Int ms = m - n1;
#ifndef NDEBUG
    Int ns = n - n1;
#endif

    LUsym->aParent = aParent = (Int *)paru_alloc(ms + nf, sizeof(Int));
    LUsym->aChild = aChild = (Int *)paru_alloc(ms + nf + 1, sizeof(Int));
    LUsym->aChildp = aChildp = (Int *)paru_alloc(ms + nf + 2, sizeof(Int));
    LUsym->first = first = (Int *)paru_alloc(nf + 1, sizeof(Int));
    LUsym->row2atree = rM = (Int *)paru_alloc(ms, sizeof(Int));
    LUsym->super2atree = snM = (Int *)paru_alloc(nf, sizeof(Int));

    double *front_flop_bound = NULL;
    double *stree_flop_bound = NULL;
    LUsym->front_flop_bound = front_flop_bound =
        (double *)paru_alloc(nf + 1, sizeof(double));
    LUsym->stree_flop_bound = stree_flop_bound =
        (double *)paru_calloc(nf + 1, sizeof(double));

    if (aParent == NULL || aChild == NULL || aChildp == NULL || rM == NULL ||
        snM == NULL || first == NULL || front_flop_bound == NULL ||
        stree_flop_bound == NULL)
    {
        printf("Out of memory in symbolic phase");
        paru_free(m, sizeof(Int), Pinv);
        paru_freesym(&LUsym);
        return NULL;
    }
    // initialization
    paru_memset(aParent, -1, (ms + nf) * sizeof(Int));
#ifndef NDEBUG  // TODO: should it be in debug mode?
    paru_memset(aChild, -1, (ms + nf + 1) * sizeof(Int));
    PR = 1;
#endif
    paru_memset(aChildp, -1, (ms + nf + 2) * sizeof(Int));
    paru_memset(first, -1, (nf + 1) * sizeof(Int));

    aChildp[0] = 0;
    Int offset = 0;  // number of rows visited in each iteration orig front+
    // rows
    Int lastChildFlag = 0;
    Int childpointer = 0;

    for (Int f = 0; f < nf; f++)
    {
        PRLEVEL(PR, ("%% Front %ld\n", f));
        PRLEVEL(PR, ("%% pivot columns [ %ld to %ld ] n: %ld \n", Super[f],
                     Super[f + 1] - 1, ns));

        // computing works in each front
        Int fp = Super[f + 1] - Super[f];  // k
        Int fm = LUsym->Fm[f];             // m
        Int fn = LUsym->Cm[f];             // n Upper bound number of cols of f
        front_flop_bound[f] = (double)(fp * fm * fn + fp * fm + fp * fn);
        stree_flop_bound[f] += front_flop_bound[f];
        for (Int i = Childp[f]; i <= Childp[f + 1] - 1; i++)
        {
            stree_flop_bound[f] += stree_flop_bound[Child[i]];
            PRLEVEL(PR, ("%% child=%ld fl=%lf ", Child[i],
                         front_flop_bound[Child[i]]));
        }
        PRLEVEL(PR, ("%% flops bound= %lf\n ", front_flop_bound[f]));
        PRLEVEL(PR, ("%% %ld %ld %ld\n ", fp, fm, fn));
        PRLEVEL(PR, ("%% stree bound= %lf\n ", stree_flop_bound[f]));

        ASSERT(Super[f + 1] <= ns);
        Int numRow = Sleft[Super[f + 1]] - Sleft[Super[f]];

        Int numoforiginalChild = 0;
        if (lastChildFlag)
        {  // the current node is the parent
            PRLEVEL(PR, ("%% Childs of %ld: ", f));
            numoforiginalChild = Childp[f + 1] - Childp[f];

            for (Int i = Childp[f]; i <= Childp[f + 1] - 1; i++)
            {
                PRLEVEL(PR, ("%ld ", Child[i]));
                Int c = Child[i];
                ASSERT(snM[c] < ms + nf + 1);
                aParent[snM[c]] = offset + numRow;
                PRLEVEL(PR, ("%% aParent[%ld] =%ld\n", aParent[snM[c]],
                             offset + numRow));
                ASSERT(childpointer < ms + nf + 1);
                aChild[childpointer++] = snM[c];
            }
        }

        PRLEVEL(1, ("%% numRow=%ld ", numRow));
        PRLEVEL(1, ("#offset=%ld\n", offset));
        for (Int i = offset; i < offset + numRow; i++)
        {
            ASSERT(aChildp[i + 1] == -1);
            aChildp[i + 1] = aChildp[i];
            PRLEVEL(1, ("%% @i=%ld aCp[%ld]=%ld aCp[%ld]=%ld", i, i + 1,
                        aChildp[i + 1], i, aChildp[i]));
        }

        for (Int i = Sleft[Super[f]]; i < Sleft[Super[f + 1]]; i++)
        {
            // number of rows
            ASSERT(i < ms);

            rM[i] = i + f;
            ASSERT(i + f < ms + nf + 1);
            aParent[i + f] = offset + numRow;
            ASSERT(childpointer < ms + nf + 1);
            aChild[childpointer++] = i + f;
        }

        offset += numRow;
        snM[f] = offset++;
        ASSERT(offset < ms + nf + 1);
        ASSERT(aChildp[offset] == -1);
        aChildp[offset] = aChildp[offset - 1] + numRow + numoforiginalChild;
        PRLEVEL(
            1, ("\n %% f=%ld numoforiginalChild=%ld\n", f, numoforiginalChild));

        if (Parent[f] == f + 1)
        {  // last child due to staircase
            PRLEVEL(1, ("%% last Child =%ld\n", f));
            lastChildFlag = 1;
        }
        else
            lastChildFlag = 0;
    }

    // Initialize first descendent of the etree
    PRLEVEL(PR, ("%% computing first of\n "));
    for (Int i = 0; i < nf; i++)
    {
        for (Int r = i; r != -1 && first[r] == -1; r = Parent[r])
        {
            PRLEVEL(PR, ("%% first of %ld is %ld\n", r, i));
            first[r] = i;
        }
    }

#ifndef NDEBUG
    PR = 1;
    PRLEVEL(PR, ("%% super node mapping ,snM (and rows): "));
    for (Int f = 0; f < nf; f++)
    {
        ASSERT(snM[f] != -1);
        PRLEVEL(PR, ("%ld (%ld) ", snM[f], Sleft[Super[f]]));
    }
    PRLEVEL(PR, ("- (%ld) ", Sleft[Super[nf]]));
    PRLEVEL(PR, ("\n"));

    PR = 1;
    PRLEVEL(PR, ("%% row mapping (rM): "));
    for (Int i = 0; i < ms; i++)
    {
        ASSERT(rM[i] != -1);
        PRLEVEL(PR, ("%ld ", rM[i]));
    }
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("%% aParent: "));
    for (Int i = 0; i < ms + nf; i++) PRLEVEL(PR, ("%ld ", aParent[i]));
    PRLEVEL(PR, ("%% \n"));

    PRLEVEL(PR, ("%% aChildp: "));
    for (Int i = 0; i < ms + nf + 2; i++) PRLEVEL(PR, ("%ld ", aChildp[i]));
    PRLEVEL(PR, ("\n"));

    PRLEVEL(PR, ("%% aChild: "));
    for (Int i = 0; i < ms + nf + 1; i++) PRLEVEL(PR, ("%ld ", aChild[i]));
    PRLEVEL(PR, ("\n"));

    // HERE
    PR = 1;
    for (Int i = 0; i < ms + nf; i++)
    {
        if (aChildp[i] == aChildp[i + 1]) continue;
        PRLEVEL(PR, ("%% anode:%ld", i));
        for (Int c = aChildp[i]; c < aChildp[i + 1]; c++)
            PRLEVEL(PR, (" %ld,", aChild[c]));
        PRLEVEL(PR, ("\n"));
    }
    PR = 1;

    PRLEVEL(PR, ("%% first: "));
    for (Int i = 0; i < nf + 1; i++)
        PRLEVEL(PR, ("first[%ld]=%ld ", i, first[i]));
    PRLEVEL(PR, ("\n"));

#endif
    LUsym->my_time = omp_get_wtime() - my_start_time;
    paru_free(m, sizeof(Int), Pinv);
    return (LUsym);
}
