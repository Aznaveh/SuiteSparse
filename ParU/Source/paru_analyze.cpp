/* =========================================================================   /
 * ============================== paru_analyze =============================   /
 * =========================================================================   /
 * @brief Computing etree and do the symbolic analysis. In this file I am going
 * to use umfpack symbolic analysis instaed of spqrsym. However, I will keep the
 * style of spqr mostly.
 *         
 *
 * Example: ./Matrix/problem.mtx
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
#include "Parallel_LU.hpp"
paru_symbolic *paru_analyze
(
 // inputs, not modified
 cholmod_sparse *A,
 // workspace and parameters
 cholmod_common *cc ){   

    DEBUGLEVEL(0);
    paru_symbolic *LUsym;

    LUsym = (paru_symbolic*) paru_alloc (1, sizeof(paru_symbolic), cc);
    // ... check for LUsym NULL ...
    if(LUsym == NULL){
        //out of memory
        return NULL;
    }


    Int  anz, rjsize ; 

    Int *Ap = (Int*) A->p;
    Int *Ai = (Int*) A->i;
    double *Ax = (double*) A->x;
    Int m = A->nrow;
    Int n = A->ncol;

    // Initializaing pointers with NULL; just in case for an early release
    // not to free an uninitialized space
    LUsym->Chain_start = LUsym->Chain_maxrows = LUsym->Chain_maxcols = NULL;
    LUsym->Parent = LUsym->Super = LUsym->Child = LUsym->Childp = NULL;
    LUsym->Qfill = LUsym->Pinv = NULL;
    LUsym->Sp = LUsym->Sj = LUsym->Sleft= NULL; LUsym->Sx = NULL;
    LUsym->Fm= LUsym->Cm= LUsym->Rj= LUsym->Rp= NULL;
    LUsym->aParent = LUsym->aChildp = LUsym->aChild = LUsym->row2atree = NULL;
    LUsym->super2atree = LUsym->first = NULL;
    
    

    //############  Calling UMFPACK and retrieving data structure ##############

    /* ---------------------------------------------------------------------- */
    /*    The varialbes are needed for the UMFPACK symbolic analysis phase    */
    /* ---------------------------------------------------------------------- */

    Int nr, nc,   // A is nrxnc, I will use mxn; they should be the same anyway

        n1,       // The number of pivots with zero Markowitz cost.
        // Info[UMFPACK_COL_SINGLETONS]+Info[UMFPACK_ROW_SINGLETONS]
        // They apper first in the output permutations P and Q
        //        _______________
        //       |\**************r    
        //       |  \************r -> UMFPACK_ROW_SINGLETONS
        //       |   \***********r                              
        //       |    *\         |
        //       |    ***\       |
        //       |    ***xx\xxxxx|            |                  
        //       |    ***xxxx\xxx|            +   = n1
        //       -----ccc--------            /
        //             |
        //            UMFPACK_COL_SINGLETONS

        nfr,      // The number of frontam matrices; nf in SPQR analysis

        nchains,  // The frontal matrices are related to one another by the 
        // supernodal column elimination tree. Each nod in this tree
        // is one frontal matrix. The tree is partitioned into a set
        // of disjoint paths, and a frontal matrix chaing is one path
        // in this tree.  UMFPACK uses unifrontal technique to 
        // factroize chains, with a single working array that holds
        // each frontal matrix in the chain, one at a time. nchains is
        // in the range 0 to nfr

        *Pinit,   // The inital row permutation. If P [k] = i, then this means
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

        *Qinit,   // The inital column permutation. If Q [k] = j, then this 
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

        *Front_npivcol, // size = n_col +1;  actual size = nfr+1
        // NOTE: This is not the case for SPQR 
        // I think SPQR is easier:
        // Front_npivcol [f] = Super [f] - Super [f+1]
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
        // Aznaveh: I am not sure if I need to keep it

        *Front_leftmostdesc, // size = n_col +1;  actual size = nfr+1
        // Aznaveh: I have a module computing leftmostdesc
        // for my augmented tree; so maybe do not need it

        *Chain_start,   // size = n_col +1;  actual size = nfr+1
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

    double status,   // Info [UMFPACK_STATUS] 
    Info[UMFPACK_INFO],// Contains statistics about the symbolic analysis
    
    Control [UMFPACK_CONTROL]; // it is set in umfpack_dl_defaults and
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
    umfpack_dl_defaults (Control) ;
    Control [UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;

#ifndef NDEBUG
    /* print the control parameters */
    Int p = 1;
    if (p <= 0)  umfpack_dl_report_control (Control) ;
#endif

    /* performing the symbolic analysis */ 
    //status = umfpack_dl_symbolic (m, n, Ap, Ai, Ax, &Symbolic, Control, Info);

   void *SW; 
   //const Int Quser[];
   status = umfpack_dl_azn_symbolic
      (m, n, Ap, Ai, Ax,
       NULL,    // user provided ordering
       FALSE,
       NULL,    // user params
       &Symbolic,
       &SW,
       Control, Info);



    if (status < 0){
        umfpack_dl_report_info (Control, Info);
        umfpack_dl_report_status (Control, status);
        printf ("umfpack_dl_symbolic failed");
        paru_freesym (&LUsym , cc);
        return NULL;
    }

    Int cs1 = Info[UMFPACK_COL_SINGLETONS];
    Int rs1 = Info[UMFPACK_ROW_SINGLETONS];

//    Int maxnrows = Symbolic->maxnrows;  // I can not access to inside of 
//    Int maxncols = Symbolic->maxncols;  // SymbolicType it is internal in
                                        //  UMFPACK
                                        //      
    SWType *mySW = (SWType *)SW;
    Int *Front_nrows = (Int *) mySW->Front_nrows;
    mySW->Front_nrows = NULL; //TODO is it OK? or I should copy

#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("\n%%Symbolic factorization of A: "));
    if (p <= 0) (void) umfpack_dl_report_symbolic (Symbolic, Control);
    PRLEVEL(p, ("\n %%colsingleton = %ld, rowsingleton=%ld",cs1,rs1));
    //    PRLEVEL(p, ("\n %%maxnrow= %ld, maxncol=%ld",maxnrows, maxncols));
#endif

    /* ---------------------------------------------------------------------- */
    /*    Copy the contents of Symbolic in my data structure                  */
    /* ---------------------------------------------------------------------- */
    Pinit = (Int *) paru_alloc ((m+1), sizeof (Int), cc);
    Qinit = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Front_npivcol = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Front_1strow = (Int *) paru_alloc ((n+1), sizeof (Int) , cc);
    Front_leftmostdesc = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Front_parent = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Chain_start = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Chain_maxrows = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Chain_maxcols = (Int *) paru_alloc ((n+1), sizeof (Int), cc);

    if (!Pinit || !Qinit || !Front_npivcol || !Front_parent || !Chain_start ||
            !Chain_maxrows || !Chain_maxcols || 
            !Front_1strow || !Front_leftmostdesc){
        printf("out of memory") ;
    }



    status = umfpack_dl_get_symbolic (&nr, &nc, &n1, &anz, &nfr, &nchains,
            Pinit, Qinit, Front_npivcol, Front_parent, Front_1strow,
            Front_leftmostdesc, Chain_start, Chain_maxrows, Chain_maxcols,
            Symbolic);
    if (status < 0)    {
        printf ("symbolic factorization invalid");

        //free memory
        paru_free ((m+1), sizeof (Int), Pinit, cc);
        paru_free ((n+1), sizeof (Int), Qinit, cc);
        paru_free ((n+1), sizeof (Int), Front_npivcol, cc);
        paru_free ((n+1), sizeof (Int), Front_1strow, cc);
        paru_free ((n+1), sizeof (Int), Front_leftmostdesc, cc);
        paru_free ((n+1), sizeof (Int), Front_parent, cc);
        paru_free ((n+1), sizeof (Int), Chain_start, cc);
        paru_free ((n+1), sizeof (Int), Chain_maxrows, cc);
        paru_free ((n+1), sizeof (Int), Chain_maxcols, cc);

        paru_freesym (&LUsym , cc);

        return NULL; 
    }
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%%%% n1 is %d\n", n1 ));
    PRLEVEL (p, ("From the Symbolic object,\
                C is of dimension %ld-by-%ld\n", nr, nc));
    PRLEVEL (p, ("   with nz = %ld, number of fronts = %ld,\n", anz, nfr));
    p = 1;
    PRLEVEL (p, ("   number of frontal matrix chains = %ld\n", nchains));

    PRLEVEL (1, ("\nPivot columns in each front, and parent of each front:\n"));
    Int k = 0;
    for (Int i = 0 ; i < nfr ; i++)
    {
        Int fnpiv = Front_npivcol [i];
        PRLEVEL (p, ("Front %ld: parent front: %ld number of pivot cols: %ld\n",
                    i, Front_parent [i], fnpiv));
        for (Int j = 0 ; j < fnpiv ; j++)
        {
            Int col = Qinit [k];
            PRLEVEL (p, ("%ld-th pivot column is column %ld\
                        in original matrix\n", k, col));
            k++;
        }
    }

    PRLEVEL (p, ("\nTotal number of pivot columns\
                in frontal matrices: %ld\n", k));

    PRLEVEL (p, ("\nFrontal matrix chains:\n"));
    for (Int j = 0 ; j < nchains ; j++){
        PRLEVEL (p, ("Frontal matrices %ld to %ld are factorized in a single\n",
                    Chain_start [j], Chain_start [j+1] - 1));
        PRLEVEL (p, ("        working array of size %ld-by-%ld\n",
                    Chain_maxrows [j], Chain_maxcols [j]));
    }
#endif

    umfpack_dl_free_symbolic (&Symbolic);
    UMFPACK_azn_free_sw (&SW);

    paru_free ((n+1), sizeof (Int), Front_1strow, cc);
    paru_free ((n+1), sizeof (Int), Front_leftmostdesc, cc);

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    ASSERT ( m == nr); 
    ASSERT ( n == nc); 

    LUsym->m = m;
    LUsym->n = n;
    LUsym->n1 = n1;
    LUsym->rs1= rs1;
    LUsym->cs1= cs1;
    LUsym->anz = anz;
    Int nf =  LUsym->nf = nfr;
    LUsym->Chain_start = Chain_start;
    LUsym->Chain_maxrows = Chain_maxrows;
    LUsym->Chain_maxcols = Chain_maxcols;
    Int *Qfill =  LUsym->Qfill = Qinit;
    //    rjsize =  LUsym->rjsize = QRsym->rjsize;
    //

    PRLEVEL (0, ("%% A  is  %ld x %ld \n",m, n ));
    PRLEVEL (-1, ("LU = zeros(%ld,%ld);\n",m, n ));
    PRLEVEL (-1, ("npivots =[]; \n" ));
    PRLEVEL (-1, ("S = zeros(%ld,%ld); %% n1 = %ld\n",m, n, n1 ));
    PRLEVEL (0, ("%% nf=%ld\n",nf ));
    //    PRLEVEL (0, ("%% anz = %ld  rjsize=%ld\n", anz, rjsize));
    //
    /* ---------------------------------------------------------------------- */
    /*           Fixing Parent and computing Children datat structure         */
    /* ---------------------------------------------------------------------- */

    // Parent size is nf+1 potentially smaller than what UMFPACK allocate
    Int size = n + 1;
    Int *Parent = 
        (Int*) paru_realloc (nf+1, sizeof(Int), Front_parent, &size, cc);
    ASSERT (size < n+1);
    if (Parent == NULL){    // should not happen anyway it is always shrinking
        printf ("memory problem");
        //free memory
        paru_free ((n+1), sizeof (Int), Front_npivcol, cc);
        paru_free ((n+1), sizeof (Int), Front_parent, cc);

        paru_freesym (&LUsym , cc);
        return NULL;  
    }
    LUsym->Parent = Parent; 

    // Making Super data structure
    // like SPQR: Super[f]<= pivotal columns of (f) < Super[f+1]
    Int *Super =  LUsym->Super = (Int *) paru_alloc ((nf+1), sizeof (Int), cc); 
    if (Super == NULL){   
        printf ("memory problem");
        paru_freesym (&LUsym , cc);
        return NULL;  
    }
    Super [0] = 0;
    for(Int k = 1; k <= nf ; k++){
        Super[k] = Front_npivcol[k-1];
    }
    paru_cumsum (nf+1, Super);
    paru_free ((n+1), sizeof (Int), Front_npivcol, cc);

    //Making Children list
    Int *Childp = (Int *) paru_calloc ((nf+2), sizeof (Int), cc);
    LUsym->Childp =  Childp;
    if (Childp== NULL){   
        printf ("memory problem");
        paru_freesym (&LUsym , cc);
        return NULL;  
    }

    for (Int f = 0; f < nf; f++){
        if (Parent [f] > 0)
            Childp[Parent[f]+1]++;
    }
    paru_cumsum (nf+2, Childp);
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%%%%-Chidlp-----"));
    for (Int f = 0; f < nf+2; f++){
        PRLEVEL (p, ("%ld ", Childp[f]));
    }
    PRLEVEL (p, ("\n"));
#endif 
    Int *Child = (Int *) paru_calloc ((nf+1), sizeof (Int), cc); 
    LUsym->Child  =  Child;
    if (Child == NULL){   
        printf ("memory problem");
        paru_freesym (&LUsym , cc);
        return NULL;  
    }

    //copy of Childp using Work for other places also
    Int *Work = (Int *) paru_alloc ((MAX(m,n)+2), sizeof (Int), cc); 
    Int *cChildp = Work;
    if (cChildp == NULL){   
        printf ("memory problem");
        paru_freesym (&LUsym , cc);
        return NULL;  
    }


    memcpy(cChildp, Childp, (nf+2)*sizeof(Int) );

    for (Int f = 0; f < nf; f++){
        if (Parent [f] > 0)
            Child[cChildp[Parent[f]]++] = f;
    }
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%%%%_cChidlp_____"));
    for (Int f = 0; f < nf+2; f++){
        PRLEVEL (p, ("%ld ",  cChildp[f]));
    }
    PRLEVEL (p, ("\n"));
#endif 



    /* ---------------------------------------------------------------------- */
    /*                   computing the Staircase structures                   */
    /* ---------------------------------------------------------------------- */

    Int *Sp = LUsym->Sp = (Int*) paru_calloc (m+1-n1, sizeof(Int),cc);
    Int *Sleft = LUsym->Sleft = (Int*) paru_alloc (n+2-n1, sizeof(Int),cc);
    Int *Pinv =  LUsym->Pinv =  (Int *) paru_alloc (m+1, sizeof (Int), cc);

    if (Sp == NULL || Sleft == NULL || Pinv == NULL ){   
        printf ("memory problem");
        paru_freesym (&LUsym , cc);
        return NULL;  
    }

    //-------- computing the inverse permutation for P
    for(Int i = 0; i < m ; i++){
        Pinv[Pinit[i]] = i;
    }
    
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("Qinit =\n"));
    for (Int j = 0; j < m; j++){
        PRLEVEL (p, ("%ld ", Qinit[j]));
    }
    PRLEVEL (p, ("\n"));

    PRLEVEL (p, ("Pinit =\n"));
    for (Int i = 0; i < m; i++){
        PRLEVEL (p, ("%ld ", Pinit[i]));
    }
    PRLEVEL (p, ("\n"));


    PRLEVEL (p, ("Pinv =\n"));
    for (Int i = 0; i < m; i++){
        PRLEVEL (p, ("%ld ", Pinv[i]));
    }
    PRLEVEL (p, ("\n"));

#endif



    Int *Ps; // new row permutation for just the Submatrix part
    Ps = (Int*) paru_calloc (m-n1, sizeof(Int),cc); 
    if (Ps == NULL ){   
        printf ("memory problem");
        paru_freesym (&LUsym , cc);
        return NULL;  
    }

    Int unz = 0;  // U nnz: singlteton nnzero of s
    Int snz = 0;  //s nonzero: nnz in submatrix excluding singletons
    int rowcount = 0;
    Sleft[0] = 0;
    //counting number of entries in each row of submatrix Sp
    PRLEVEL (1, ("Computing Staircase Structure\n"));
    for (Int newcol = n1; newcol < n ; newcol++){
        Int oldcol = Qinit [newcol];
        PRLEVEL (1, ("newcol = %ld oldcol=%ld\n",newcol, oldcol));
        for(Int p = Ap[oldcol] ; p < Ap[oldcol+1]; p++){
            Int oldrow = Ai[p];
            Int newrow = Pinv[oldrow]; 
            Int srow = newrow - n1;
            PRLEVEL (1, ("\tnewrow=%ld oldrow=%ld srow=%ld\n",
                        newrow, oldrow, srow));
            if (srow >= 0){  // it is insdie S otherwise it is part of singleton
                if(Sp[srow+1] == 0){ //first time seen
                    PRLEVEL (1, ("\tPs[%ld]= %ld\n",rowcount, srow));
                    Ps[rowcount] = srow;
                    Pinit[n1+rowcount] = oldrow;
                    rowcount++;
                }
                snz++;
                Sp[srow+1]++;
            } 
            else { // inside the upart
               unz++; 
            }
        }
        Sleft[newcol-n1+1] = rowcount;
    }
    Sleft[n-n1+1] = m-n1-rowcount;  //empty rows of S if any
    LUsym->snz = snz;

#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("Ps =\n"));
    for (Int k = 0; k < rowcount; k++){
        PRLEVEL (p, ("%ld ", Ps[k]));
    }
    PRLEVEL (p, ("\n"));
#endif

    if (rowcount < m-n1){
        //I think that must not happen anyway while umfpack finds it
        printf("Empty rows in submatrix\n"); 
#ifndef NDEBUG
        PRLEVEL (1, ("m = %ld, n1 = %ld, rowcount = %ld, snz = %ld\n",
                    m , n1 ,rowcount, snz));
        for(Int srow = 0; srow < m-n1; srow++){
            if(Sp[srow] == 0){
                PRLEVEL (1, ("Row %ld is empty\n",srow));
            }
        }
#endif
        paru_freesym (&LUsym , cc);
        return NULL; //Free memory
    }
    ASSERT (rowcount == m-n1);

    //Update Sp based on new permutation for Sp[1...m+1]
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("Before permutation Sp =\n"));
    for (Int i = n1; i <= m; i++){
        PRLEVEL (p, ("%ld ", Sp[i-n1]));
    }
    PRLEVEL (p, ("\n"));
#endif

    // update Pinv
    for (Int i = n1; i < m; i++){
        Pinv[Pinit[i]] = i;
    }

    

///////////////////////////////////////////////////////////////
    Int *cSp = Work;
    for (Int i = n1; i < m; i++){
        Int row = i - n1 ;
        PRLEVEL (1, ("Permutation row = %ld, Ps[row]= %ld, Sp[row]=%ld\n",
                    row, Ps[row], Sp[row]));
        cSp[row+1] = Sp[Ps[row]+1];
    }

#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("After permutation Sp =\n"));
    for (Int i = n1; i <= m; i++){
        PRLEVEL (p, ("%ld ", cSp[i-n1]));
    }
    PRLEVEL (p, ("\n"));
#endif
    memcpy(Sp, cSp, (m+1-n1)*sizeof(Int) );
    paru_cumsum (m+1-n1, Sp);
    memcpy(cSp, Sp, (m+1-n1)*sizeof(Int) );

#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("After cumsum  Sp =\n"));
    for (Int i = n1; i <= m; i++){
        PRLEVEL (p, ("%ld ", Sp[i-n1]));
    }
    PRLEVEL (p, ("\n"));
#endif


#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("After Stair case Pinv =\n"));
    for (Int i = 0; i < m; i++){
        PRLEVEL (p, ("%ld ", Pinv[i]));
    }
    PRLEVEL (p, ("\n"));
#endif

    paru_free ((m+1), sizeof (Int), Pinit, cc);
///////////////////////////////////////////////////////////////
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("Before Sj Sp =\n"));
    for (Int i = n1; i <= m; i++){
        PRLEVEL (p, ("%ld ", cSp[i-n1]));
    }
    PRLEVEL (p, ("\n"));
#endif

    //Updating Sj and Sx using copy of Sp
    Int *Sj = (Int*) paru_alloc (snz, sizeof(Int),cc); 
    double *Sx = (double*) paru_alloc (snz, sizeof(double),cc); 
    LUsym->Sj = Sj; LUsym->Sx = Sx;

    if (Sj == NULL || Sx == NULL){   
        printf ("memory problem");
        paru_freesym (&LUsym , cc);
        return NULL;  
    }


    //construct Sj
    PRLEVEL (1, ("Constructing Sj\n"));
    for (Int newcol = n1; newcol < n ; newcol++){
        Int oldcol = Qinit [newcol];
        PRLEVEL (1, ("newcol = %ld oldcol=%ld\n",newcol, oldcol));
        for(Int p = Ap[oldcol] ; p < Ap[oldcol+1]; p++){
            Int oldrow = Ai[p];
            Int newrow = Pinv[oldrow]; 
            Int srow = newrow - n1;
            Int scol = newcol - n1;
            PRLEVEL (1, ("\tnewrow=%ld oldrow=%ld srow=%ld \n",
                        newrow, oldrow, srow));
            if (srow >= 0){  // it is insdie S otherwise it is part of singleton
                Sj[cSp[srow]] = scol;
                Sx[cSp[srow]++] = Ax[p];
            } 
        }
    }

    paru_free ( (m-n1), sizeof (Int), Ps, cc);
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("Sp =\n"));
    for (Int i = n1; i <= m; i++){
        PRLEVEL (p, ("%ld ", Sp[i-n1]));
    }
    PRLEVEL (p, ("\n"));

    PRLEVEL (p, ("Sj =\n"));
    for (Int k = 0; k < snz; k++){
        PRLEVEL (p, ("%ld ", Sj[k]));
    }
    p = 1;
    PRLEVEL (p, ("\n"));
    for (Int i = 0; i < m-n1; i++){
        PRLEVEL (p, (" i=%ld cSp=%ld Sp=%ld\n",i, cSp[i],Sp[i+1]));
        ASSERT (cSp[i] == Sp[i+1]);
    }
#endif
    paru_free ((MAX(m,n)+2), sizeof (Int), Work, cc);

    /* ---------------------------------------------------------------------- */
    /*         Finding the Upper bound of rows and cols                       */
    /* ---------------------------------------------------------------------- */
    //TODO: it must computed somehow from UMFPACK

    LUsym->Fm= NULL;
    //    LUsym->Fm = QRsym->Fm; 
    //    QRsym->Fm = NULL;


    LUsym->Cm= NULL;    //It seems I have never used it anyway
    //    LUsym->Cm = QRsym->Cm; 
    //    QRsym->Cm = NULL;


    LUsym->Rj= NULL;
    //    LUsym->Rj = QRsym->Rj; 
    //    QRsym->Rj = NULL;


    LUsym->Rp= NULL;
    //    LUsym->Rp = QRsym->Rp; 
    //    QRsym->Rp = NULL;


    /* ---------------------------------------------------------------------- */
    /*    computing the augmented tree and the data structures related        */
    /* ---------------------------------------------------------------------- */

    /*Computing augmented tree */
    Int *aParent = LUsym->aParent = NULL; //augmented tree size m+nf
    Int *aChildp = LUsym->aChildp = NULL;  // size m+nf+2
    Int *aChild = LUsym->aChild = NULL;// size m+nf+1
    Int *rM = LUsym->row2atree = NULL; // row map 
    Int *snM = LUsym->super2atree = NULL; //and supernode map
    Int *first= LUsym->first = NULL; //first descendent in augmented tree
    // augmented tree size m+nf


#ifndef NDEBUG
    p = 1;
    /* print fronts*/
    for (Int f = 0; f < nf; f++){
        //       Int fm = LUsym->Fm[f];
        //       Int fn = LUsym->Rp[f+1]-LUsym->Rp[f];
        Int col1 = Super[f]; 
        Int col2 = Super[f+1]; 
        Int fp = Super[f+1]-Super[f];
        PRLEVEL (p,("%% Front=%ld Parent=%ld\n", f,  Parent[f]));
        PRLEVEL (p,("%% %ld...%ld npivotal=%ld\n", col1, col2, fp));
        PRLEVEL (p,("%% list of %ld children: ", Childp[f+1]-Childp[f]));
        for (Int i = Childp[f]; i <= Childp[f+1]-1; i++) 
            PRLEVEL (p,("%ld ",Child[i]));
        PRLEVEL (p,("\n"));
    }
#endif 

    Int ms = m-n1;  Int ns = n-n1; //submatrix is msxns



    LUsym->aParent = aParent = (Int*) paru_alloc (ms+nf, sizeof(Int),cc);
    LUsym->aChild = aChild =  (Int*) paru_alloc (ms+nf+1, sizeof(Int),cc);
    LUsym->aChildp = aChildp = (Int*) paru_alloc (ms+nf+2, sizeof(Int),cc);
    LUsym->first = first= (Int*) paru_alloc (ms+nf, sizeof(Int),cc);
    LUsym->row2atree = rM =  (Int*) paru_alloc(ms  ,sizeof(Int), cc);
    LUsym->super2atree = snM = (Int*) paru_alloc(nf ,sizeof(Int), cc);


    if(aParent == NULL || aChild == NULL || aChildp == NULL ||
            rM == NULL  || snM == NULL || free == NULL){
        printf ("Out of memory in symbolic phase");
        paru_freesym (&LUsym , cc);
        return NULL;
    }
    //initialization
    memset (aParent, -1, (ms+nf)*sizeof(Int));
    //memset (aChild, -1, (ms+nf+1)*sizeof(Int));  //doens't need
    memset (aChildp, -1, (ms+nf+2)*sizeof(Int));
    memset (first, -1, (ms+nf)*sizeof(Int));

    aChildp[0] = 0;
    Int offset = 0; //number of rows visited in each iteration orig front+ rows
    Int lastChildFlag = 0;
    Int childpointer = 0;

    for (Int f = 0; f < nf; f++) {
        PRLEVEL (1,("%% Front %ld\n", f)) ;
        PRLEVEL (1,("%% pivot columns [ %ld to %ld ] n: %ld \n",
                    Super [f], Super [f+1]-1, ns)) ;
        ASSERT(Super[f+1] <= ns);
        Int numRow =Sleft[Super[f+1]]-Sleft[Super[f]] ;

        Int numoforiginalChild=0;
        if (lastChildFlag){  // the current node is the parent
            PRLEVEL (1,("%% Childs of %ld: ",f)) ;
            numoforiginalChild= Childp[f+1] - Childp[f];

            for (Int i = Childp[f]; i <= Childp[f+1]-1; i++) {
                PRLEVEL (1,("%ld ",Child[i]));
                Int c= Child [i];
                ASSERT(snM[c] < ms+nf+1);
                aParent[ snM[c]]=offset+numRow;
                PRLEVEL (1, ("%% aParent[%ld] =%ld\n", 
                            aParent[snM [c]], offset+numRow));
                ASSERT(childpointer < ms+nf+1);
                aChild[childpointer++] = snM[c];
            }
        }

        PRLEVEL (1,("%% numRow=%ld ",numRow));
        PRLEVEL (1,("#offset=%ld\n",offset));
        for(Int i = offset ; i < offset+numRow ; i++){
            ASSERT (aChildp [i+1] == -1);
            aChildp[i+1] = aChildp[i];
            PRLEVEL (1, ("%% @i=%ld aCp[%ld]=%ld aCp[%ld]=%ld",
                        i, i+1, aChildp[i+1], i, aChildp[i]));
        }

        for (Int i = Sleft[Super[f]]; i < Sleft[Super[f+1]]; i++){
            // number of rows
            ASSERT(i < ms);

            rM[i] = i+f;
            ASSERT(i+f < ms+nf+1);
            aParent[i+f] = offset+numRow;
            ASSERT(childpointer < ms+nf+1);
            aChild[childpointer++] = i+f;
        }

        offset += numRow;
        snM[f] = offset++;
        ASSERT(offset < ms+nf+1);
        ASSERT (aChildp [offset] == -1);
        aChildp[offset] = aChildp[offset-1]+numRow+numoforiginalChild;
        PRLEVEL (1, ("\n %% f=%ld numoforiginalChild=%ld\n",
                    f, numoforiginalChild));

        if( Parent[f] == f+1){  //last child due to staircase
            PRLEVEL (1, ("%% last Child =%ld\n", f));
            lastChildFlag = 1;  
        }else
            lastChildFlag = 0;  
    }

    //Initialize first descendent of augmented tree
    for(Int i=0 ; i < ms+nf; i++){
        for (Int r = i; r!= -1 && first[r] == -1; r = aParent[r])
            first[r] = i;
    }


#ifndef NDEBUG
    p = 1;
    PRLEVEL (p,("%% super node mapping (snM): ")); 
    for (Int f = 0; f < nf; f++){
        ASSERT (snM [f] != -1);
        PRLEVEL (p,("%ld ",snM[f]));     
    }
    PRLEVEL (p,("\n"));

    PRLEVEL (p,("%% row mapping (rM): "));  
    for (Int i = 0; i < ms; i++){
        ASSERT (rM [i] != -1);
        PRLEVEL (p,("%ld ",rM[i]));        
    }
    PRLEVEL (p,("\n"));

    PRLEVEL (p,("%% aParent: ")); 
    for (Int i = 0; i < ms+nf; i++) PRLEVEL (p,("%ld ",aParent[i]));
    PRLEVEL (p,("%% \n"));

    PRLEVEL (p,("%% aChildp: "));
    for (Int i = 0; i < ms+nf+2; i++) PRLEVEL (p,("%ld ",aChildp[i]));
    PRLEVEL (p,("\n"));

    PRLEVEL (p,("%% aChild: ")); 
    for (Int i = 0; i < ms+nf+1; i++) PRLEVEL (p,("%ld ",aChild[i]));
    PRLEVEL (p,("\n"));

    for(Int i=0; i< ms+nf; i++){
        PRLEVEL (p,("%% anode:%ld",i));
        for(Int c=aChildp[i]; c< aChildp[i+1]; c++)
            PRLEVEL (p,(" %ld,",aChild[c]));  
        PRLEVEL (p,("\n"));
    }

    PRLEVEL (p,("%% first: ")); 
    for (Int i = 0; i < ms+nf; i++) 
        PRLEVEL (p,("first[%ld]=%ld ",i,first[i]));
    PRLEVEL (p,("\n"));

#endif
    return (LUsym) ;
}
