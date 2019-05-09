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

    DEBUGLEVEL(1);
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
                  //       |    ***xx\xxxxx|            \                  
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
                  // NOTE: SPQR uses PLinv for stairecase structure and that is
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
                                      //
    /* ---------------------------------------------------------------------- */
    /*    Setting up umfpakc symbolic analysis and do the  analysis phase     */
    /* ---------------------------------------------------------------------- */

    /* get the default control parameters */
    umfpack_dl_defaults (Control) ;

    /* print the control parameters */
    umfpack_dl_report_control (Control) ;

    /* performing the symbolic analysis */ 
    status = umfpack_dl_symbolic (m, n, Ap, Ai, Ax, &Symbolic, Control, Info);

    if (status < 0){
        umfpack_dl_report_info (Control, Info);
        umfpack_dl_report_status (Control, status);
        printf ("umfpack_dl_symbolic failed");
        return NULL;
    }

    printf ("\nSymbolic factorization of A: ");
    (void) umfpack_dl_report_symbolic (Symbolic, Control);

    /* ---------------------------------------------------------------------- */
    /*    Copy the contents of Symbolic in my data structure                  */
    /* ---------------------------------------------------------------------- */
    Pinit = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Qinit = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Front_npivcol = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Front_1strow = (Int *) paru_alloc ((n+1), sizeof (Int) , cc);
    Front_leftmostdesc = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Front_parent = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Chain_maxrows = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
    Chain_start = (Int *) paru_alloc ((n+1), sizeof (Int), cc);
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
        return NULL;
    }
#ifndef NDEBUG
    Int p = 1;

    PRLEVEL (p, ("%%%% n1 is %d\n", n1 ));
    PRLEVEL (p, ("From the Symbolic object,\
                C is of dimension %ld-by-%ld\n", nr, nc));
    PRLEVEL (p, ("   with nz = %ld, number of fronts = %ld,\n", anz, nfr));
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

    paru_free ((n+1), sizeof (Int), Pinit, cc);
    paru_free ((n+1), sizeof (Int), Qinit, cc);
    paru_free ((n+1), sizeof (Int) , Front_1strow, cc);
    paru_free ((n+1), sizeof (Int), Front_leftmostdesc, cc);
    paru_free ((n+1), sizeof (Int), Chain_start, cc) ;
    paru_free ((n+1), sizeof (Int), Chain_maxrows, cc);
    paru_free ((n+1), sizeof (Int), Chain_maxcols, cc);


    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    ASSERT ( m == nr); 
    ASSERT ( n == nc); 

    LUsym->m = m;
    LUsym->n = n;
    LUsym->anz = anz;
    Int nf =  LUsym->nf = nfr;
    //    rjsize =  LUsym->rjsize = QRsym->rjsize;
    //

    PRLEVEL (0, ("%% A  is  %ld x %ld \n",m, n ));
    PRLEVEL (-1, ("LU = zeros(%ld,%ld);\n",m, n ));
    PRLEVEL (-1, ("npivots =[]; \n" ));
    PRLEVEL (-1, ("S = zeros(%ld,%ld);\n",m, n ));
    PRLEVEL (0, ("%% nf=%ld\n",nf ));
    //    PRLEVEL (0, ("%% anz = %ld  rjsize=%ld\n", anz, rjsize));
    //

    //    LUsym->maxfn = QRsym->maxfn;

    Int  *Qfill, *PLinv;
    //brain transplant

    // Parent size is nf+1 potentially smaller than what UMFPACK allocate
    Int size = n + 1;
    Int *Parent = 
        (Int*) paru_realloc (nf+1, sizeof(Int), Front_parent, &size, cc);
    LUsym->Parent = Parent  = Front_parent;

    // Making Super like SPQR: Super[f]<= pivotal columns of (f) < Super[f+1]
    Int *Super =  LUsym->Super = (Int *) paru_alloc ((nf+1), sizeof (Int), cc); 
    Super [0] = 0;
    for(Int k = 1; k <= nf ; k++){
        Super[k] = Front_npivcol[k-1];
    }
    paru_cumsum (nf+1, Super);
    paru_free ((n+1), sizeof (Int), Front_npivcol, cc);

    //Making Children list
    Int *Childp = (Int *) paru_calloc ((nf+2), sizeof (Int), cc);
    LUsym->Childp =  Childp;
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
    //copy of Childp
    Int *cChildp = (Int *) paru_alloc ((nf+2), sizeof (Int), cc); 
    memcpy(cChildp, Childp, nf+2);
    for (Int f = 0; f < nf+2; f++){
   //     ASSERT (cChildp[f] == Childp[f]);
        cChildp[f] = Childp[f];
    }
#ifndef NDEBUG
    p = 1;
    PRLEVEL (p, ("%%%%After Copy_cChidlp_____"));
    for (Int f = 0; f < nf+2; f++){
        PRLEVEL (p, ("%ld ",  cChildp[f]));
    }
    PRLEVEL (p, ("\n"));
#endif


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

    for (Int f = 0; f < nf; f++){
       if (cChildp[f] != Childp[f+1])
           PRLEVEL (p, ("%% cChildp [%ld] == %ld & Chldp[%ld+1] = %ld",
                       f, cChildp[f], f , Childp[f+1]));
       ASSERT (cChildp[f] == Childp[f+1]);
    }
#endif 
    paru_free ((nf+2), sizeof (Int), cChildp, cc);
       
    
    LUsym->Qfill = NULL;
    //    Qfill =  LUsym->Qfill =  QRsym->Qfill;  
    //    QRsym->Qfill = NULL;
    
    LUsym->PLinv = NULL;
    //    PLinv =  LUsym->PLinv =  QRsym->PLinv;  
    //    QRsym->PLinv = NULL;


    LUsym->Fm= NULL;
    //    LUsym->Fm = QRsym->Fm; 
    //    QRsym->Fm = NULL;


    LUsym->Cm= NULL;
    //    LUsym->Cm = QRsym->Cm; 
    //    QRsym->Cm = NULL;


    LUsym->Rj= NULL;
    //    LUsym->Rj = QRsym->Rj; 
    //    QRsym->Rj = NULL;


    LUsym->Rp= NULL;
    //    LUsym->Rp = QRsym->Rp; 
    //    QRsym->Rp = NULL;
    //
    //Staircase structure
    // Int *Sp, *Sj, *Sleft;
    LUsym->Sp= NULL;
    // Sp =  LUsym->Sp = QRsym->Sp;     
    // QRsym->Sp = NULL;
    LUsym->Sj= NULL;
    // Sj =  LUsym->Sj = QRsym->Sj;     
    // QRsym->Sj = NULL;
    LUsym->Sleft= NULL;
    // Sleft = LUsym->Sleft = QRsym->Sleft;  
    // QRsym->Sleft = NULL;


    /*Computing augmented tree */
    Int *aParent = LUsym->aParent = NULL; //augmented tree size m+nf
    Int *aChildp = LUsym->aChildp = NULL;  // size m+nf+2
    Int *aChild = LUsym->aChild = NULL;// size m+nf+1
    Int *rM = LUsym->row2atree = NULL; // row map 
    Int *snM = LUsym->super2atree = NULL; //and supernode map
    Int *first= LUsym->first = NULL; //first in augmented tree
    // augmented tree size m+nf

    /*! Check if there exist empty row or column in square part*/
    // TODO: after making the staircase
    //    for (Int row = 0; row < m; row++){
    //        PRLEVEL (1,("Sprow[%ld]=%ld\n", row, Sp[row]));
    //        if (Sp [row] == Sp[row+1] ){
    //            printf("%%Empty Row\n");
    //            paru_freesym (&LUsym , cc);
    //            return NULL;
    //        }
    //    }
    //    for (Int col = 0; col < n; col++){
    //        Int *Ap =(Int*) A->p;
    //        if (Ap [col] == Ap [col+1]){
    //            printf("%%Empty Column\n");
    //            paru_freesym (&LUsym , cc);
    //            return NULL;
    //        }
    //    }
    //

#ifndef NDEBUG
    /* print fronts*/
    for (Int f = 0; f < nf; f++){
        //       Int fm = LUsym->Fm[f];
        //       Int fn = LUsym->Rp[f+1]-LUsym->Rp[f];
        Int col1 = Super[f]; 
        Int col2 = Super[f+1]; 
        Int fp = Super[f+1]-Super[f];

        Int p=1;
        PRLEVEL (p,("%% Front=%ld Parent=%ld\n", f,  Parent[f]));
        PRLEVEL (p,("%% %ld...%ld npivotal=%ld\n", col1, col2, fp));
        PRLEVEL (p,("%% list of %ld children: ", Childp[f+1]-Childp[f]));
        for (Int i = Childp[f]; i <= Childp[f+1]-1; i++) 
            PRLEVEL (p,("%ld ",Child[i]));
        PRLEVEL (p,("\n"));
    }
#endif 

    aParent = (Int*) paru_alloc (m+nf, sizeof(Int),cc);
    aChild =  (Int*) paru_alloc (m+nf+1, sizeof(Int),cc);
    aChildp = (Int*) paru_alloc (m+nf+2, sizeof(Int),cc);
    first= (Int*) paru_alloc (m+nf, sizeof(Int),cc);


    rM =  (Int*) paru_alloc(m  ,sizeof(Int), cc);
    snM = (Int*) paru_alloc(nf ,sizeof(Int), cc);

    if(aParent == NULL || aChild == NULL || aChildp == NULL ||
            rM == NULL  || snM == NULL || free == NULL){
        printf ("Out of memory in symbolic phase");
        paru_freesym (&LUsym , cc);
        return NULL;
    }

    /* ---------------------------------------------------------------------- */
    /*    computing the augmented tree and the data structures related        */
    /* ---------------------------------------------------------------------- */

    //    //initialization
    //    /*aParent needs initialization*/
    //    //memset (snM, -1, nf*sizeof(Int));
    //    //memset (rM, -1, m*sizeof(Int));
    //    memset (aParent, -1, (m+nf)*sizeof(Int));
    //    //memset (aChild, -1, (m+nf+1)*sizeof(Int));
    //    memset (aChildp, -1, (m+nf+2)*sizeof(Int));
    //    memset (first, -1, (m+nf)*sizeof(Int));
    //
    //    aChildp[0] = 0;
    //    Int offset = 0; //number of rows visited in each iteration orig front+ rows
    //    Int lastChildFlag = 0;
    //    Int childpointer = 0;
    //
    //    for (Int f = 0; f < nf; f++) {
    //        PRLEVEL (1,("%% Front %ld\n", f)) ;
    //        PRLEVEL (1,("%% pivot columns [ %ld to %ld ] n: %ld \n",
    //                    Super [f], Super [f+1]-1, n)) ;
    //        ASSERT(Super[f+1] <= n);
    //        Int numRow =Sleft[Super[f+1]]-Sleft[Super[f]] ;
    //
    //        Int numoforiginalChild=0;
    //        if (lastChildFlag){  // the current node is the parent
    //            PRLEVEL (1,("%% Childs of %ld: ",f)) ;
    //            numoforiginalChild= Childp[f+1] - Childp[f];
    //
    //            for (Int i = Childp[f]; i <= Childp[f+1]-1; i++) {
    //                PRLEVEL (1,("%ld ",Child[i]));
    //                Int c= Child [i];
    //                ASSERT(snM[c] < m+nf+1);
    //                aParent[ snM[c]]=offset+numRow;
    //                PRLEVEL (1, ("%% aParent[%ld] =%ld\n", 
    //                            aParent[snM [c]], offset+numRow));
    //                ASSERT(childpointer < m+nf+1);
    //                aChild[childpointer++] = snM[c];
    //            }
    //        }
    //
    //        PRLEVEL (1,("%% numRow=%ld ",numRow));
    //        PRLEVEL (1,("#offset=%ld\n",offset));
    //        for(Int i = offset ; i < offset+numRow ; i++){
    //            ASSERT (aChildp [i+1] == -1);
    //            aChildp[i+1] = aChildp[i];
    //            PRLEVEL (1, ("%% @i=%ld aCp[%ld]=%ld aCp[%ld]=%ld",
    //                        i, i+1, aChildp[i+1], i, aChildp[i]));
    //        }
    //
    //        for (Int i = Sleft[Super[f]]; i < Sleft[Super[f+1]]; i++){
    //            // number of rows
    //            ASSERT(i < m);
    //
    //            rM[i] = i+f;
    //            ASSERT(i+f < m+nf+1);
    //            aParent[i+f] = offset+numRow;
    //            ASSERT(childpointer < m+nf+1);
    //            aChild[childpointer++] = i+f;
    //        }
    //
    //        offset += numRow;
    //        snM[f] = offset++;
    //        ASSERT(offset < m+nf+1);
    //        ASSERT (aChildp [offset] == -1);
    //        aChildp[offset] = aChildp[offset-1]+numRow+numoforiginalChild;
    //        PRLEVEL (1, ("\n %% f=%ld numoforiginalChild=%ld\n",
    //                    f, numoforiginalChild));
    //
    //        if( Parent[f] == f+1){  //last child due to staircase
    //            PRLEVEL (1, ("%% last Child =%ld\n", f));
    //            lastChildFlag = 1;  
    //        }else
    //            lastChildFlag = 0;  
    //    }
    //
    //    LUsym->aParent = aParent;
    //    LUsym->aChildp = aChildp;
    //    LUsym->aChild = aChild;
    //    LUsym->row2atree = rM;
    //    LUsym->super2atree = snM;
    //
    //    //Initialize first descendent of augmented tree
    //    for(Int i=0 ; i<m+nf; i++){
    //        for (Int r = i; r!= -1 && first[r] == -1; r = aParent[r])
    //            first[r] = i;
    //    }
    //
    //
    //    LUsym->first= first;
    //
    //#ifndef NDEBUG
    //    Int p = 1;
    //    PRLEVEL (p,("%% super node mapping (snM): ")); 
    //    for (Int f = 0; f < nf; f++){
    //        ASSERT (snM [f] != -1);
    //        PRLEVEL (p,("%ld ",snM[f]));     
    //    }
    //    PRLEVEL (p,("\n"));
    //
    //    PRLEVEL (p,("%% row mapping (rM): "));  
    //    for (Int i = 0; i < m; i++){
    //        ASSERT (rM [i] != -1);
    //        PRLEVEL (p,("%ld ",rM[i]));        
    //    }
    //    PRLEVEL (p,("\n"));
    //
    //    PRLEVEL (p,("%% aParent: ")); 
    //    for (Int i = 0; i < m+nf; i++) PRLEVEL (p,("%ld ",aParent[i]));
    //    PRLEVEL (p,("%% \n"));
    //
    //    PRLEVEL (p,("%% aChildp: "));
    //    for (Int i = 0; i < m+nf+2; i++) PRLEVEL (p,("%ld ",aChildp[i]));
    //    PRLEVEL (p,("\n"));
    //
    //    PRLEVEL (p,("%% aChild: ")); 
    //    for (Int i = 0; i < m+nf+1; i++) PRLEVEL (p,("%ld ",aChild[i]));
    //    PRLEVEL (p,("\n"));
    //
    //    for(Int i=0; i< m+nf; i++){
    //        PRLEVEL (p,("%% anode:%ld",i));
    //        for(Int c=aChildp[i]; c< aChildp[i+1]; c++)
    //            PRLEVEL (p,(" %ld,",aChild[c]));  
    //        PRLEVEL (p,("\n"));
    //    }
    //
    //    PRLEVEL (p,("%% first: ")); 
    //    for (Int i = 0; i < m+nf; i++) 
    //        PRLEVEL (p,("first[%ld]=%ld ",i,first[i]));
    //    PRLEVEL (p,("\n"));
    //
    //#endif
    return (LUsym) ;
}
