/* ========================================================================== */
/* =======================  paru_init_rowFronts  ============================ */
/* ========================================================================== */
/* Assembling row fronts; other fronts will be assembled from their childs.
 *      Initializing Row and column list: 
 *          allocating memory and updating lists and initializing matrix 
 *          structre
 *      Assemble each front:
 *          Adding numerical values, allocating data
 *          updating the list
 * 
 * Updating row and column list and allocating memory for fronts will
 *  be implemnted in different files and will be reused for other fronts      */

#include "Parallel_LU.hpp"
void paru_init_rowFronts(
        // inputs, not modified
        cholmod_sparse *A,
        //symbolic analysis
        paru_symbolic *LUsym,
        // workspace and parameters
        cholmod_common *cc
){

    paru_matrix *Amat;
    Amat = (paru_matrix*) paralloc (1,sizeof(paru_matrix),cc);
    if (Amat == NULL){   //out of memory
        return;
    }

    Int m,n;     Int slackRow=2,slackCol=2;
    m= Amat->m= LUsym->m;   n= Amat->n= LUsym->n; 
    tupleList *RowList,*ColList;
    RowList=Amat->RowList =
        (tupleList*) paralloc (slackRow, m*sizeof(tupleList), cc);
    ColList=Amat->ColList =
        (tupleList*) paralloc (slackCol, n*sizeof(tupleList), cc);
    Element **elementList; Int nf = LUsym->nf;
    elementList = Amat->elementList = // Initialize with NULL
        (Element**) parcalloc (1, (m+nf+1)*sizeof(Element), cc);

    /// -------------------------------------------------------------------------
    // create S = A (p,q)', or S=A(p,q) if S is considered to be in row-form
    // -------------------------------------------------------------------------
    Int *Qfill = LUsym->Qfill;
    Int *PLinv = LUsym->PLinv;
    Int anz = LUsym->anz;
    Int *Wi = (Int *) cc->Iwork ;   // size m, aliased with the rest of Iwork

    // create numeric values of S = A(p,q) in row-form in Sx
    /*! TODO: dont forget to free all allocated memories	 */
    double *Sx = (double*) cholmod_l_malloc (anz, sizeof (double), cc) ;
    Int *W = (Int*) cholmod_l_malloc (m, sizeof (Int), cc) ;
    Int *Sp = LUsym->Sp;
    Int *Ai = (Int*) A->i;
    Int *Ap = (Int*) A->p;
    double *Ax = (double*) A->x;

    /*  I can either run the function or copy it here */
    //    if (cc->status == CHOLMOD_OK)
    //    {
    //        // use Wi as workspace (Iwork (0:m-1)) [
    //        spqr_stranspose2 (A, Qfill, Sp, PLinv, Sx, Wi) ;
    //        // Wi no longer needed ]
    //    }
    for (Int row = 0 ; row < m ; row++)
    {
        W [row] = Sp [row] ;
    }

    for (Int col = 0 ; col < n ; col++)     // for each column of A(:,Qfill)
    {
        Int j = Qfill ? Qfill [col] : col ; // col of S is column j of A
        Int pend = Ap [j+1] ;
        for (Int p = Ap [j] ; p < pend ; p++)
        {
            Int i = Ai [p] ;                // the entry A(i,j)
            Int row = PLinv [i] ;           // row of S is row i of A
            Int s = W [row]++ ;             // place S(row,col) in position
            Sx [s] = Ax [p] ;
        }
    }


    /*! TODO: Substitute C matrix with S
     * Sj?*/
    Int *Sj= LUsym->Sj;
    cholmod_sparse *C = cholmod_l_transpose (A, 1, cc);
    double *Cx;
    Int *Cp, *Ci, *Cnz;
    Int p, ncolC, xtype;
    ncolC = C->ncol;        Cp = (Int*) C->p;
    Ci = (Int*) C->i;       Cx = (double*) C->x;
    Cnz = (Int*) C->nz;     xtype = C->xtype;
    ASSERT(ncolC == m);
    if(xtype != CHOLMOD_REAL){
        printf("Error: Input matrix is not double\n");
        return; // Error: Just working with real for now
    }

    for(Int i = 0; i < m ; i++){
        Int e = LUsym->row2atree[i]; //element number in augmented tree
        Int nrows = 1,
            ncols = Cnz[i]; 

        Element *curEl = elementList[e] =
            (Element*) paralloc(1,
                    sizeof(Element)+sizeof(Int)*(nrows+ncols)+ 
                    sizeof(double)*nrows*ncols, cc);
        curEl->nrowsleft = curEl->nrows = nrows;
        curEl->ncolsleft = curEl->ncols = ncols;

        Int *colrowIndex = (Int*)(curEl+1);
        double *colrowNum = (double*)(colrowIndex+nrows+ncols);
        Int j = 0;
        /* TODO: Initializing tuple list */
        p = Cp [i];
        Int pend = p + Cnz[i];
        for ( ; p < pend ; p++){
            colrowIndex[j] = Ci[p];
            colrowNum[j++] =   Cx[p];
            /* TODO Add col tuples:  <06-06-17, Me> */
        }
        colrowIndex[j++] = i;  //initializing row one item  
        /*! TODO: add Row tuple	 */

    }
}
