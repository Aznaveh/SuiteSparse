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

    paru_matrix *paruMatInfo;
    paruMatInfo = (paru_matrix*) paru_alloc (1,sizeof(paru_matrix),cc);
    if (paruMatInfo == NULL){   //out of memory
        printf("Out of memory: paruMatInfo\n");
        return;
    }

    
    Int m,n;  
    m = paruMatInfo->m = LUsym->m;   
    n = paruMatInfo->n = LUsym->n; 

    tupleList *RowList,*ColList;
    RowList= paruMatInfo->RowList =
        (tupleList*) paru_alloc (1, m*sizeof(tupleList), cc);
    if (RowList == NULL){   //out of memory
        printf("Out of memory: RowList\n");
        return;
    }

    ColList= paruMatInfo->ColList =
        (tupleList*) paru_alloc (1, n*sizeof(tupleList), cc);
    if (ColList== NULL){   //out of memory
        printf("Out of memory: ColList\n");
        return;
    }


    Element **elementList; Int nf = LUsym->nf;
    elementList = paruMatInfo->elementList = // Initialize with NULL
        (Element**) paru_calloc (1, (m+nf+1)*sizeof(Element), cc);
    if (elementList == NULL){   //out of memory
        printf("Out of memory: elementList\n");
        return;
    }


    /// -------------------------------------------------------------------------
    // create S = A (p,q)', or S=A(p,q) if S is considered to be in row-form
    // -------------------------------------------------------------------------
    Int *Qfill = LUsym->Qfill;
    Int *PLinv = LUsym->PLinv;
    Int anz = LUsym->anz;
    Int *Wi = (Int *) cc->Iwork ;   // size m, aliased with the rest of Iwork

    // create numeric values of S = A(p,q) in row-form in Sx
    double *Sx = (double*) cholmod_l_malloc (anz, sizeof (double), cc) ;
    if (Sx == NULL){   //out of memory
        printf("Out of memory: Sx\n");
        return;
    }
    Int *W = (Int*) cholmod_l_malloc (m, sizeof (Int), cc) ;
    if (W == NULL){   //out of memory
        printf("Out of memory: W\n");
        return;
    }
    Int *Sp = LUsym->Sp;
    Int *Ai = (Int*) A->i;
    Int *Ap = (Int*) A->p;
    double *Ax = (double*) A->x;

    /*  I can either run the function or copy the content here */
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

    Int *Sj= LUsym->Sj;

    //      I am not using CHOLMOD for transposing
    //    cholmod_sparse *C = cholmod_l_transpose (A, 1, cc);
    //    double *Cx;
    //    Int *Cp, *Ci, *Cnz;
    //    Int p, ncolC, xtype;
    //    ncolC = C->ncol;        Cp = (Int*) C->p;
    //    Ci = (Int*) C->i;       Cx = (double*) C->x;
    //    Cnz = (Int*) C->nz;     xtype = C->xtype;
    //    ASSERT(ncolC == m);
    //    if(xtype != CHOLMOD_REAL){
    //        printf("Error: Input matrix is not double\n");
    //        return; // Error: Just working with real for now
    //    }

    //constants for initialzing lists
    Int slackRow = paruMatInfo->slackRow = 2;
    Int slackCol = paruMatInfo->slackCol = 2;

        /* TODO: Initializing tuple list */
    /* allocating column list must happend beforehand
     *  I think I should use Qfill inverse 
     *  I need number of column members of S
     *  while S is row based I am using A so 
     */
    for (Int col = 0 ; col < n ; col++){     //allocating Column list tuple 
        Int j = Qfill ? Qfill [col] : col ;  // col of S is column j of A
        Int ncols = Ap[j+1]-Ap[j];
         ColList[col].list = 
             (Tuple*) paru_alloc (1, slackCol*ncols*sizeof(Tuple), cc);
         if (ColList[col].list == NULL){   //out of memory
             printf("Out of memory: ColList[col].list\n");
             return;
         }

    }


    for(Int row = 0; row < m ; row++){
        Int e = LUsym->row2atree[row]; //element number in augmented tree
        Int nrows = 1,
            // ncols = Cnz[row]; 
            ncols = Sp[row+1]-Sp[row];

        Element *curEl = elementList[e] =
            (Element*) paru_alloc(1,
                    sizeof(Element)+sizeof(Int)*(nrows+ncols)+ 
                    sizeof(double)*nrows*ncols, cc);
         if (curEl == NULL){   //out of memory
             printf("Out of memory: curEl\n");
             return;
         }


        curEl->nrowsleft = curEl->nrows = nrows;
        curEl->ncolsleft = curEl->ncols = ncols;

        Int *colrowIndex = (Int*)(curEl+1);
        double *colrowNum = (double*)(colrowIndex+nrows+ncols);
        Int j = 0;

        RowList[row].list =
            (Tuple*) paru_alloc (1, slackRow*nrows*sizeof(Tuple), cc);
         if (RowList[row].list == NULL){   //out of memory
             printf("Out of memory: RowList[row].list \n");
             return;
         }


        Tuple curTuple;
        curTuple.e = e;  
        curTuple.f = 1;  
        RowList[row].list[1] = curTuple; /* structure assignment          */
        /* there is no pointer in Tuple  */
        RowList[row].numTuple = 1;
        RowList[row].len = slackRow;

        Int p = Sp [row];
        Int pend = p + ncols;
        for ( ; p < pend ; p++){
            colrowIndex[j] = Sj[p];
            colrowNum[j++] =   Sx[p];
            /* TODO Add col tuples:  <06-06-17, Me> */
        }
        colrowIndex[j++] = row;  //initializing row one item  
        /*! TODO: add Row tuple	 */

    }

    paru_free (anz, sizeof (double), Sx , cc) ;
    paru_free (m, sizeof (Int), W, cc) ;
}
