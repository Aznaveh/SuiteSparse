/** =========================================================================  /
 * =======================  paru_init_rowFronts  ============================  /
 * ==========================================================================  /
 * Initializing row fronts; fronts will be assembled later.
 *      Initializing Row and column tuple lists: 
 *          allocating memory and updating lists and initializing matrix 
 *          structre
 *      Assemble each front:
 *          Adding numerical values, allocating data
 *          updating the list
 *          
 *  be implemnted in different files and will be reused for other fronts     
 *  */
#include "Parallel_LU.hpp"

void paru_init_rowFronts (
        // inputs, not modified
        cholmod_sparse *A,
        //symbolic analysis
        paru_symbolic *LUsym,
        // workspace and parameters
        cholmod_common *cc
        ){
    if (!A->packed){
        printf ("A is not packed; Wrong format \n");
        return;
    }

    paru_matrix *paruMatInfo;
    paruMatInfo = (paru_matrix*) paru_alloc (1,sizeof(paru_matrix),cc);
    if (paruMatInfo == NULL){   //out of memory
        printf ("Out of memory: paruMatInfo\n");
        return;
    }


    Int m,n;  
    m = paruMatInfo->m = LUsym->m;   
    n = paruMatInfo->n = LUsym->n; 

    // RowList, ColList and elementList are place holders 
    // pointers to pointers that are allocated

    tupleList *RowList= paruMatInfo->RowList =
        (tupleList*) paru_alloc (1, m*sizeof(tupleList), cc);
    if (RowList == NULL){   //out of memory
        paru_freemat (&paruMatInfo, cc);
        printf ("Out of memory: RowList\n");
        return;
    }

    tupleList *ColList= paruMatInfo->ColList =
        (tupleList*) paru_alloc (1, n*sizeof(tupleList), cc);
    if (ColList== NULL){   //out of memory
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory: ColList\n");
        return;
    }

    Element **elementList; 
    Int nf = LUsym->nf;
    elementList = paruMatInfo->elementList = // Initialize with NULL
        (Element**) paru_calloc (1, (m+nf+1)*sizeof(Element), cc);
    if (elementList == NULL){   //out of memory
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory: elementList\n");
        return;
    }


    /// ------------------------------------------------------------------------
    // create S = A (p,q)', or S=A(p,q) if S is considered to be in row-form
    // -------------------------------------------------------------------------
    Int *Qfill = LUsym->Qfill;
    Int *PLinv = LUsym->PLinv;
    Int anz = LUsym->anz;
    Int *Wi = (Int *) cc->Iwork ;   // size m, aliased with the rest of Iwork

    // create numeric values of S = A(p,q) in row-form in Sx
    double *Sx = (double*) cholmod_l_malloc (anz, sizeof (double), cc) ;
    if (Sx == NULL){   //out of memory
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory: Sx\n");
        return;
    }
    Int *W = (Int*) cholmod_l_malloc (m, sizeof (Int), cc) ;
    if (W == NULL){   //out of memory
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory: W\n");
        return;
    }

    Int *Sp = LUsym->Sp;
    if (cc->status == CHOLMOD_OK){
        // use Wi as workspace (Iwork (0:m-1)) [
        spqr_stranspose2 (A, Qfill, Sp, PLinv, Sx, Wi) ;
        // Wi no longer needed ]
    }

    Int *Sj= LUsym->Sj;
    Int *Ap = (Int*) A->p;
    //constants for initialzing lists
    Int slackRow = 2;
    Int slackCol = 2;

    /* allocating column list must happend beforehand
     *  I need number of column members of S
     *  while S is row based I am using A so 
     */
    for (Int col = 0 ; col < n ; col++){     //allocating Column list tuple 
        Int j = Qfill ? Qfill [col] : col ;  // col of S is column j of A
        Int ncols = Ap[j+1]-Ap[j];
        ColList[col].numTuple = 0;
        ColList[col].len = slackCol*ncols;
        ColList[col].list = 
            (Tuple*) paru_alloc (1, slackCol*ncols*sizeof(Tuple), cc);
        if (ColList[col].list == NULL){   //out of memory
            paru_freemat (&paruMatInfo, cc);
            printf("Out of memory: ColList[col].list\n");
            return;
        }
    }

    // allocating row tuples, elements and updating column tuples
    for(Int row = 0; row < m ; row++){  // row is number of row in A
        Int e = LUsym->row2atree[row]; //element number in augmented tree
        Int nrows = 1, 
            ncols = Sp[row+1]-Sp[row]; //nrows and ncols of current front/row

        Element *curEl = elementList[e] =
            (Element*) paru_alloc(1, sizeof(Element)+sizeof(Int)*(nrows+ncols)+
                    sizeof(double)*nrows*ncols, cc);
        if (curEl == NULL){   //out of memory
            paru_freemat (&paruMatInfo, cc);
            printf("Out of memory: curEl\n");
            return;
        }
        curEl->nrowsleft = curEl->nrows = nrows;
        curEl->ncolsleft = curEl->ncols = ncols;

        // Allocating Rowlist and updating its tuples
        RowList[row].list =
            (Tuple*) paru_alloc (1, slackRow*nrows*sizeof(Tuple), cc);
        if (RowList[row].list == NULL){   //out of memory
            paru_freemat (&paruMatInfo, cc);
            printf("Out of memory: RowList[row].list \n");
            return;
        }
        RowList[row].numTuple = 0;
        RowList[row].len = slackRow;

        Tuple rowTuple;
        rowTuple.e = e;
        rowTuple.f = 1;
        if ( !paru_add_rowTuple (RowList, row, rowTuple, cc) ){
            paru_freemat (&paruMatInfo, cc);
            printf("Out of memory: add_rowTuple \n");
            return;
        }

        //Allocating elements, and updating column tuple list

        Int *el_colrowIndex = (Int*)(curEl+1);     // pointers to element index 
        double *el_colrowNum = (double*)(el_colrowIndex+nrows+ncols); //and values
        Int j = 0;  //Index inside an element

        for ( Int p = Sp [row]; p < Sp [row+1]; p++){
            el_colrowIndex[j] = Sj[p];
            el_colrowNum[j++] =   Sx[p];
            // adding column tuple
            Tuple colTuple;
            colTuple.e = e;
            colTuple.f = j;
            if ( !paru_add_colTuple (ColList, Sj [p], colTuple, cc) ){
                paru_freemat (&paruMatInfo, cc);
                printf("Out of memory: add_colTuple \n");
                return;
            }
        }
        el_colrowIndex[j++] = 1;  //initializing element row index to 1
    }

    paru_free (anz, sizeof (double), Sx , cc) ;
    paru_free (m, sizeof (Int), W, cc) ;
}
