/** =========================================================================  /
 * =======================  paru_init_rowFronts  ============================  /
 * ========================================================================== */
/*!  @brief  Initializing row fronts; fronts will be assembled later.
 *         Initializing Row and column tuple lists: 
 *            allocating memory and updating lists and initializing matrix 
 *            structre
 *        Assemble each front:
 *            Adding numerical values, allocating data
 *            updating the list
 *          
 * @author Aznaveh
 */
#include "Parallel_LU.hpp"


paru_matrix *paru_init_rowFronts (
        // inputs, not modified
        cholmod_sparse *A,
        int scale, //scales the matrix if > 0
        //symbolic analysis
        paru_symbolic *LUsym,
        // workspace and parameters
        cholmod_common *cc
        )
{

    DEBUGLEVEL(-2);
    if (!A->packed)
    {
        printf ("A is not packed; Wrong format \n");
        return NULL;
    }


    if (LUsym == NULL )
    {
        printf("LUsym is NULL\n");
        return NULL;
    }

    paru_matrix *paruMatInfo;
    paruMatInfo = (paru_matrix*) paru_alloc (1,sizeof(paru_matrix),cc);
    if (paruMatInfo == NULL)
    {   //out of memory
        printf ("Out of memory: paruMatInfo\n");
        return NULL;
    }




    Int m, n, nf;  

    paruMatInfo->LUsym = LUsym;

    m = paruMatInfo->m = LUsym->m - LUsym->n1;   
    n = paruMatInfo->n = LUsym->n - LUsym->n1; 
    nf =  LUsym->nf; 
    paruMatInfo->panel_width = 32;

    Int *row_degree_bound = (Int*) paru_alloc (m, sizeof (Int), cc);
    if (row_degree_bound == NULL)
    {   //out of memory

        printf ("Out of memory: row_degree_bound\n");
        return NULL;
    }

    paruMatInfo->row_degree_bound = row_degree_bound;
 

    Int *rowSize= (Int*) paru_alloc (m, sizeof (Int), cc);
    if (rowSize == NULL)
    {   //out of memory
        printf ("Out of memory: Work\n");
        return NULL;
    }
    memset (rowSize, -1, m*sizeof(Int));
    PRLEVEL (1, ("%% rowSize pointer=%p size=%ld \n", rowSize,m*sizeof(Int) ));


    Int *elRow = (Int*) paru_alloc (m+nf, sizeof (Int), cc);
    if (elRow == NULL)
    {   //out of memory
        printf ("Out of memory: Work\n");
        return NULL;
    }
    memset (elRow, -1, (m+nf)*sizeof(Int));
    PRLEVEL (1, ("%% elRow=%p\n",elRow));

    Int *elCol = (Int*) paru_alloc (m+nf, sizeof (Int), cc);
    if (elCol == NULL)
    {   //out of memory
        printf ("Out of memory: Work\n");
        return NULL;
    }
    memset (elCol, -1, (m+nf)*sizeof(Int));
    PRLEVEL (1, ("%% elCol=%p\n",elCol));



    work_struct *Work= (work_struct*) paru_alloc (1, sizeof (work_struct), cc);
    if (Work == NULL)
    {   //out of memory
        printf ("Out of memory: Work\n");
        return NULL;
    }


    Work->rowSize = rowSize;
    Work->rowMark = 0;


    Work->elRow= elRow;
    Work->elRMark= 0;
    Work->elCol= elCol;
    Work->elCMark= 0;

    PRLEVEL (1, ("%% Work =%p\n ", Work));
    paruMatInfo->Work = Work;

#ifdef COUNT_FLOPS
    //flop count info init
    paruMatInfo->flp_cnt_dgemm = 
    paruMatInfo->flp_cnt_trsm =
    paruMatInfo->flp_cnt_dger = 0.0;
#endif


    PRLEVEL (1, ("%% m=%ld, n=%ld\n",m,n));
    // RowList, ColList and elementList are place holders 
    // pointers to pointers that are allocated
    if (m == 0 || n == 0) 
    {
        printf("%%The dimension of matrix is zero: %ld x %ld \n",m,n);
        paru_free (1, sizeof(paru_matrix), paruMatInfo, cc);
        return NULL;

    }

    tupleList *RowList= paruMatInfo->RowList =
        (tupleList*) paru_alloc (1, m*sizeof(tupleList), cc);
    if (RowList == NULL)
    {   //out of memory
        paru_freemat (&paruMatInfo, cc);
        printf ("Out of memory: RowList\n");
        return NULL;
    }

    PRLEVEL (1, ("%% $RowList =%p\n", RowList));

    tupleList *ColList= paruMatInfo->ColList =
        (tupleList*) paru_alloc (1, n*sizeof(tupleList), cc);
    if (ColList== NULL)
    {   //out of memory
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory: ColList\n");
        return NULL;
    }
    PRLEVEL (1, ("%% $ColList =%p\n", ColList));

    paru_Element **elementList; 
    elementList = paruMatInfo->elementList = // Initialize with NULL
        (paru_Element**) paru_calloc (1, (m+nf+1)*sizeof(paru_Element), cc);
    if (elementList == NULL)
    {   //out of memory
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory: elementList\n");
        return NULL;
    }



    paruMatInfo->frowCount = (Int*) paru_alloc (1, nf*sizeof(Int), cc);
     if(paruMatInfo->frowCount == NULL)
     {
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory frowCount\n");
        return NULL;
    }
    paruMatInfo->fcolCount = (Int*) paru_alloc (1, nf*sizeof(Int), cc);
    if(paruMatInfo->fcolCount == NULL)
    {
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory fcolCount\n");
        return NULL;
    }
 
    paruMatInfo->frowList = (Int**) paru_calloc (1, nf*sizeof(Int*), cc);
     if(paruMatInfo->frowList == NULL)
     {
        paru_freemat (&paruMatInfo, cc);
        printf("Out of memory frowList \n");
        return NULL;
     }

     paruMatInfo->fcolList = (Int**) paru_calloc (1, nf*sizeof(Int*), cc);
     if(paruMatInfo->fcolList == NULL)
     {
         paru_freemat (&paruMatInfo, cc);
         printf("Out of memory fcolList \n");
         return NULL;
     }


     paruMatInfo->partial_Us = //Initialize with NULL
         (paru_fac*) paru_calloc (1, nf*sizeof(paru_fac), cc);
     if (paruMatInfo->partial_Us == NULL)
     {   //out of memory
         paru_freemat (&paruMatInfo, cc);
         printf("Out of memory: Us\n");
         return NULL;
     }

     paruMatInfo->partial_LUs= //Initialize with NULL
         (paru_fac*) paru_calloc (1, nf*sizeof(paru_fac), cc);
     if (paruMatInfo->partial_LUs == NULL)
     {   
         printf("Out of memory: LUs\n");
         return NULL;
     }


     paruMatInfo->time_stamp = 
         (Int*) paru_alloc (1, nf*sizeof(Int), cc);
     if (paruMatInfo->time_stamp == NULL)
     {   //out of memory
         paru_freemat (&paruMatInfo, cc);
         printf("Out of memory: time_stamp\n");
         return NULL;
     }

     paruMatInfo->heapList =  (std::vector<Int>**) 
         paru_calloc (1, (m+nf+1)*sizeof(std::vector<Int>*), cc);
     if (paruMatInfo->heapList == NULL)
     {   //out of memory
         paru_freemat (&paruMatInfo, cc);
         printf("Out of memory: heapList\n");
         return NULL;
     }



   //~~~~~~~~~~ scaling the A matrix
   //TODO change this part using S matrix
   if (scale)
   {
       Int *Ap = (Int*) A->p;
       Int *Ai = (Int*) A->i;
       double *Ax = (double*) A->x;
       double *max_row = (double*) paru_calloc (m, sizeof (double), cc);
       if (max_row == NULL)
       {   //out of memory
           paru_freemat (&paruMatInfo, cc);
           printf("of memory: max_row\n");
           return NULL;
       }

      for(Int col = 0; col < n ; col++)
      {  // finding the max in a row
           for ( Int p = Ap [col]; p < Ap [col+1]; p++)
           {
               if ( fabs(Ax[p]) > max_row[Ai[p]] )
                   max_row[Ai[p]] = fabs(Ax[p]);
           }
       }
       for(Int col = 0; col < n ; col++)
       { // dividing by the max of row 
           for ( Int p = Ap [col]; p < Ap [col+1]; p++)
           {
               ASSERT (max_row[Ai[p]] > 0 ); //TODO 0 or epsilon
               Ax[p] /= max_row[Ai[p]];
           }
       }


       paruMatInfo->scale_row = max_row;
   }
   else
       paruMatInfo->scale_row = NULL;




     /// ------------------------------------------------------------------------
     // create S = A (p,q)', or S=A(p,q) if S is considered to be in row-form
     // -------------------------------------------------------------------------
     Int *Qfill = LUsym->Qfill;
     // Int *Pinv = LUsym->Pinv;
     Int snz = LUsym->snz; 
     double *Sx = LUsym->Sx;
     Int *Sp = LUsym->Sp;
     Int *Sj= LUsym->Sj;

     Int *Ap = (Int*) A->p;
#ifndef NDEBUG
     Int p = 1;
     PRLEVEL (p, ("\n%% Insid init row fronts\n"));
     PRLEVEL (p, ("%% Sp =\n%%"));
     for (Int i = 0; i <= m; i++)
         PRLEVEL (p, ("%ld ", Sp[i]));
     PRLEVEL (p, ("\n"));

     PRLEVEL (p, ("Sj =\n"));
     for (Int k = 0; k < snz; k++)
         PRLEVEL (p, ("%ld ", Sj[k]));
     PRLEVEL (p, ("\n"));

#endif


     //constants for initialzing lists
     Int slackRow = 2;
     Int slackCol = 2;

     /* allocating column list must happend beforehand
      *  I need number of column members of S
      *  while S is row based and column permutation did nt change 
      *  I am using A so 
      */
     for (Int col = 0 ; col < n ; col++)
     {     //allocating Column list tuple 
         Int j = Qfill ? Qfill [col] : col ;  // col of S is column j of A
         Int ncols = Ap[j+1]-Ap[j];

         PRLEVEL (2, ("%% ncols[%ld]=%ld\n",col,ncols));

         ColList[col].numTuple = 0;
         ColList[col].len = slackCol*ncols;
         ColList[col].list = 
             (paru_Tuple*) paru_alloc (slackCol*ncols, sizeof(paru_Tuple), cc);
         if (ColList[col].list == NULL)
         {   //out of memory
             paru_freemat (&paruMatInfo, cc);
             printf("Out of memory: ColList[col].list\n");
             return NULL;
         }
     }

     PRLEVEL (0, ("InMatrix=[\n") ); //MATLAB matrix, 
     //Activating comments after this parts will break the matlab input matrix



     // allocating row tuples, elements and updating column tuples
     for(Int row = 0; row < m ; row++)
     {  
         Int e = LUsym->row2atree[row]; 
         Int nrows = 1, 
             ncols = Sp[row+1]-Sp[row]; //nrows and ncols of current front/row

         PRLEVEL (1, ("%% element %ld = %ld x %ld\n", e, nrows, ncols));

         row_degree_bound [row] = ncols; //Initialzing row degree

         paru_Element *curEl = elementList[e] = 
             paru_create_element (nrows, ncols,0 ,cc);
         if (curEl == NULL)
         {   //out of memory
             paru_freemat (&paruMatInfo, cc);
             printf("Out of memory: curEl\n");
             return NULL;
         }
         
         std::vector<Int>* curHeap;
         curHeap = paruMatInfo->heapList[e] = new std::vector<Int>;
         curHeap->push_back(e);

//         curEl->next = -1;
//         curEl->prev = -1;
//         curEl->tail= e;

#ifndef NDEBUG  // Printing the pointers info
         Int p=1;
         PRLEVEL (p, ("%% curEl = %p ", curEl));
         Int size= sizeof(paru_Element)+
             sizeof(Int)*(2*(nrows+ncols))+
             sizeof(double)*nrows*ncols;
         PRLEVEL (p, ("size= %ld", size));
         PRLEVEL (p, ("\n"));
#endif

         // Allocating Rowlist and updating its tuples
         RowList[row].list =
             (paru_Tuple*) paru_alloc (slackRow*nrows, sizeof(paru_Tuple), cc);
         if (RowList[row].list == NULL)
         {   //out of memory
             paru_freemat (&paruMatInfo, cc);
             printf("Out of memory: RowList[row].list \n");
             return NULL;
         }
         RowList[row].numTuple = 0;
         RowList[row].len = slackRow;

         paru_Tuple rowTuple;
         rowTuple.e = e;
         rowTuple.f = 0;
         if (paru_add_rowTuple (RowList, row, rowTuple, cc) )
         {
             paru_freemat (&paruMatInfo, cc);
             printf("Out of memory: add_rowTuple \n");
             return NULL;
         }

         //Allocating elements, and updating column tuple list
         Int *el_colrowIndex = colIndex_pointer (curEl); 
         double *el_colrowNum = numeric_pointer (curEl);

         PRLEVEL (1, ("el_colrowIndex =%p, el_colrowNum = %p \n", 
                     el_colrowIndex, el_colrowNum));

         Int j = 0;  //Index inside an element
         for ( Int p = Sp [row]; p < Sp [row+1]; p++)
         {
             // adding column tuple
             paru_Tuple colTuple;
             colTuple.e = e;
             colTuple.f = j;

             el_colrowIndex[j] = Sj[p];
             el_colrowNum[j++] =   Sx[p];
             PRLEVEL (1, ("Sj[%ld] =%ld Sx[%ld]=%lf\n", p, Sj[p], p, Sx[p] ));
             //for Matlab
             PRLEVEL (0, ("%ld,%ld, %.16lf;\n", row+1,Sj[p]+1, Sx[p]) );


             if (paru_add_colTuple (ColList, Sj [p], colTuple, cc) )
             {
                 paru_freemat (&paruMatInfo, cc);
                 printf("Out of memory: add_colTuple \n");
                 return NULL;
             }
         }
         el_colrowIndex[j++] = row;  //initializing element row index 
     }

     PRLEVEL (0, ("];\n") );
     PRLEVEL (0, ("I = InMatrix(:,1);\n") );
     PRLEVEL (0, ("J = InMatrix(:,2);\n") );
     PRLEVEL (0, ("X = InMatrix(:,3);\n") );
     PRLEVEL (0, ("S = sparse(I,J,X);\n") );

     // Free here or if not wil be freed in paru_mem anyway
     paru_free (snz, sizeof (double), Sx , cc) ;
     paru_free (snz, sizeof (Int), Sj, cc);
     LUsym->Sx = NULL; LUsym->Sj = NULL;
     


     return paruMatInfo;
}
