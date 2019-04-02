/** =========================================================================  /
 * =======================  paru_factorize  =================================  /
 * ========================================================================== 
 * @brief Doing the BLAS factorization in different panels and call degree
 * update it it is necessary
 *   
 *         
 * @author Aznaveh
 **/
 

#include "Parallel_LU.hpp"
#define TOLER .1  //pivot tolerance

extern "C" void dgetrf_( BLAS_INT *dim1, BLAS_INT *dim2, double *a, 
        BLAS_INT *lda, BLAS_INT *ipiv, BLAS_INT *info);
//dger is already defined in ~/SuiteSparse/CHOLMOD/Include/cholmod_blas.h
void inline swap_int(Int *a, Int *b){
    Int tmp = *a; *a = *b; *b = tmp;
}

template <class T> void inline swap (T &a,T &b){
    T c(a); a=b; b=c;
}

void swap_rows(double *F, Int *fsRowList, Int m, Int n, Int r1, Int r2){
    //This function also swap rows r1 and r2 wholly and indices 
    if(r1 == r2) return;
    swap(fsRowList[r1], fsRowList[r2]);
    for (Int j=0; j < n; j++){   //each column
        swap(F[j*m+r1],F[j*m+r2]);
    }
}

Int paru_panel_factorize (double *F, Int *fsRowList, Int m, Int n, 
        const Int panel_width, Int panel_num, Int row_end, 
        paru_matrix *paruMatInfo) {
    // works like dgetf2f.f in netlib v3.0  here is a link:
    // https://github.com/xianyi/OpenBLAS/blob/develop/reference/dgetf2f.f
    DEBUGLEVEL(0);
    PRLEVEL (1, ("%% Inside panel factorization %ld \n",panel_num));


    Int *row_degree_bound = paruMatInfo->row_degree_bound;
    Int j1 = panel_num*panel_width; // panel starting column


    //  j1 <= panel columns < j2
    //     last panel might be smaller
    Int j2 = (j1 + panel_width < n ) ? 
        j1 + panel_width : n ;

    Int num_col_panel=  j2 - j1 ;
    PRLEVEL (1, ("%% j1= %ld j2 =%ld \n", j1, j2));
    PRLEVEL (1, ("%% row_end= %ld\n", row_end));


    ASSERT ( row_end >= j2);

#ifndef NDEBUG  // Printing the panel
    Int p = 1;
    PRLEVEL (p, ("%% Starting the factorization\n"));
    PRLEVEL (p, ("%% This Panel:\n"));
    for (Int r = j1; r < row_end; r++){
        PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
        for (Int c = j1; c < j2; c++){
            PRLEVEL (p, (" %2.5lf\t", F[c*m+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
 
    //column jth of the panel
    for (Int j = j1; j < j2 ; j++){ 

        PRLEVEL (1, ("%% j = %ld\n", j));

        //Initializing maximum element in the column
        Int row_max = j;
        Int row_deg_max = row_degree_bound[row_max];
        double maxval = F[ j*m + row_max ];
        PRLEVEL (1, ("%% before search max value= %2.4lf\n", maxval));

        //find max
        for (Int i = j+1 ; i < row_end; i++){ 
            PRLEVEL (1, ("%%i=%ld value= %2.4lf\n", i, F[j*m+i]));
            if (fabs (maxval) < fabs(F[j*m+i])){
                row_max = i; 
                row_deg_max = row_degree_bound[i];
                maxval = F[j*m + i];
            }
        }
        PRLEVEL (1, ("%% max value= %2.4lf\n", maxval));

        //initialzing pivot as max numeric value
        Int row_sp = row_max;
        Int row_deg_sp= row_degree_bound[row_max];
        double piv= maxval;


        //find sparsest between accepteble ones
        for (Int i = j; i < row_end; i++) 
            if ( fabs(TOLER*maxval) < fabs(F[j*m+i]) &&//numerically acceptalbe
                    row_degree_bound[i] < row_deg_sp){     // and sparser
                piv = F[j*m+i];
                row_deg_sp = row_degree_bound[i];
                row_sp = i;
            }

        PRLEVEL (1, ("%% piv value= %2.4lf\n", piv));
        //swap rows
        PRLEVEL (1, ("%% Swaping rows j=%ld, spr=%ld\n", j, row_sp));
        swap_rows (F, fsRowList, m , n, j, row_sp);

#ifndef NDEBUG  // Printing the pivotal front
        p = 1;
        if (row_sp != row_max)
            PRLEVEL (p, ("%% \n"));
        PRLEVEL (p, ("%% After Swaping\n"));
        PRLEVEL (p, (" ;\n"));
        for (Int r = 0; r < row_end; r++){
            PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
            for (Int c = 0; c < num_col_panel; c++){
                PRLEVEL (p, (" %2.5lf\t", F[c*m+ r]));
            }
            PRLEVEL (p, ("\n"));
        }
#endif


        //dscal //TODO?: loop unroll is also possible

        if ( j < row_end-1 && piv !=0 ){

            PRLEVEL (1, ("%% dscal\n"));
            for (Int i = j +1 ; i < row_end; i++){ 
                PRLEVEL (1, ("%%i=%ld before cal value= %2.4lf", i, F[j*m+i]));
                F[j*m + i] /= piv;
                PRLEVEL (1, (" after dscal value= %2.4lf\n", i, F[j*m+i]));
            }
        }


        //dger
        /*               dgemm   A := alpha *x*y**T + A
         *
         *
         *                <----------fp------------------->
         *                        j1 current j2    
         *                         ^  panel  ^
         *                         |         |
         *                         |   j     |
         *             F           |   ^     |
         *              \  ____..._|___|_____|__________...___
         * ^              |\      |         |                 |
         * |              | \     |<--panel | rest of         |
         * |              |  \    | width-> |  piv front      |
         * |              |___\...|_______ _|_________ ... ___|
         * |   ^    j1--> |       |\* *|    |                 |
         * | panel        |       |**\*|    |                 |
         * | width        |       |___|Pyyyy|                 |
         * |   v          |____...|___|xAAAA__________...____ |
         * |        j2--> |       |   |xAAAA|                 |
         * rowCount       |       |   |xAAAA|                 |
         * |              |       |   |xAAAA|                 |
         * |              |       |   |xAAAA|                 |
         * |              |       |row_end  |                 |
         * |              .       .         .                 . 
         * |              .       .         .                 . 
         * |              .       .         .                 . 
         * v              |___....____________________..._____|              
         *                         
         */


        if ( j < j2 - 1 && maxval != 0 ){
            BLAS_INT M = (BLAS_INT) row_end - 1 - j ; 
            BLAS_INT N = (BLAS_INT) j2 - 1 - j;
            double alpha = -1.0;
            double *X = F + j*m+j+ 1;
            BLAS_INT Incx = (BLAS_INT) 1;
            double *Y = F + j*m+j + m;
            BLAS_INT Incy = (BLAS_INT) m;
            double *A =  F + j*m+j+ m + 1;
            BLAS_INT lda = (BLAS_INT) m;

#ifndef NDEBUG  // Printing dger input
            Int p = 1;
            PRLEVEL (p, ("%% M =%d ",  M));
            PRLEVEL (p, ("N =%d \n %% x= ",  N));
            for (Int i=0; i<M; i++)
                PRLEVEL (p, (" %lf ",  X[i] ));
            PRLEVEL (p, ("\n %% y= ",  N));
            for (Int j=0; j<N; j++)
                PRLEVEL (p, (" %lf ",  Y[j*m] ));
            PRLEVEL (p, ("\n"));

#endif
            BLAS_DGER(&M, &N, &alpha, X ,  &Incx, Y, &Incy , A, &lda);
        }

#ifndef NDEBUG  // Printing the pivotal front
        Int p = 2;
        PRLEVEL (p, ("%% After dger\n"));
        for (Int r = j1; r < row_end; r++){
            PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
            for (Int c = j1; c < j2; c++){
                PRLEVEL (p, (" %2.5lf\t", F[c*m+ r]));
            }
            PRLEVEL (p, ("\n"));
        }
#endif

    }
    return 1;
}

Int paru_dgetrf (double *F, Int *fsRowList, Int lm, Int ln,
        BLAS_INT *ipiv){
    DEBUGLEVEL(0);

    BLAS_INT m = (BLAS_INT) lm;
    BLAS_INT n = (BLAS_INT) ln;


    PRLEVEL (1, (" %d x %d\n", m, n));
#ifndef NDEBUG  // Printing the pivotal front before computation
    Int p = 1;
    PRLEVEL (1, ("Befor factorization:\n"));
    for (Int r = 0; r < m; r++){
        for (Int c = 0; c < n; c++){
            PRLEVEL (p, (" %3.4lf\t", F[c*m+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
    BLAS_INT info;
    BLAS_INT lda = m;

    PRLEVEL (1, ("ipiv =%p\n", ipiv));
    PRLEVEL (1, ("F=%p\n", F));


#ifndef NDEBUG  // Initializing permutation; just for debug
    for (int i = 0; i < lda ; i++){
        ipiv [i] = -1;
    }
    if (m < n ){
        PRLEVEL (0, ("%%!!!!! FAIL m= %d  n= %d\n", m, n));
        return -1;
    }
#endif

#ifndef NDEBUG  // Printing the list of rows
    p = 1;
    PRLEVEL (p, ("Befor factorization (inside factorize): \n"));
    for (int i = 0; i < m ; i++){
        PRLEVEL (p, ("fsRowList [%d] =%d\n",i, fsRowList [i]));
    }
    PRLEVEL (p, ("\n"));
#endif



    dgetrf_(&m, &n, F, &lda, ipiv, &info);


    ASSERT (m >= n);

    /* changing swap permutation to a real permutation */


#ifndef NDEBUG  // Printing the swap permutation
    p = 1;
    // ATTENTION: ipiv is 1 based
    PRLEVEL (p, ("swap permutation:\n"));
    for (int i = 0; i < m; i++){
        PRLEVEL (p, ("ipiv[%d] =%d\n",i, ipiv[i]));
    }
    PRLEVEL (p, ("\n"));
#endif 

    PRLEVEL (1, (" m=%d n=%d\n", m, n));


    // swap (fsRowList[ipiv [i]], fsRowList[i] ) and it is off by one
    for (Int i = 0; i < n; i++){
        PRLEVEL (1, ("ipiv[%d] =%d\n", i, ipiv[i]));
        ASSERT (ipiv [i] > 0);
        ASSERT (ipiv [i] <= m);
        Int tmp =  fsRowList[ipiv [i]-1];
        PRLEVEL (1, ("tmp =%ld\n", tmp));
        ASSERT (tmp >= 0);

        fsRowList[ipiv [i]-1] = fsRowList [i];
        fsRowList [i] = tmp;
    }

#ifndef NDEBUG   // Printing the LU decomposition
    p = 1;
    PRLEVEL (p, ("After factorization:\n"));
    for (Int r = 0; r < m; r++){
        for (Int c = 0; c < n; c++){
            PRLEVEL (p, (" %3.1lf\t", F[c*m+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif

    PRLEVEL (1, ("info = %ld\n", info));
    if (info != 0 ){
        printf("%%Some problem in factorization\n");
        return info;
    }
    return 0;
}

Int paru_factorize(double *F, Int *fsRowList, Int rowCount, Int f, 
        Int *panel_row, paru_matrix *paruMatInfo){
    DEBUGLEVEL (0);

    Int *Super = paruMatInfo->LUsym->Super;
    Int col1 = Super [f];       /* fornt F has columns col1:col2-1 */
    Int col2 = Super [f+1];
    Int fp = col2 - col1;   /* first fp columns are pivotal */ 


    work_struct *Work =  paruMatInfo->Work;
    Int *row_degree_bound = paruMatInfo->row_degree_bound;
    BLAS_INT *ipiv = (BLAS_INT*) (Work->scratch+rowCount);
    //    return paru_panel_factorize ( F, fsRowList, rowCount, fp, 
    //                fp, 0, rowCount, paruMatInfo);
    //return paru_dgetrf (F , fsRowList, rowCount, fp, ipiv);

    Int panel_width = paruMatInfo->panel_width;
    for(Int panel_num = 0; ; panel_num++){

#ifndef NDEBUG  // Printing the pivotal front
    Int p = 1;
    PRLEVEL (p, ("%%Pivotal Front Before %ld\n",panel_num));
    
    for (Int r = 0; r < rowCount; r++){
        PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
        for (Int c = 0; c < fp; c++){
            PRLEVEL (p, (" %2.5lf\t", F[c*rowCount+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
 
        Int row_end = panel_row [panel_num];
        Int j1 = panel_num*panel_width;
        Int j2 = (panel_num+1)*panel_width;
        // factorize current panel
        paru_panel_factorize ( F, fsRowList, rowCount, fp, 
                panel_width, panel_num, row_end, paruMatInfo);
        if ( j2 >= fp) //if it is not the last panel
            break;
        /*               trsm
         *
         *        F = fully summed part of the pivotal front
         *           op( A ) * B = alpha*B 
         *
         *                <----------fp------------------->
         *                        j1 current j2    
         *    F                    ^  panel  ^
         *     \           ____..._|_________|__________...___
         * ^              |\      |         |                 |
         * |              | \     |<--panel | rest of         |
         * |              |  \    | width-> |  piv front      |
         * |              |___\...|_______ _|_________ ... ___|
         * |   ^    j1--> |       |\        |                 |
         * | panel        |       |**\ A    |   B(In out)     |
         * | width        |       |*L**\    |                 |
         * |   v          |____...|******\ _|_________...____ |
         * |        j2--> |       |         |                 |
         * rowCount       |       |         |                 |
         * |              .       .         .                 . 
         * |              .       .         .                 . 
         * |              .       .         .                 . 
         * v              |___....____________________..._____|              
         *                         
         */
        ASSERT (j2 < fp);

        {
            BLAS_INT M = (BLAS_INT) panel_width;
            BLAS_INT N = (BLAS_INT) fp - j2;
            double alpha = 1.0;
            double *A = F + j1*rowCount + j1;
            BLAS_INT lda = (BLAS_INT) rowCount;
            double *B = F + j2*rowCount + j1;
            BLAS_INT ldb = (BLAS_INT) rowCount;
#ifndef NDEBUG  
            Int p = 1;
            PRLEVEL (p, ("%% M =%d N = %d alpha = %f \n", M, N, alpha));
            PRLEVEL (p, ("%% lda =%d ldb =%d\n", lda, ldb));
            PRLEVEL (p, ("%% Pivotal Front Before Trsm: %ld x %ld\n",
                        fp, rowCount));
            for (Int r = 0; r < rowCount; r++){
                PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
                for (Int c = 0; c < fp; c++){
                    PRLEVEL (p, (" %2.5lf\t", F[c*rowCount+ r]));
                }
                PRLEVEL (p, ("\n"));
            }

#endif
            BLAS_DTRSM ("L" ,"L" ,"N" ,"U", &M, &N, &alpha, A, &lda, 
                    B, &ldb);


#ifndef NDEBUG  
            PRLEVEL (p, ("%% Pivotal Front After Trsm: %ld x %ld\n %%", 
                        fp, rowCount));
            for (Int r = 0; r < rowCount; r++){
                PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
                for (Int c = 0; c < fp; c++){
                    PRLEVEL (p, (" %2.5lf\t", F[c*rowCount+ r]));
                }
                PRLEVEL (p, ("\n"));
            }
#endif
        }




        /*               dgemm   C := alpha*op(A)*op(B) + beta*C
         *
         *        F = fully summed part of the pivotal front
         *
         *                <----------fp------------------->
         *                        j1 current j2    
         *    F                    ^  panel  ^
         *     \           ____..._|_________|__________...___
         * ^              |\      |         |                 |
         * |              | \     |<--panel | rest of         |
         * |              |  \    | width-> |  piv front      |
         * |              |___\...|_______ _|_________ ... ___|
         * |   ^    j1--> |       |\        |**************** |
         * | panel        |       |  \      |****In***B****** |
         * | width        |       |    \    |**************** |
         * |   v          |____...|_______\___________...____ |
         * |        j2--> |       |******** |ccccccccccccccccc|
         * rowCount       |       |******** |ccccccccccccccccc|
         * |              .       .******** .ccccccCcccccccccc. 
         * |              .       .***In*** .ccccOutcccccccccc. 
         * |              .       .***A**** .ccccccccccccccccc. 
         * |              |       |******** |ccccccccccccccccc|              
         * |              |       |******** |ccccccccccccccccc|              
         * |              |       |row_end  |                 |              
         * |              |       |         |                 |              
         * v              |___....|_______ _|__________...____|              
         *                         
         */


        {

            BLAS_INT M = (BLAS_INT) (row_end - j2); 

            BLAS_INT N = (BLAS_INT) fp - j2;
            BLAS_INT K = (BLAS_INT) panel_width;
            double alpha = -1;
            double *A = F+ j1*rowCount + j2;
            BLAS_INT lda = (BLAS_INT) rowCount;
            double *B = F+ j2*rowCount + j1;
            BLAS_INT ldb = (BLAS_INT) rowCount;
            double beta = 1; //keep current values
            double *C = F+ j2*rowCount + j2;
            BLAS_INT ldc = (BLAS_INT) rowCount ;
#ifndef NDEBUG  
            Int p = 1;
            PRLEVEL (p, ("%% DGEMM "));
            PRLEVEL (p, ("%% M =%d K = %d N = %d alpha = %f \n",
                        M, K, N,  alpha));
            PRLEVEL (p, ("%% lda =%d ldb =%d\n", lda, ldb));
            PRLEVEL (p, ("%% j2 =%ld j1=%ld\n", j2, j1));
            PRLEVEL (p, ("\n %%"));
#endif


            BLAS_DGEMM ("N", "N", &M, &N, &K, &alpha, A, &lda, B, &ldb, 
                    &beta, C, &ldc);


        }

#ifndef NDEBUG  
        PRLEVEL (p, ("%% Pivotal Front After Dgemm: %ld x %ld\n %%", 
                    fp, rowCount));
        for (Int r = 0; r < rowCount; r++){
            PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
            for (Int c = 0; c < fp; c++){
                PRLEVEL (p, (" %2.5lf\t", F[c*rowCount+ r]));
            }
            PRLEVEL (p, ("\n"));
        }

#endif


        paruMatInfo->time_stamp[f]++;
        //        paru_update_rowDeg ( panel_num, row_end, colCount, f, 
        //                CBColList, paruMatInfo);
        paruMatInfo->time_stamp[f]+= 2;


    }


    return 0;
}
