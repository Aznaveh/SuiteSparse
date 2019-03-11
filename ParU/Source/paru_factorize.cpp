/** =========================================================================  /
 * =======================  paru_factorize  =================================  /
 * ========================================================================== 
 * @brief Doing the BLAS factorization 
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
    DEBUGLEVEL(1);
    PRLEVEL (1, ("%% Inside panel factorization \n"));


    Int *row_degree_bound = paruMatInfo->row_degree_bound;
    Int col_st = panel_num*panel_width; // panel starting column

    PRLEVEL (1, ("%% col_st= %ld\n", col_st));
    PRLEVEL (1, ("%% row_end= %ld\n", row_end));


    //  col_st <= panel columns < col_end
    //     last panel might be smaller
    Int col_end = (col_st + panel_width < n ) ? 
        col_st + panel_width : n ;

    Int num_col_panel=  col_end - col_st ;



#ifndef NDEBUG  // Printing the pivotal front
    Int p = 1;
    PRLEVEL (p, ("%% Starting the factorization\n"));
    PRLEVEL (p, (" ;\n"));
    for (Int r = 0; r < row_end; r++){
        PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
        for (Int c = 0; c < num_col_panel; c++){
            PRLEVEL (p, (" %2.5lf\t", F[c*row_end+ r]));
        }
        PRLEVEL (p, ("\n"));
    }
#endif
 
    //column ith of the panel
    for (Int j = col_st; j < col_end ; j++){ 

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
        if (maxval == 0){
            if ( j == col_end -1)
                continue;
            else{
                printf ("Singular submatrix\n");
                return -1;
            }
        }

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
                PRLEVEL (p, (" %2.5lf\t", F[c*row_end+ r]));
            }
            PRLEVEL (p, ("\n"));
        }
#endif


        //dscal //TODO?: loop unroll is also possible

        if ( j < row_end-1){

            PRLEVEL (1, ("%% dscal\n"));
            for (Int i = j +1 ; i < row_end; i++){ 
                PRLEVEL (1, ("%%i=%ld before cal value= %2.4lf", i, F[j*m+i]));
                F[j*m + i]= F[j*m + i]/piv;
                PRLEVEL (1, (" after dscal value= %2.4lf\n", i, F[j*m+i]));
            }
        }


        //dger
        if ( j < col_end - 1){
            BLAS_INT M = (BLAS_INT) row_end - 1 - j ; 
            BLAS_INT N = (BLAS_INT) col_end - 1 - j;
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
        Int p = 1;
        PRLEVEL (p, ("%% After dger\n"));
        PRLEVEL (p, (" ;\n"));
        for (Int r = 0; r < row_end; r++){
            PRLEVEL (p, ("%% %ld\t", fsRowList [r]));
            for (Int c = 0; c < num_col_panel; c++){
                PRLEVEL (p, (" %2.5lf\t", F[c*row_end+ r]));
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

Int paru_factorize(double *F, Int *fsRowList, Int rowCount, Int fp, 
        paru_matrix *paruMatInfo){

    work_struct *Work =  paruMatInfo->Work;
    Int *row_degree_bound = paruMatInfo->row_degree_bound;
    BLAS_INT *ipiv = (BLAS_INT*) (Work->scratch+rowCount);
    return paru_panel_factorize ( F, fsRowList, rowCount, fp, 
                fp, 0, rowCount, paruMatInfo);
//    return paru_dgetrf (F , fsRowList, rowCount, fp, ipiv);
    return 0;
}
