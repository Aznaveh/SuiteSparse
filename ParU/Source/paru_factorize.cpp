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
    for (Int i=0; i < n; i++){
        swap(F[r1*m+i],F[r2*m+i]);
    }
}

Int paru_panel_factorize (double *F, Int *fsRowList, Int m, Int n, 
        const Int panel_width, Int panel_num, Int num_rows, 
        paru_matrix *paruMatInfo) {
    // works like dgetf2f.f in netlib v3.0  here is a link:
    // https://github.com/xianyi/OpenBLAS/blob/develop/reference/dgetf2f.f

    Int *row_degree_bound = paruMatInfo->row_degree_bound;
    Int panel_point = panel_num*panel_width;

    ASSERT (panel_point < m);

    for (Int i=0; i < panel_width; i++){ //column ith of the panel

        Int row_max = panel_point;
        Int row_deg_max = row_degree_bound[row_max];
        double maxval = F[i*m+row_max];

        //find max
        for (Int j = panel_point+1 ; j < num_rows; j++){ 
            if (maxval < F[i*m+j]){
                row_max = j; row_deg_max = row_degree_bound[j];
                maxval = F[i*m+j];
            }
        }

        //initialzing pivot as max numeric value
        Int row_sp = row_max;
        Int row_deg_sp= row_degree_bound[row_max];
        double piv= maxval;


        //find sparsest between accepteble ones
        for (Int j = panel_point; j < num_rows; j++){ 
            if ( TOLER*maxval < F[i*m+j] ) //numerically acceptalbe
                if (row_degree_bound[j] < row_deg_sp){
                    piv = F[i*m+j];
                    row_deg_sp = row_degree_bound[j];
                    row_sp = j;
            }

        }

        //swap rows
        swap_rows (F, fsRowList, m , n, panel_point, row_sp);

        //dscal //TODO?: loop unroll is also possible
        for (Int j = panel_point+1 ; j < num_rows; j++){ 
            F[i*m+j]= F[i*m+j]/piv;
        }

        //dger
        BLAS_INT M = (BLAS_INT) num_rows-1;
        BLAS_INT N = (BLAS_INT) n-panel_point;
        double alpha = -1.0;
        double *X = F+panel_point*panel_point+1;
        BLAS_INT Incx = (BLAS_INT) 1;
        double *Y = F+panel_point*panel_point+m;
        BLAS_INT Incy = (BLAS_INT) m;
        double *A =  F+panel_point*panel_point+m+1;
        BLAS_INT lda = (BLAS_INT) m;

        BLAS_DGER(&M, &N, &alpha, X ,  &Incx, Y, &Incy , A, &lda);
        

    }
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
    return paru_dgetrf (F , fsRowList, rowCount, fp, ipiv);
    return 0;
}
