#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>
#include <malloc.h>

#include <memory_resource>

#include "Parallel_LU.hpp"

int main () 
{

   // mallopt (M_MMAP_MAX, 0) ;           // disable mmap; it's too slow
   // mallopt (M_TRIM_THRESHOLD, -1) ;    // disable sbrk trimming
   // mallopt (M_TOP_PAD, 16*1024*1024) ; // increase padding to speedup malloc
    
    // Create a text string, which is used to output the text file
    std::string oneLine;

    // Read from the text file
    std::ifstream InputFile("filename.txt");


    ///// Finding the max
    long max = 0;
    int i = 0;
    // Use a while loop together with the getline()
    // function to read the file line by line
    while (getline (InputFile, oneLine)) {
        std::stringstream ss(oneLine);

        int tmp; double dtmp;
        ss >> tmp; ss>>dtmp;
        int m,n,k,lda,ldb,ldc;
        ss >> m;     max = max > m ? max: m;
        ss >> n;     max = max > n ? max: n;
        ss >> k;     max = max > k ? max: k;
        //ss >> lda;    max = max > lda ? max: lda;
        //ss >> ldb;    max = max > ldb ? max: ldb;
        //ss >> ldc;    max = max > ldc ? max: ldc;
        i++;
    }
    std::cout << "nol=" <<i;
    // Close the file
    InputFile.close();
    std::cout << std::endl << "max matrix size = " << max*max << std::endl;

    double alpha = -1;
    double beta= 0; 
    double flops= 0; 

    // Read from the text file
    InputFile.open("filename.txt");
    i = 0;

    //std::pmr::pool_options::largest_required_pool_block = 40000000;
    // std::pmr::synchronized_pool_resource pool;

    double my_start_time = omp_get_wtime();

    // double *A = new double[max*max];
    // double *B = new double[max*max];

    //double *A = (double*) malloc (max*max*sizeof(double));
    //double *B = (double*) malloc (max*max*sizeof(double));




    // Use a while loop together with the getline() 
    // function to read the file line by line
    while (getline (InputFile, oneLine) ) 
    {
        std::stringstream ss(oneLine);


        int tmp; double dtmp;
        ss >> tmp; ss>>dtmp;

        BLAS_INT m,n,k,lda,ldb,ldc;
        ss >> m;     
        ss >> n;    
        ss >> k;   
        ss >> lda; 
        ss >> ldb;  
        ss >> ldc;   

        if (m*n*k == 0) 
            continue;
        
        double *A = (double*) malloc (lda*k*sizeof(double));
        double *B = (double*) malloc (k*ldc*sizeof(double));
        double *C = (double*) malloc (lda*n*sizeof(double));

        //double *A = (double*) pool.allocate(m*k*sizeof(double),8);
        //double *B = (double*) pool.allocate(k*n*sizeof(double),8);
        //double *C = (double*) pool.allocate(m*n*sizeof(double),8);

       // std::cout<< "m = " << m;
       // std::cout<< " n = " << n;
       // std::cout<< " k = " << k; 
       // std::cout<< " lda = " << lda;
       // std::cout<< " ldb = " << ldb;
       // std::cout<< " ldc = " << ldc << std::endl;;

        //m= 144; n=3; k = 13;
        // C = A*B
        // 
        BLAS_DGEMM ("N" ,"N" , &m, &n, &k, &alpha, A, &lda, B,
                &ldb, &beta, C, &ldc);
        ++i;
        flops += m*n*k;

        //pool.deallocate( (void *) A, m*k*sizeof(double), 8);
        //pool.deallocate( (void *) B, k*n*sizeof(double), 8);
        //pool.deallocate( (void *) C, m*n*sizeof(double), 8);
        
        free(A);
        free(B);
        free(C);
    }
    double my_time = omp_get_wtime() - my_start_time;

    std::cout << "nol=" << i <<" flops = " << flops << std::endl;
    std::cout  << my_time << std::endl;

    // Close the file
    InputFile.close();
}
