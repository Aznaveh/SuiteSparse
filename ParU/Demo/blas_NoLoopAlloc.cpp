#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>

#include "Parallel_LU.hpp"

int main () {

    // Create a text string, which is used to output the text file
    std::string oneLine;

    // Read from the text file
    std::ifstream InputFile("filename.txt");


    ///// Finding the max
    long max = 0;
    int i = 0;
    // Use a while loop together with the getline() function to read the file line by line
    while (getline (InputFile, oneLine)) {
        std::stringstream ss(oneLine);

        int tmp; double dtmp;
        ss >> tmp; ss>>dtmp;
        int m,n,k,lda,ldb,ldc;
        ss >> m;     max = max > m ? max: m;
        ss >> n;     max = max > n ? max: n;
        ss >> k;     max = max > k ? max: k;
        ss >> lda;    max = max > lda ? max: lda;
        ss >> ldb;    max = max > ldb ? max: ldb;
        ss >> ldc;    max = max > ldc ? max: ldc;
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

    double my_start_time = omp_get_wtime();

    double *C = (double*) malloc (max*max*sizeof(double));
    double *A = (double*) malloc (max*max*sizeof(double));
    double *B = (double*) malloc (max*max*sizeof(double));

    // Use a while loop together with the getline() 
    // function to read the file line by line
    
    while (getline (InputFile, oneLine)) {
        std::stringstream ss(oneLine);
        int f; double gemm_time;
        ss >> f; ss>>gemm_time;

        BLAS_INT m,n,k,lda,ldb,ldc;
        ss >> m;     
        ss >> n;    
        ss >> k;   
        ss >> lda; 
        ss >> ldb;  
        ss >> ldc;   
        // C = A*B
        BLAS_DGEMM ("N" ,"N" , &m, &n, &k, &alpha, A, 
                &lda, B, &ldb, &beta, C, &ldc);
        ++i;
        flops += m*n*k;
        }

    free(A);
    free(B);
    free(C);

    std::cout << "nol=" << i <<" flops = " << flops << std::endl;

    double my_time = omp_get_wtime() - my_start_time;
    std::cout  << my_time << std::endl;

    // Close the file
    InputFile.close();
}
