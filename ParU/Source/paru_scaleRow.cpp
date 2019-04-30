/** =========================================================================  /
 * =======================  paru_init_rowFronts  ============================  /
 * ==========================================================================  /
 * @brief scaling rows of the matrix
 *        finding the max in each row and then dividing all elements by that 
 * @author Aznaveh
 * */

#include "Parallel_LU.hpp"
paru_scaleRow (paruMatInfo){
    

    for(Int row = 0; row < m ; row++){  
        for ( Int p = Sp [row]; p < Sp [row+1]; p++){

            = Sj[p];
            double X = Sx[p];
            PRLEVEL (1, ("Sj[%ld] =%ld Sx[%ld]=%lf\n", p, Sj[p], p, Sx[p] ));
            //for Matlab
            PRLEVEL (0, ("%ld,%ld, %.16lf;\n", row+1,Sj[p]+1, Sx[p]) );


        }
    }
}

