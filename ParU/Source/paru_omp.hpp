////////////////////////////////////////////////////////////////////////////////
//////////////////////////  paru_omp.hpp ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#ifndef PARU_OMP_H
#define PARU_OMP_H
//!
// definitions of using OpenMP ins ParU 
//  @author Aznaveh
//

#if defined ( _OPENMP )
    #include <omp.h>
    #define PARU_OPENMP_MAX_THREADS       omp_get_max_threads ( )
    #define PARU_OPENMP_GET_NUM_THREADS   omp_get_num_threads ( )
    #define PARU_OPENMP_GET_WTIME         omp_get_wtime ( )
    #define PARU_OPENMP_GET_THREAD_ID     omp_get_thread_num ( )
    #define PARU_OPENMP_SET_DYNAMIC(d)   omp_set_dynamic(d)
    #define PARU_OPENMP_SET_MAX_ACTIVE_LEVELS(l)   omp_set_max_active_levels(l)

#else

    #define PARU_OPENMP_MAX_THREADS       (1)
    #define PARU_OPENMP_GET_NUM_THREADS   (1)
    #define PARU_OPENMP_GET_WTIME         (0)
    #define PARU_OPENMP_GET_THREAD_ID     (0)
    #define PARU_OPENMP_SET_DYNAMIC       (0)
    #define PARU_OPENMP_SET_MAX_ACTIVE_LEVELS(l)   ()
#endif

#endif
