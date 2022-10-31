// ============================================================================/
// ======================= ParU_definitions.h =================================/
// ============================================================================/
//
#ifndef PARU_DEFINITIONS_H
#define PARU_DEFINITIONS_H


// ============================================================================/
// ======================= ParU version =======================================/
// ============================================================================/

#define PARU_DATE "Apr 15, 2022"
#define PARU_VERSION_MAJOR 0
#define PARU_VERSION_MINOR 1
#define PARU_VERSION_UPDATE 2


#define CHOLMOD_BLAS_H
#ifdef BLAS_INT
#undef BLAS_INT
#endif

#ifdef MKLROOT
// MKL BLAS
#define CHOLMOD_BLAS_H
#include <mkl.h>
#define BLAS_set_num_threads(n) mkl_set_num_threads(n)
#define BLAS_INT MKL_INT
#else
// assume OpenBLAS
#include <cblas.h>
#define BLAS_set_num_threads(n) openblas_set_num_threads(n)
#define BLAS_INT int
#endif

extern "C"
{
#include "cholmod.h"
#include "umfpack.h"
}

#ifdef Int  // defined in amd
#undef Int
#endif
#define Int int64_t

//  the same values as UMFPACK_STRATEGY defined in UMFPACK/Include/umfpack.h
#define PARU_STRATEGY_AUTO 0         // decided to use sym. or unsym. strategy
#define PARU_STRATEGY_UNSYMMETRIC 1  // COLAMD(A), metis, ...
#define PARU_STRATEGY_SYMMETRIC 3    // prefer diagonal

//copied from umfpack.h
// Control [UMFPACK_ORDERING] and Info [UMFPACK_ORDERING_USED] are one of 
#define UMFPACK_ORDERING_CHOLMOD 0      // use CHOLMOD (AMD/COLAMD then METIS)
#define UMFPACK_ORDERING_AMD 1          // use AMD/COLAMD 
//#define UMFPACK_ORDERING_GIVEN 2        // user-provided Qinit 
#define UMFPACK_ORDERING_METIS 3        // use METIS 
#define UMFPACK_ORDERING_BEST 4         // try many orderings, pick best
#define UMFPACK_ORDERING_NONE 5         // natural ordering 
//#define UMFPACK_ORDERING_USER 6         // user-provided function 
// AMD/COLAMD means: use AMD for symmetric strategy, COLAMD for unsymmetric 
//
enum ParU_Ret
{
    PARU_SUCCESS,
    PARU_OUT_OF_MEMORY,
    PARU_INVALID,
    PARU_SINGULAR
};

#endif
