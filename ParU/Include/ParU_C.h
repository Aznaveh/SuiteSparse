// ============================================================================/
// ======================= ParU_C.h ===========================================/
// ============================================================================/


//------------------------------------------------------------------------------
#ifndef PARU_C_H
#define PARU_C_H

#ifdef __cplusplus
extern "C" {
#endif 



// =============================================================================
// ========================= ParU_C_Control ====================================
// =============================================================================
// Just like ParU_Control the only difference is the initialization whic is 
// handled in the C interface
typedef struct ParU_C_Control_struct
{
    Int mem_chunk;  // chunk size for memset and memcpy

    // Symbolic controls
    Int umfpack_ordering;
    Int umfpack_strategy;  // symmetric or unsymmetric
    Int umfpack_default_singleton; //filter singletons if true

    Int relaxed_amalgamation_threshold;
    // symbolic analysis tries that each front have more pivot columns
    // than this threshold

    // Numeric controls
    Int scale;         // if 1 matrix will be scaled using max_row
    Int panel_width;  // width of panel for dense factorizaiton
    Int paru_strategy;  // the same strategy umfpack used

    double piv_toler;     // tolerance for accepting sparse pivots
    double diag_toler;  // tolerance for accepting symmetric pivots
    Int trivial;  // dgemms with sizes less than trivial doesn't call BLAS
    Int worthwhile_dgemm;  // dgemms bigger than worthwhile are tasked
    Int worthwhile_trsm;  // trsm bigger than worthwhile are tasked
    Int paru_max_threads;    // It will be initialized with omp_max_threads
    // if the user do not provide a smaller number
} ParU_C_Control;

//  // =========================================================================
//  // ========================= ParU_C_Symbolic ===============================
//  // =========================================================================
//  // just a carrier for the C++ data structure
//  typedef struct ParU_C_Symbolic_struct
//  {
//      void *sym_handle;
//  } ParU_C_Symbolic;
//  
//  // =========================================================================
//  // ========================= ParU_C_Numeric ================================
//  // =========================================================================
//  // just a carrier for the C++ data structure
//  typedef struct ParU_C_Numeric_struct
//  {
//      void *num_handle;
//  } ParU_C_Numeric;
#define ParU_C_Symbolic void
#define ParU_C_Numeric void


// ParU_Version: 
// print out the version
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Version (int ver [3], char date [128]);
//------------------------------------------------------------------------------
// ParU_Analyze: Symbolic analysis is done in this routine. UMFPACK is called
// here and after that some more speciallized symbolic computation is done for
// ParU. ParU_Analyze can be called once and can be used for differen 
// ParU_Factorize calls. 
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Analyze(
        // input:
        cholmod_sparse *A,  // input matrix to analyze ...
        // output:
        ParU_C_Symbolic **Sym_handle,  // output, symbolic analysis
        // control:
        ParU_C_Control *Control);

//------------------------------------------------------------------------------
// ParU_Factorize: Numeric factorization is done in this routine. Scaling and
// making Sx matrix, computing factors and permutations is here. ParU_C_Symbolic
// structure is computed ParU_Analyze and is an input in this routine.
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Factorize(
        // input:
        cholmod_sparse *A, ParU_C_Symbolic *Sym,
        // output:
        ParU_C_Numeric **Num_handle,
        // control:
    ParU_C_Control *Control);

//------------------------------------------------------------------------------
//--------------------- Solve routines -----------------------------------------
//------------------------------------------------------------------------------
// In all the solve routines Num structure must come with the same Sym struct
// that comes from ParU_Factorize
//-------- Ax = b (x is overwritten on b)---------------------------------------
ParU_Ret ParU_C_Solve_Axx(
    // input:
    ParU_C_Symbolic *Sym, ParU_C_Numeric *Num,
    // input/output:
    double *b,
    // control:
    ParU_C_Control *Control);
//-------- Ax = b --------------------------------------------------------------
ParU_Ret ParU_C_Solve_Axb(
    // input:
    ParU_C_Symbolic *Sym, ParU_C_Numeric *Num, double *b,
    // output
    double *x,
    // control:
    ParU_C_Control *user_Control);
//-------- AX = B  (X is overwritten on B, multiple rhs)------------------------
ParU_Ret ParU_C_Solve_AXX(
    // input
    ParU_C_Symbolic *Sym, ParU_C_Numeric *Num, Int nrhs,
    // input/output:
    double *B,  // m(num_rows of A) x nrhs
    // control:
    ParU_C_Control *Control);
//-------- AX = B  (multiple rhs)-----------------------------------------------
ParU_Ret ParU_C_Solve_AXB(
    // input
    ParU_C_Symbolic *Sym, ParU_C_Numeric *Num, Int nrhs, double *B,
    // output:
    double *X,
    // control:
    ParU_C_Control *Control);

//------------------------------------------------------------------------------
//-------------- computing residual --------------------------------------------
//------------------------------------------------------------------------------
// The user provide both x and b
// resid = norm1(b-A*x) / norm1(A)
ParU_Ret ParU_C_Residual_bAx(
    // inputs:
    cholmod_sparse *A, double *x, double *b, Int m,
    // output:
    double &resid, double &anorm,
    // control:
    ParU_C_Control *Control);

// resid = norm1(B-A*X) / norm1(A) (multiple rhs)
ParU_Ret ParU_C_Residual_BAX(
    // inputs:
    cholmod_sparse *A, double *X, double *B, Int m, Int nrhs,
    // output:
    double &resid, double &anorm,
    // control:
    ParU_C_Control *Control);

//------------------------------------------------------------------------------
//------------ Free routines----------------------------------------------------
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Freenum(ParU_C_Numeric **Num_handle, ParU_C_Control *Control);

ParU_Ret ParU_C_Freesym(ParU_C_Symbolic **Sym_handle, ParU_C_Control *Control);

#ifdef __cplusplus
}
#endif 

#endif //PARU_C_H
