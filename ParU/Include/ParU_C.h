// ParU_C.h

// C wrappers defined here
// TODO: I am not sure if I should have a defintion file or not
// but it might be a good idea
//
// I am not sure how to initialize ParU_Control and even how to control it in C
// Unsure how to control ParU_Ret
//
//
// Also unsure if I can send Symbolic data structure in C
//
//------------------------------------------------------------------------------
#ifndef PARU_C_H
#define PARU_C_H

extern "C" {
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
        ParU_Symbolic **Sym_handle,  // output, symbolic analysis
        // control:
        ParU_Control *Control);

//------------------------------------------------------------------------------
// ParU_Factorize: Numeric factorization is done in this routine. Scaling and
// making Sx matrix, computing factors and permutations is here. ParU_Symbolic
// structure is computed ParU_Analyze and is an input in this routine.
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Factorize(
    // input:
    cholmod_sparse *A, ParU_Symbolic *Sym,
    // output:
    ParU_Numeric **Num_handle,
    // control:
    ParU_Control *Control);

//------------------------------------------------------------------------------
//--------------------- Solve routines -----------------------------------------
//------------------------------------------------------------------------------
// In all the solve routines Num structure must come with the same Sym struct
// that comes from ParU_Factorize
//-------- Ax = b (x is overwritten on b)---------------------------------------
ParU_Ret ParU_C_Solve_axx(
    // input:
    ParU_Symbolic *Sym, ParU_Numeric *Num,
    // input/output:
    double *b,
    // control:
    ParU_Control *Control);
//-------- Ax = b --------------------------------------------------------------
ParU_Ret ParU_C_Solve_axb(
    // input:
    ParU_Symbolic *Sym, ParU_Numeric *Num, double *b,
    // output
    double *x,
    // control:
    ParU_Control *user_Control);
//-------- AX = B  (X is overwritten on B, multiple rhs)------------------------
ParU_Ret ParU_C_Solve_AXX(
    // input
    ParU_Symbolic *Sym, ParU_Numeric *Num, Int nrhs,
    // input/output:
    double *B,  // m(num_rows of A) x nrhs
    // control:
    ParU_Control *Control);
//-------- AX = B  (multiple rhs)-----------------------------------------------
ParU_Ret ParU_C_Solve_AXB(
    // input
    ParU_Symbolic *Sym, ParU_Numeric *Num, Int nrhs, double *B,
    // output:
    double *X,
    // control:
    ParU_Control *Control);

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
    ParU_Control *Control);

// resid = norm1(B-A*X) / norm1(A) (multiple rhs)
ParU_Ret ParU_C_Residual_BAX(
    // inputs:
    cholmod_sparse *A, double *X, double *B, Int m, Int nrhs,
    // output:
    double &resid, double &anorm,
    // control:
    ParU_Control *Control);

//------------------------------------------------------------------------------
//------------ Free routines----------------------------------------------------
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Freenum(ParU_Numeric **Num_handle, ParU_Control *Control);

ParU_Ret ParU_C_Freesym(ParU_Symbolic **Sym_handle, ParU_Control *Control);
}
#endif
