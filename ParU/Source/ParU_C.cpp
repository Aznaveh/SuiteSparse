////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// ParU_C.cpp ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.
// SPDX-License-Identifier: GNU GPL 3.0
//
/*! @brief  This C++ file provides a set of C-callable wrappers so that a C
 *  program can call ParU.
 *
 *  @author Aznaveh
 */
#include "ParU_C.h"

#include "paru_internal.hpp"

extern "C"
{
// ParU_Version: 
// print out the version
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Version (int ver [3], char date [128]){}
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
        ParU_C_Control *Control)
{ 
    //TODO: copy the inside of the Control and then 
    // make a new symbolic object
    // call the Cpp ParU_Analyze
}

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
    ParU_C_Control *Control)
{ 
    //TODO: copy the inside of the Control and then 
    // make a new numeric object
    // call the Cpp ParU_Factorize
}

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
    ParU_C_Control *user_Control)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}


//-------- AX = B  (X is overwritten on B, multiple rhs)------------------------
ParU_Ret ParU_C_Solve_AXX(
    // input
    ParU_C_Symbolic *Sym, ParU_C_Numeric *Num, Int nrhs,
    // input/output:
    double *B,  // m(num_rows of A) x nrhs
    // control:
    ParU_C_Control *Control)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}


//-------- AX = B  (multiple rhs)-----------------------------------------------
ParU_Ret ParU_C_Solve_AXB(
    // input
    ParU_C_Symbolic *Sym, ParU_C_Numeric *Num, Int nrhs, double *B,
    // output:
    double *X,
    // control:
    ParU_C_Control *Control)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}

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
    ParU_C_Control *Control)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}


// resid = norm1(B-A*X) / norm1(A) (multiple rhs)
ParU_Ret ParU_C_Residual_BAX(
    // inputs:
    cholmod_sparse *A, double *X, double *B, Int m, Int nrhs,
    // output:
    double &resid, double &anorm,
    // control:
    ParU_C_Control *Control)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}

//------------------------------------------------------------------------------
//------------ Free routines----------------------------------------------------
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Freenum(ParU_C_Numeric **Num_handle, ParU_C_Control *Control)
{ 
    //TODO: // No need to copy the inside of the Control and then 
    // free the numeric object
}


ParU_Ret ParU_C_Freesym(ParU_C_Symbolic **Sym_handle, ParU_C_Control *Control)
{ 
    //TODO: // No need to copy the inside of the Control and then 
    // free the symblic object
}


} //extern c
