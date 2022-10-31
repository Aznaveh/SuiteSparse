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
ParU_Ret ParU_C_Version (int ver [3], char date [128])
    {return ParU_Version (ver ,date);}
// copy the inside of the C structrue to the Cpp structure
void  cp_control (ParU_Control *Control, ParU_C_Control *Control_C)
{
    Control->mem_chunk = Control_C->mem_chunk;

    Control->umfpack_ordering = Control_C->umfpack_ordering;
    Control->umfpack_strategy = Control_C->umfpack_strategy;
    Control->umfpack_default_singleton = Control_C->umfpack_default_singleton;

    Control->relaxed_amalgamation_threshold = 
        Control_C->relaxed_amalgamation_threshold;

    Control->scale = Control_C->scale;
    Control->panel_width = Control_C->panel_width;
    Control->paru_strategy = Control_C->paru_strategy;


    Control->piv_toler = Control_C->piv_toler;
    Control->diag_toler = Control_C->diag_toler;
    Control->trivial = Control_C->trivial;
    Control->worthwhile_dgemm = Control_C->worthwhile_dgemm;
    Control->worthwhile_trsm = Control_C->worthwhile_trsm;
    Control->paru_max_threads = Control_C->paru_max_threads;
}
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
        ParU_C_Symbolic **Sym_handle_C,  // output, symbolic analysis
        // control:
        ParU_C_Control *Control_C)
{ 
    //TODO: copy the inside of the Control and then 
    ParU_Control Control;
    //TODO:  make a new symbolic object

    ParU_Ret info;
    return info = ParU_Analyze(A, &Sym, &Control);
}

//------------------------------------------------------------------------------
// ParU_Factorize: Numeric factorization is done in this routine. Scaling and
// making Sx matrix, computing factors and permutations is here. ParU_C_Symbolic
// structure is computed ParU_Analyze and is an input in this routine.
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Factorize (
        // input:
        cholmod_sparse *A, ParU_C_Symbolic *Sym_C,
        // output:
        ParU_C_Numeric **Num_handle_C,
        // control:
    ParU_C_Control *Control_C)
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
ParU_Ret ParU_C_Solve_Axx (
    // input:
    ParU_C_Symbolic *Sym_C, ParU_C_Numeric *Num_C,
    // input/output:
    double *b,
    // control:
    ParU_C_Control *Control_C);
//-------- Ax = b --------------------------------------------------------------
ParU_Ret ParU_C_Solve_Axb (
    // input:
    ParU_C_Symbolic *Sym_C, ParU_C_Numeric *Num_C, double *b,
    // output
    double *x,
    // control:
    ParU_C_Control *user_Control_C)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}


//-------- AX = B  (X is overwritten on B, multiple rhs)------------------------
ParU_Ret ParU_C_Solve_AXX (
    // input
    ParU_C_Symbolic *Sym_C, ParU_C_Numeric *Num_C, Int nrhs,
    // input/output:
    double *B,  // m(num_rows of A) x nrhs
    // control:
    ParU_C_Control *Control_C)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}


//-------- AX = B  (multiple rhs)-----------------------------------------------
ParU_Ret ParU_C_Solve_AXB (
    // input
    ParU_C_Symbolic *Sym_C, ParU_C_Numeric *Num_C, Int nrhs, double *B,
    // output:
    double *X,
    // control:
    ParU_C_Control *Control_C)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}

//------------------------------------------------------------------------------
//-------------- computing residual --------------------------------------------
//------------------------------------------------------------------------------
// The user provide both x and b
// resid = norm1(b-A*x) / norm1(A)
ParU_Ret ParU_C_Residual_bAx (
    // inputs:
    cholmod_sparse *A, double *x, double *b, Int m,
    // output:
    double &resid, double &anorm,
    // control:
    ParU_C_Control *Control_C)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}


// resid = norm1(B-A*X) / norm1(A) (multiple rhs)
ParU_Ret ParU_C_Residual_BAX (
    // inputs:
    cholmod_sparse *A, double *X, double *B, Int m, Int nrhs,
    // output:
    double &resid, double &anorm,
    // control:
    ParU_C_Control *Control_C)
{ 
    //TODO: copy the inside of the Control and then 
    // call the Cpp ParU_Solve
}

//------------------------------------------------------------------------------
//------------ Free routines----------------------------------------------------
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Freenum (
        ParU_C_Numeric **Num_handle_C, ParU_C_Control *Control_C)
{ 
    //TODO: // No need to copy the inside of the Control and then 
    // free the numeric object
}

 
ParU_Ret ParU_C_Freesym (
        ParU_C_Symbolic **Sym_handle_C, ParU_C_Control *Control_C)
{ 
    //TODO: // No need to copy the inside of the Control and then 
    // free the symblic object
}


} //extern c
