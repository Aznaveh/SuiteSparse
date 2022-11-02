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
ParU_Ret ParU_C_Version(int ver[3], char date[128])
{
    return ParU_Version(ver, date);
}

// initialize C_Control with the default values
void init_control(ParU_C_Control *Control_C)
{
    Control_C->mem_chunk = 1024 * 1024;  // chunk size for memset and memcpy

    Control_C->umfpack_ordering = UMFPACK_ORDERING_METIS;
    Control_C->umfpack_strategy =
        UMFPACK_STRATEGY_AUTO;  // symmetric or unsymmetric
    Control_C->umfpack_default_singleton = 1;

    Control_C->relaxed_amalgamation_threshold = 32;

    Control_C->scale = 1;
    Control_C->panel_width = 32;
    Control_C->paru_strategy = PARU_STRATEGY_AUTO;

    Control_C->piv_toler = .1;
    Control_C->diag_toler = .001;
    Control_C->trivial = 4;
    Control_C->worthwhile_dgemm = 512;
    Control_C->worthwhile_trsm = 4096;
    Control_C->paru_max_threads = 0;
}
// copy the inside of the C structrue to the Cpp structure
void cp_control(ParU_Control *Control, ParU_C_Control *Control_C)
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
    ParU_Control Control;
    cp_control(Control, Control_C);
    ParU_Symbolic *Sym;
    ParU_Ret info;
    info = ParU_Analyze(A, &Sym, &Control);
    *Sym_handle_C = (ParU_C_Symbolic *)Sym;

    // If ParU_C_Symbolic is not defined as void what I can do is to alloc a new
    // ParU_C_Symbolic and store the pointers in the write place
    //   ParU_C_Symbolic p_Sym = paru_alloc(sizeof(ParU_C_Symbolic))
    //   //update inside the data structrue
    //   p_Sym->Sym = Sym; ...
    //   *Sym_handle_C = (ParU_C_Symbolic *) Sym;
    //
    return info;
}

//------------------------------------------------------------------------------
// ParU_Factorize: Numeric factorization is done in this routine. Scaling and
// making Sx matrix, computing factors and permutations is here. ParU_C_Symbolic
// structure is computed ParU_Analyze and is an input in this routine.
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Factorize(
    // input:
    cholmod_sparse *A, ParU_C_Symbolic *Sym_C,
    // output:
    ParU_C_Numeric **Num_handle_C,
    // control:
    ParU_C_Control *Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    ParU_Symbolic *Sym = (ParU_C_Symbolic *)Sym_C;
    ParU_Numeric *Num;
    ParU_Ret info;
    info = ParU_Factorize(A, Sym, &Num, &Control);
    *Num_handle_C = (ParU_C_Numeric *)Num;
    return info;
}

//------------------------------------------------------------------------------
//--------------------- Solve routines -----------------------------------------
//------------------------------------------------------------------------------
// In all the solve routines Num structure must come with the same Sym struct
// that comes from ParU_Factorize
//-------- Ax = b (x is overwritten on b)---------------------------------------
ParU_Ret ParU_C_Solve_Axx(
    // input:
    ParU_C_Symbolic *Sym_C, ParU_C_Numeric *Num_C,
    // input/output:
    double *b,
    // control:
    ParU_C_Control *Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    return ParU_Solve((ParU_C_Symbolic *)Sym_C, (ParU_Numeric *)ParU_C_Numeric,
                      b, Control);
}
//-------- Ax = b --------------------------------------------------------------
ParU_Ret ParU_C_Solve_Axb(
    // input:
    ParU_C_Symbolic *Sym_C, ParU_C_Numeric *Num_C, double *b,
    // output
    double *x,
    // control:
    ParU_C_Control *user_Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    return ParU_Solve((ParU_C_Symbolic *)Sym_C, (ParU_Numeric *)ParU_C_Numeric,
                      b, x, Control);
}

//-------- AX = B  (X is overwritten on B, multiple rhs)------------------------
ParU_Ret ParU_C_Solve_AXX(
    // input
    ParU_C_Symbolic *Sym_C, ParU_C_Numeric *Num_C, Int nrhs,
    // input/output:
    double *B,  // m(num_rows of A) x nrhs
    // control:
    ParU_C_Control *Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    return ParU_Solve((ParU_C_Symbolic *)Sym_C, (ParU_Numeric *)ParU_C_Numeric,
                      B, Control);
}

//-------- AX = B  (multiple rhs)-----------------------------------------------
ParU_Ret ParU_C_Solve_AXB(
    // input
    ParU_C_Symbolic *Sym_C, ParU_C_Numeric *Num_C, Int nrhs, double *B,
    // output:
    double *X,
    // control:
    ParU_C_Control *Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    return ParU_Solve((ParU_C_Symbolic *)Sym_C, (ParU_Numeric *)ParU_C_Numeric,
                      B, X, Control);
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
    ParU_C_Control *Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    return ParU_Residual(A, x, b, m, resid, anorm, Control);
}

// resid = norm1(B-A*X) / norm1(A) (multiple rhs)
ParU_Ret ParU_C_Residual_BAX(
    // inputs:
    cholmod_sparse *A, double *X, double *B, Int m, Int nrhs,
    // output:
    double &resid, double &anorm,
    // control:
    ParU_C_Control *Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    return ParU_Residual(A, X, B, m, nrhs, resid, anorm, Control);
}

//------------------------------------------------------------------------------
//------------ Free routines----------------------------------------------------
//------------------------------------------------------------------------------
ParU_Ret ParU_C_Freenum(ParU_C_Numeric **Num_handle_C,
                        ParU_C_Control *Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    return ParU_Freenum((ParU_C_Numeric **)Num_handle_C, Control);
}

ParU_Ret ParU_C_Freesym(ParU_C_Symbolic **Sym_handle_C,
                        ParU_C_Control *Control_C)
{
    ParU_Control Control;
    cp_control(Control, Control_C);
    return ParU_Freesym((ParU_C_Symbolic **)Sym_handle_C, Control);
}

} //extern c
