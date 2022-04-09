# SuiteSparse:ParU

ParU, Mohsen Aznaveh and Timothy A. Davis, (c) 2022, All Rights Reserved.

SPDX-License-Identifier: GNU GPL 3.0

--------------------------------------------------------------------------------

## Introduction

Paru: is a set of routings for solving sparse linear system via parallel
multifrontal LU factorization algorithms.  Requires OpenMP 4.0+, BLAS, CHOLMOD,
UMFPACK, AMD, COLAMD and METIS.

##  How to use

You should include ParU.hpp in your C++ project. Then for solving Ax=b
In which A is a sparse matrix in matrix market format with double entries
and b is a dense vector of double 
(or a dense matrix B for multiple rhs)

     // you can have different Controls for each
     info = ParU_Analyze(A, &Sym, &Control);
     // you can have multiple different factorization with a single ParU_Analyze
     info = ParU_Factorize(A, Sym, &Num, &Control);
     info = ParU_Solve(Sym, Num, b, x, &Control);
     ParU_Freenum(Sym, &Num, &Control);
     ParU_Freesym(&Sym, &Control);

See Demo fore more examples


--------------------------------------------------------------------------------
## License
Copyright (C) 2022 Mohsen Aznaveh and Timothy A. Davis

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------
