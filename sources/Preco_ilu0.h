// ======================================================================================
//                                 P R E C O _ I L U 0
// ======================================================================================
// This class implements a preconditioner for iterative solvers.
// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.2004     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once ...

#ifndef PRECO_ILU0_INCL
#define PRECO_ILU0_INCL

#include "Precon.h"

class PRECO_ILU0 : public PRECON
{
  public:
    PRECO_ILU0();
    ~PRECO_ILU0();

    void Factor( PROJECT* project, EQS* eqs, CRSMAT* crsm );
    void Solve( PROJECT* project, EQS* eqs, double* B, double* X );
};
#endif
