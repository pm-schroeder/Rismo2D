// ======================================================================================
//                                 P R E C O _ I L U T
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

#ifndef PRECO_ILUT_INCL
#define PRECO_ILUT_INCL

#include "Precon.h"

class FROMAT;


class PRECO_ILUT : public PRECON
{
  public:
    FROMAT* fromat;

  public:
    PRECO_ILUT();
    ~PRECO_ILUT();

    void Factor( PROJECT* project, EQS* eqs, CRSMAT* crsm );
    void Solve( PROJECT* project, EQS* eqs, double* B, double* X );
};
#endif
