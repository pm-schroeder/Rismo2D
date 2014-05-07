// ======================================================================================
//                                       B I C G S T A B
// ======================================================================================
// This class implements an iterative solver for a system of linear equations.
// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software.
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.1994     sc     first implementation
//
// ======================================================================================

// read this header file only once

#ifndef BICGSTAB_INCL
#define BICGSTAB_INCL

#include "Defs.h"
#include "Solver.h"


class BICGSTAB : public SOLVER
{
  public:
    BICGSTAB();
    virtual ~BICGSTAB();

   virtual int Iterate( PROJECT* prj, CRSMAT* M, double* B, double* X, PRECON* P );
};
#endif
