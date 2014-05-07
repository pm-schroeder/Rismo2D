// ======================================================================================
//                                     F R O N T M
// ======================================================================================
// This class implements the frontal solver for assembled matrices.
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
// 01.01.1994     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef FRONTM_INCL
#define FRONTM_INCL

#include "Defs.h"
#include "Front.h"

class EQL;
class FROMAT;

class FRONTM : public FRONT
{
  public:
    FRONTM();
    virtual ~FRONTM();

    virtual void Direct( PROJECT* prj, MODEL* m, EQS* eqs, double* vec )
    { };
    virtual void Direct( PROJECT* prj, CRSMAT* M, double* rhs, double* x );

    EQL* Eliminate( FROMAT* fromat, EQL* eqtop, int neq_up, CRSMAT* solve,
                    int* width, int** index, REALPR** A, double* vec );

    int BiCGStab( PROJECT* project, FROMAT* fromat, double* B, double* X );
};
#endif
