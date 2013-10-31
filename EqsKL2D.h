// ======================================================================================
//                                   E Q S _ K L 2 D
// ======================================================================================
// This class implements a differential equation system for the k-equation.
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
// 01.01.1998     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef EQS_KL2D_INCL
#define EQS_KL2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_KL2D : public EQS
{
  protected:
    double relaxThdt_KD;

  public:
    EQS_KL2D();
    ~EQS_KL2D();

    void Execute( PROJECT*, int );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );
    void Dissipation( PROJECT* );
};

#endif
