// ======================================================================================
//                                   E Q S _ S L 2 D
// ======================================================================================
// This class implements a differential equation system for the
// advection-diffusion equation.
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
// 01.01.1992     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef EQS_SL2D_INCL
#define EQS_SL2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_SL2D : public EQS
{
  public:
    EQS_SL2D();
    ~EQS_SL2D();

    void Execute( PROJECT* );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Bound( ELEM*, PROJECT*, double**, double* );
    void Region( ELEM*, PROJECT*, double**, double* );
};

#endif
