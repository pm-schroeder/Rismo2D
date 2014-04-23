// ======================================================================================
//                                  E Q S _ P P E 2 D
// ======================================================================================
// This class implements a differential equation system for pressure correction.
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
// 01.01.1995     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef EQS_PPE2D_INCL
#define EQS_PPE2D_INCL

#include "Eqs.h"


class  NODE;
class  ELEM;
class  MODEL;
class  PROJECT;
class  TIMEINT;
class  DRYREW;


class EQS_PPE2D : public EQS
{
  public:
    EQS_PPE2D();
    ~EQS_PPE2D();

    void Execute( PROJECT* project );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );
    void Bound( ELEM*, PROJECT*, double**, double* );
};
#endif
