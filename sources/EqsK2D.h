// =======================================================================================
//                                   E Q S _ K 2 D
// =======================================================================================
// This class implements a differential equation system for the (turbulent) k-equation.
// =======================================================================================
//
// Copyright (C) 1992-2010  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// =======================================================================================
//
//    date               description
// ----------   ------   -----------------------------------------------------------------
// 09.04.2010     sc     first implementation / first concept
//
// =======================================================================================

// read this header file only once

#ifndef EQS_K2D_INCL
#define EQS_K2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_K2D : public EQS
{
  protected:
    double  relaxThdt_KD;

  public:
    EQS_K2D();
    ~EQS_K2D();

    void Execute( PROJECT*, int, int );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );

    void Validate( int, NODE**, PROJECT* );
};

#endif
