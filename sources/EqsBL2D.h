// =======================================================================================
//                                     E Q S _ B L
// =======================================================================================
// This class implements a differential equation system for bed load equations.
// =======================================================================================
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
// =======================================================================================
//
//    date               description
// ----------   ------   -----------------------------------------------------------------
// 11.04.2005     sc     first implementation / first concept
//
// =======================================================================================

// read this header file only once

#ifndef EQS_BL2D_INCL
#define EQS_BL2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  SHAPE;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_BL2D : public EQS
{
  public:
    enum
    {
      kLinear = 1, kQuadratic = 2
    };
    int coefs;

  public:
    EQS_BL2D();
    ~EQS_BL2D();

    void Execute( PROJECT* project, int shape );

  protected:
    int  Coefs( ELEM*, PROJECT*, double**, double* );

    void Bound( ELEM*, PROJECT*, double**, double*, SHAPE* );
    void Region( ELEM*, PROJECT*, double**, double*, SHAPE* );

//    void BoundDiffusion( ELEM*, PROJECT*, double**, double*, SHAPE* );
//    void RegionDiffusion( ELEM*, PROJECT*, double**, double*, SHAPE* );
};

#endif
