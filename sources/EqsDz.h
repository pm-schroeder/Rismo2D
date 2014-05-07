// ======================================================================================
//                                     E Q S _ D Z
// ======================================================================================
// This class implements a differential equation system for bottom evolution.
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
// 11.05.2005     sc     first implementation / first concept
// 05.05.2006     sc     splitting bottom evolution from bed load module
//
// ======================================================================================

// read this header file only once

#ifndef EQS_DZ_INCL
#define EQS_DZ_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  SHAPE;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_DZ : public EQS
{
  public:
    enum
    {
      kLinear = 1, kQuadratic = 2
    };
    enum
    {
      kBottomEvol, kBottomEvol_L
    };
    int     coefs;

    int     initDzQb;
    double  factDzQb;

    double* etaQb;
    double* dzdt;

  public:
    EQS_DZ();
    ~EQS_DZ();

    void   Execute( PROJECT* project, int shape );
    double MorphTime( PROJECT* project, int n, double* dz );
    void   QbDiff( PROJECT* project, int shape );

  protected:
    int  Coefs( ELEM*, PROJECT*, double**, double* );
    void Bound( ELEM*, PROJECT*, double**, double*, SHAPE* );
    void Region( ELEM*, PROJECT*, double**, double*, SHAPE* );
};

#endif
