// ======================================================================================
//                                      V A R S
// ======================================================================================
// This class implements the variables at nodes.
// ======================================================================================
//
// Copyright (C) 1992-2010  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version)
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

// read this header file only once ...

#ifndef VARS_INCL
#define VARS_INCL

#include "Defs.h"


#define kVarU        0       // defines for dependent variables
#define kVarV        1
#define kVarPhi      2
#define kVarH        3
#define kVarS        4
#define kVarK        5
#define kVarD        6
#define kVarC        7
#define kVarQb       8
#define kVarDz       9

#define kVarUH      10       // defines for virtual variables
#define kVarVH      11


class VARS
{
  public:
    double U, V;             // velocities
    double dUdt, dVdt;       // and their time derivatives

    double S;                // water elevation
    double dSdt;

    double K;                // turbulent kinetic energy (TKE)
    double D;                // dissipation of TKE
    double dKdt, dDdt;

    double C;                // suspended load transport
    double Qb;               // bed load transport

  public:
    VARS()
    {
      U    = 0.0;
      V    = 0.0;
      dUdt = 0.0;
      dVdt = 0.0;

      S    = 0.0;
      dSdt = 0.0;

      K    = 0.0;
      D    = 0.0;
      dKdt = 0.0;
      dDdt = 0.0;

      C    = 0.0;
      Qb   = 0.0;
    };

    VARS operator =( const VARS& v )
    {
      U  = v.U;
      V  = v.V;
      dUdt = v.dUdt;
      dVdt = v.dVdt;

      S  = v.S;
      dSdt = v.dSdt;

      K  = v.K;
      D  = v.D;
      dKdt = v.dKdt;
      dDdt = v.dDdt;

      C  = v.C;
      Qb = v.Qb;

      return *this;
    };
};

#endif
