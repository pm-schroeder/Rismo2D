// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// V A R S
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Vars.h   : definition file of the class.
// Vars.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the variables at nodes.
//
// -------------------------------------------------------------------------------------------------
//
// COPYRIGHT (C) 2011 - 2014  by  P.M. SCHROEDER  (sc)
//
// This program is free software; you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with this program; if
// not, write to the
//
// Free Software Foundation, Inc.
// 59 Temple Place
// Suite 330
// Boston
// MA 02111-1307 USA
//
// -------------------------------------------------------------------------------------------------
//
// P.M. Schroeder
// Walzbachtal / Germany
// michael.schroeder@hnware.de
//
// -------------------------------------------------------------------------------------------------
//
// HISTORY
//
//    date              changes
// ------------  ----  -----------------------------------------------------------------------------
//  01.01.1992    sc    First implementation of C++ class.
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

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
