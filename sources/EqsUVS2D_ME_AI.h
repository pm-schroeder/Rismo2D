// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ U V S 2 D _ M E _ A I
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsUVS2D_ME_AI.h     : definition file of the class.
// EqsUVS2D_ME_AI.cpp   : implementation file of the class.
//
// CoefsUVS2D_ME_AI.cpp : methods EQS_UVS2D_ME_AI::Coefs()
//                                EQS_UVS2D_ME_AI::Bound()
//                                EQS_UVS2D_ME_AI::Region()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for shallow water flow. The element
// shape functions are based on the MINI-element with a bubble shape functions for velocities
// and linear shape functions for the water elevation.
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
//  01.08.2014    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_UVS2D_ME_AI_INCL
#define EQS_UVS2D_ME_AI_INCL

#include "EqsUVS2D_ME.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_UVS2D_ME_AI : public EQS_UVS2D_ME
{
  protected:
    double **UElimEq;
    double **VElimEq;

    double   relaxThdt_UV;
    double   relaxThdt_H;

  public:
    EQS_UVS2D_ME_AI();
    virtual ~EQS_UVS2D_ME_AI();

  protected:
    virtual void Region( ELEM*, PROJECT*, double**, double* );
};

#endif
