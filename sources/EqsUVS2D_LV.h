// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ U V S 2 D _ LV
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsUVS2D_LV.h   : definition file of the class.
// EqsUVS2D_LV.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for shallow water flow with linear
// shape functions for velocities and constant shape functions for the water elevation.
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
//  01.01.200x    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_UVS2D_LV_INCL
#define EQS_UVS2D_LV_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_UVS2D_LV : public EQS
{
  protected:
    double** UElimEq;
    double** VElimEq;
    double** PElimEq;
    double   relaxThdt_UV;
    double   relaxThdt_H;

  public:
    EQS_UVS2D_LV();
    ~EQS_UVS2D_LV();

    void Execute( PROJECT*, int );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Bound( ELEM*, PROJECT*, double**, double*, int, int, int, int, int );
    void Region( ELEM*, PROJECT*, double**, double* );
    void RegionBT( ELEM*, PROJECT*, double**, double* );
    void RegionAI( ELEM*, PROJECT*, double**, double* );
    int  PToNode( PROJECT* );
};

#endif
