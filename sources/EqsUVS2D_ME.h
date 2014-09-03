// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ U V S 2 D _ M E
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsUVS2D_ME.h     : definition file of the class.
// EqsUVS2D_ME.cpp   : implementation file of the class.
//
// CoefsUVS2D_ME.cpp : methods EQS_UVS2D_ME::Coefs()
//                             EQS_UVS2D_ME::Bound()
//                             EQS_UVS2D_ME::Region()
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
//  01.01.200x    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_UVS2D_ME_INCL
#define EQS_UVS2D_ME_INCL

#include "EqsUVS2D.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_UVS2D_ME : public EQS_UVS2D
{
  protected:
    double **UElimEq;
    double **VElimEq;

    double   relaxThdt_UV;
    double   relaxThdt_H;

  public:
    EQS_UVS2D_ME();
    virtual ~EQS_UVS2D_ME();

    virtual void Execute( PROJECT*, int );
    virtual void Eliminate( ELEM*, double**, double* , int, int, int, int );

  protected:
    virtual void Bound( ELEM*, PROJECT*, double**, double* );
    virtual void Region( ELEM*, PROJECT*, double**, double* );
};

#endif
