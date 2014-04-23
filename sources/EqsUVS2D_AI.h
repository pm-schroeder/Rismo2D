// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ U V S 2 D _ A I
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsUVS2D_AI.h   : definition file of the class.
// EqsUVS2D_AI.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for shallow water flow with an
// anisotrop turbulence closure (Elder).
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

#ifndef EQS_UVS2D_AI_INCL
#define EQS_UVS2D_AI_INCL

#include "EqsUVS2D.h"


class  ELEM;
class  PROJECT;


class EQS_UVS2D_AI : public EQS_UVS2D
{
  public:
    EQS_UVS2D_AI();
    ~EQS_UVS2D_AI();

    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );
    void Region_40306( ELEM*, PROJECT*, double**, double* );
};
#endif
