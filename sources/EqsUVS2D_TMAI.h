// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ U V S 2 D _ T M A I
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsUVS2D_TMAI.h     : definition file of the class.
// EqsUVS2D_TMAI.cpp   : implementation file of the class.
//
// CoefsUVS2D_TMAI.cpp : methods EQS_UVS2D_TMAI::Coefs()
//                               EQS_UVS2D_TMAI::Region()
//                               EQS_UVS2D_TMAI::Region_dt()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for shallow water flow with a modified
// time discretization with better convergence behaviour if time steps are very small (H-LES)
// using an anisotropic turbulence closure (Elder).
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
//  03.01.2013    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_UVS2D_TMAI_INCL
#define EQS_UVS2D_TMAI_INCL

#include "EqsUVS2D_TM.h"


class EQS_UVS2D_TMAI : public EQS_UVS2D_TM
{
  public:
    EQS_UVS2D_TMAI();
    ~EQS_UVS2D_TMAI();

  protected:
    int  Coefs( ELEM*, PROJECT*, double**, double* );
    void Region( ELEM*, PROJECT*, double**, double* );
    void Region_dt( ELEM*, PROJECT*, double**, double* );
};
#endif
