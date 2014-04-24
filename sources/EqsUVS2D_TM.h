// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ U V S 2 D _ T M
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsUVS2D_TM.h     : definition file of the class.
// EqsUVS2D_TM.cpp   : implementation file of the class.
//
// CoefsUVS2D_TM.cpp : methods EQS_UVS2D_TM::Coefs()
//                             EQS_UVS2D_TM::Bound()
//                             EQS_UVS2D_TM::Region()
//                             EQS_UVS2D_TM::Bound_dt()
//                             EQS_UVS2D_TM::Region_dt()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for shallow water flow with a modified
// time discretization with better convergence behaviour if time steps are very small (H-LES).
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

#ifndef EQS_UVS2D_TM_INCL
#define EQS_UVS2D_TM_INCL

#include "EqsUVS2D.h"


class EQS_UVS2D_TM : public EQS_UVS2D
{
  protected:
    int  timegrad;

  public:
    EQS_UVS2D_TM();
    ~EQS_UVS2D_TM();

    virtual void Predict( PROJECT*, int, double, double );
    virtual void Timegrad( PROJECT*, double, double );

  protected:
    virtual int  Coefs( ELEM*, PROJECT*, double**, double* );

    virtual void Bound( ELEM*, PROJECT*, double**, double* );
    virtual void Region( ELEM*, PROJECT*, double**, double* );

    virtual void Bound_dt( ELEM*, PROJECT*, double**, double* );
    virtual void Region_dt( ELEM*, PROJECT*, double**, double* );
};
#endif
