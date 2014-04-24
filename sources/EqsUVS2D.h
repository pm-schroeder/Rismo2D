// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ U V S 2 D
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsUVS2D.h     : definition file of the class.
// EqsUVS2D.cpp   : implementation file of the class.
//
// CoefsUVS2D.cpp : methods EQS_UVS2D::Coefs()
//                          EQS_UVS2D::Bound()
//                          EQS_UVS2D::Region()
//                          EQS_UVS2D::Bound_pinc()
//                          EQS_UVS2D::Region_pinc()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for shallow water flow.
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
//  01.09.1992    sc    first implementation / first concept
//  01.11.1999    sc
//  02.01.2010    ko    implementation of diffuse Q at Nodes (Source/Sink)
//                      in Region/RegionAI
//  14.01.2011    sc    experimental: implementation of methods Bound_pinc() and
//                      Region_pinc() with partially integrated convective terms
//  21.01.2011    sc    Error detected in Region(): Wrong initialization in formulation
//                      of the element stiffness matrix (H-derivative of y-momentum):
//                      #     ifdef kBoussinesq2D
//                            ...
//                      #     else
//                            dfxx  =  0.0;
//                            dfyx  =  0.0;   (should be dfxy = 0.0;)
//                            ...
//  18.08.2012    sc     The anisotrop method EQS_UVS2D::RegionAI is now implemented
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_UVS2D_INCL
#define EQS_UVS2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_UVS2D : public EQS
{
  protected:
    double  relaxThdt_UV;
    double  relaxThdt_H;

    // dispersion terms
    double *Duu;
    double *Dvv;
    double *Duv;

  public:
    EQS_UVS2D();
    ~EQS_UVS2D();

    virtual void Execute( PROJECT*, int );
    virtual void Predict( PROJECT*, int, double, double );
    virtual void Timegrad( PROJECT*, double, double );

  protected:
    virtual int  Coefs( ELEM*, PROJECT*, double**, double* );
    virtual void Bound( ELEM*, PROJECT*, double**, double* );
    virtual void Region( ELEM*, PROJECT*, double**, double* );

    virtual void Bound_pinc( ELEM*, PROJECT*, double**, double* );
    virtual void Region_pinc( ELEM*, PROJECT*, double**, double* );
};
#endif
