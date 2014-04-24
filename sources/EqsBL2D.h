// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ B L 2 D
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsBL2D.h     : definition file of the class.
// EqsBL2D.cpp   : implementation file of the class.
//
// CoefsBL2D.cpp : methods EQS_BL2D::Coefs()
//                         EQS_BL2D::Bound()
//                         EQS_BL2D::Region()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for bed load equations.
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
//  11.04.2005    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_BL2D_INCL
#define EQS_BL2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  SHAPE;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_BL2D : public EQS
{
  public:
    enum
    {
      kLinear = 1, kQuadratic = 2
    };
    int coefs;

  public:
    EQS_BL2D();
    ~EQS_BL2D();

    void Execute( PROJECT* project, int shape );

  protected:
    int  Coefs( ELEM*, PROJECT*, double**, double* );

    void Bound( ELEM*, PROJECT*, double**, double*, SHAPE* );
    void Region( ELEM*, PROJECT*, double**, double*, SHAPE* );

//    void BoundDiffusion( ELEM*, PROJECT*, double**, double*, SHAPE* );
//    void RegionDiffusion( ELEM*, PROJECT*, double**, double*, SHAPE* );
};

#endif
