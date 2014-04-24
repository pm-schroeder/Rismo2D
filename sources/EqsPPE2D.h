// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ P P E 2 D
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsPPE2D.h     : definition file of the class.
// EqsPPE2D.cpp   : implementation file of the class.
//
// CoefsPPE2D.cpp : methods EQS_PPE2D::Coefs()
//                          EQS_PPE2D::Bound()
//                          EQS_PPE2D::Region()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for pressure correction.
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
//  01.01.1995    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_PPE2D_INCL
#define EQS_PPE2D_INCL

#include "Eqs.h"


class  NODE;
class  ELEM;
class  MODEL;
class  PROJECT;
class  TIMEINT;
class  DRYREW;


class EQS_PPE2D : public EQS
{
  public:
    EQS_PPE2D();
    ~EQS_PPE2D();

    void Execute( PROJECT* project );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );
    void Bound( ELEM*, PROJECT*, double**, double* );
};
#endif
