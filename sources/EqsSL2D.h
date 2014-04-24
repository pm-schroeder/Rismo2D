// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ S L 2 D
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsSL2D.h     : definition file of the class.
// EqsSL2D.cpp   : implementation file of the class.
//
// CoefsSL2D.cpp : methods EQS_SL2D::Coefs()
//                         EQS_SL2D::Bound()
//                         EQS_SL2D::Region()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for the
// advection-diffusion equation (suspended load).
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
//  01.01.1992    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_SL2D_INCL
#define EQS_SL2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_SL2D : public EQS
{
  public:
    EQS_SL2D();
    ~EQS_SL2D();

    void Execute( PROJECT* );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Bound( ELEM*, PROJECT*, double**, double* );
    void Region( ELEM*, PROJECT*, double**, double* );
};

#endif
