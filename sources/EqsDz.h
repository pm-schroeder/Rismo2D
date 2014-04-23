// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ D Z
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsDz.h   : definition file of the class.
// EqsDz.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for bottom evolution.
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
//  11.05.2005    sc    first implementation / first concept
//  05.05.2006    sc    splitting bottom evolution from bed load module
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_DZ_INCL
#define EQS_DZ_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  SHAPE;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_DZ : public EQS
{
  public:
    enum
    {
      kLinear = 1, kQuadratic = 2
    };
    enum
    {
      kBottomEvol, kBottomEvol_L
    };
    int     coefs;

    int     initDzQb;
    double  factDzQb;

    double* etaQb;
    double* dzdt;

  public:
    EQS_DZ();
    ~EQS_DZ();

    void   Execute( PROJECT* project, int shape );
    double MorphTime( PROJECT* project, int n, double* dz );
    void   QbDiff( PROJECT* project, int shape );

  protected:
    int  Coefs( ELEM*, PROJECT*, double**, double* );
    void Bound( ELEM*, PROJECT*, double**, double*, SHAPE* );
    void Region( ELEM*, PROJECT*, double**, double*, SHAPE* );
};

#endif
