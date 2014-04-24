// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S _ K D 2 D
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// EqsKD2D.h     : definition file of the class.
// EqsKD2D.cpp   : implementation file of the class.
//
// CoefsKD2D.cpp : methods EQS_KD2D::Coefs()
//                         EQS_KD2D::Region()
//                         EQS_KD2D::RegionAI()
//                         EQS_KD2D::RegionL()
//                         EQS_KD2D::RegionLAI()
//                         EQS_KD2D::RegionQ()
//                         EQS_KD2D::RegionQAI()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a differential equation system for the
// k-epsilon (turbulence) equations.
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

#ifndef EQS_KD2D_INCL
#define EQS_KD2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_KD2D : public EQS
{
  public:
    int     linearShape;
    int     quarterShape;

  protected:
    double  relaxThdt_KD;
    NODE*   cbuf;
    NODE**  cent;

  public:
    EQS_KD2D();
    ~EQS_KD2D();

    void Execute( PROJECT*, int, int );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );
    void RegionAI( ELEM*, PROJECT*, double**, double* );

    void RegionL( ELEM*, PROJECT*, double**, double* );
    void RegionLAI( ELEM*, PROJECT*, double**, double* );

    void RegionQ( ELEM*, PROJECT*, double**, double*, int*, int* );
    void RegionQAI( ELEM*, PROJECT*, double**, double* );

    void Validate( PROJECT*, int, NODE**, int =0, NODE** =NULL, ELEM** =NULL );
};

#endif
