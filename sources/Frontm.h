// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// F R O N T M
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Frontm.h   : definition file of the class.
// Frontm.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the frontal solver for assembled matrices.
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
//  01.01.1994    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FRONTM_INCL
#define FRONTM_INCL

#include "Defs.h"
#include "Front.h"

class EQL;
class FROMAT;

class FRONTM : public FRONT
{
  public:
    FRONTM();
    virtual ~FRONTM();

    virtual void Direct( PROJECT* prj, MODEL* m, EQS* eqs, double* vec )
    { };
    virtual void Direct( PROJECT* prj, CRSMAT* M, double* rhs, double* x );

    EQL* Eliminate( FROMAT* fromat, EQL* eqtop, int neq_up, CRSMAT* solve,
                    int* width, int** index, REALPR** A, double* vec );

    int BiCGStab( PROJECT* project, FROMAT* fromat, double* B, double* X );
};
#endif
