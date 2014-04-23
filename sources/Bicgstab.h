// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// B I C G S T A B
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Bicgstab.h   : definition file of the class.
// Bicgstab.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements an iterative solver for a system of linear equations
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
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.1994     sc     first implementation
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef BICGSTAB_INCL
#define BICGSTAB_INCL

#include "Defs.h"
#include "Solver.h"


class BICGSTAB : public SOLVER
{
  public:
    BICGSTAB();
    virtual ~BICGSTAB();

   virtual int Iterate( PROJECT* prj, CRSMAT* M, double* B, double* X, PRECON* P );
};
#endif
