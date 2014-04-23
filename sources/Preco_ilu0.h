// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// P R E C O _ I L U 0
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Preco_ilu0.h   : definition file of the class.
// Preco_ilu0.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a preconditioner for iterative solvers.
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
//  01.01.2004    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PRECO_ILU0_INCL
#define PRECO_ILU0_INCL

#include "Precon.h"

class PRECO_ILU0 : public PRECON
{
  public:
    PRECO_ILU0();
    ~PRECO_ILU0();

    void Factor( PROJECT* project, EQS* eqs, CRSMAT* crsm );
    void Solve( PROJECT* project, EQS* eqs, double* B, double* X );
};
#endif
