// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// F R O N T
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Front.h   : definition file of the class.
// Front.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the frontal solver.
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

#ifndef FRONT_INCL
#define FRONT_INCL

#include "Defs.h"
#include "Solver.h"


class FRONT : public SOLVER
{
  public:
    // -------------------------- scratch file (frontal method) -------------------
    FILE*  id;                 // file pointers
    char   scratch[160];       // scratch file name

    int    count;              // equation counter during writing
    int    recno;              // number of record actually in buffer

    int    recLen;             // length of one record in matrix coefficients


  public:
    FRONT();
    virtual ~FRONT();

    virtual void Direct( PROJECT* prj, MODEL* m, EQS* eqs, double* vec );
    virtual void Direct( PROJECT* prj, CRSMAT* M, double* rhs, double* x )
    { };

    void ReadEq( int, size_t, long*, int*, REALPR* );
    void WriteEq( long, size_t, int*, REALPR* );

    void OpenEquation();
    void CloseEquation();
};
#endif
