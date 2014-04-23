// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// P A R M S
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Parms.h   : definition file of the class.
// Parms.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class is the base class for the PARAMS iterative solvers.
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

#ifndef PARMS_INCL
#define PARMS_INCL

#include "Defs.h"
#include "Solver.h"


class PARMS : public SOLVER
{
  public:
    PARMS()
    {}
    ~PARMS()
    {}


   // replaced BLAS functions ------------------------------------------------------------

   void daxpy( int n, double a, double* x, double* y )      // y[1,n] += a * x[1,n]
   {
     for( int i=0; i<n; i++ )  y[i] += a * x[i];
   }

   void dscal( int n, double a, double* x )                 // x[1,n] *= a
   {
     for( int i=0; i<n; i++ )  x[i] *= a;
   }

   void dcopy( int n, double* x, double* y )                // y[1,n] = x[1,n]
   {
     for( int i=0; i<n; i++ )  y[i] = x[i];
   }

   double ddot( int n, double* x, double* y )               // d = x[1,n] * y[1,n]
   {
     double d = 0.0;
     for( int i=0; i<n; i++ )  d += x[i] * y[i];
     return d;
   }
};
#endif
