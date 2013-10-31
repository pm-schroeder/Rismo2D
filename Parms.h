// ======================================================================================
//                                      P A R M S
// ======================================================================================
// This class is the base class for the PARAMS iterative solvers.
// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.2004     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once ...

#ifndef PARMS_INCL
#define PARMS_INCL

#include "Defs.h"
#include "Solver.h"


class PARMS : public SOLVER
{
  public:
    PARMS()
    {};
    ~PARMS()
    {};


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
