// ======================================================================================
//                                    P R E C O N
// ======================================================================================
// This class implements the base class for preconditioner.
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

#ifndef PRECON_INCL
#define PRECON_INCL

class PRECON
{
  public:
    CRSMAT* crsi;

  public:
    PRECON()
    {
      crsi = NULL;
    };

    virtual ~PRECON()
    {
      if( crsi )  delete crsi;
    };

    virtual void Factor( PROJECT* project, EQS* eqs, CRSMAT* crsm ) = 0;
    virtual void Solve( PROJECT* project, EQS* eqs, double* B, double* X ) = 0;
};
#endif
