// ======================================================================================
//                                 P _ B C G S T A B
// ======================================================================================
// This class implements an iterative solver for linear equation systems (PARMS).
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

#ifndef P_BCGSTABD_INCL
#define P_BCGSTABD_INCL

#include "Defs.h"
#include "Parms.h"


class P_BCGSTABD : public PARMS
{
  public:
    P_BCGSTABD();
    virtual ~P_BCGSTABD();

   virtual int Iterate( PROJECT* prj, CRSMAT* M, double* B, double* X, PRECON* P );
};
#endif
