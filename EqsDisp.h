// ======================================================================================
//                                  E Q S _ D I S P
// ======================================================================================
// This class implements a differential equation system for curvature (curv).
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
// 25.07.2006     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef EQS_DISP_INCL
#define EQS_DISP_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;
struct KDCONST;


class EQS_DISP : public EQS
{
  public:
    EQS_DISP();
    ~EQS_DISP();

    void Execute( PROJECT* );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );
};
#endif
