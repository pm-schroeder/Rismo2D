// ======================================================================================
//                                  E Q S _ U V S 2 D
// ======================================================================================
// This class implements a differential equation system for shallow water flow.
// ======================================================================================
//
// Copyright (C) 1992-2012  by  P.M. SCHROEDER
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
// 01.01.1992     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef EQS_UVS2D_AI_INCL
#define EQS_UVS2D_AI_INCL

#include "EqsUVS2D.h"


class  ELEM;
class  PROJECT;


class EQS_UVS2D_AI : public EQS_UVS2D
{
  public:
    EQS_UVS2D_AI();
    ~EQS_UVS2D_AI();

    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );
    void Region_40306( ELEM*, PROJECT*, double**, double* );
};
#endif
