// ======================================================================================
//                              E Q S _ U V S 2 D _ T M A I
// ======================================================================================
// This class implements a differential equation system for shallow water flow.
//
// This class is used for strongly time dependent flows. For this purpose, the
// time integration has been differently implemented to the class EQS_UVS2D. 
// ======================================================================================
//
// Copyright (C) 1992-2013  by  P.M. SCHROEDER
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
// 03.01.2013     sc     first implementation as derived class from EQS_UVSD_TM
//
// ======================================================================================

// read this header file only once

#ifndef EQS_UVS2D_TMAI_INCL
#define EQS_UVS2D_TMAI_INCL

#include "EqsUVS2D_TM.h"


class EQS_UVS2D_TMAI : public EQS_UVS2D_TM
{
  public:
    EQS_UVS2D_TMAI();
    ~EQS_UVS2D_TMAI();

  protected:
    int  Coefs( ELEM*, PROJECT*, double**, double* );
    void Region( ELEM*, PROJECT*, double**, double* );
    void Region_dt( ELEM*, PROJECT*, double**, double* );
};
#endif
