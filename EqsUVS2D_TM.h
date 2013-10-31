// ======================================================================================
//                                E Q S _ U V S 2 D _ T M
// ======================================================================================
// This class implements a differential equation system for shallow water flow.
//
// This class is used for strongly time dependent flows. For this purpose, the
// time integration has been differently implemented to the class EQS_UVS2D. 
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
// 06.12.2012     sc     implementation as derived class from EQS_UVSD
//
// ======================================================================================

// read this header file only once

#ifndef EQS_UVS2D_TM_INCL
#define EQS_UVS2D_TM_INCL

#include "EqsUVS2D.h"


class EQS_UVS2D_TM : public EQS_UVS2D
{
  protected:
    int  timegrad;

  public:
    EQS_UVS2D_TM();
    ~EQS_UVS2D_TM();

    virtual void Predict( PROJECT*, int, double, double );
    virtual void Timegrad( PROJECT*, double, double );

  protected:
    virtual int  Coefs( ELEM*, PROJECT*, double**, double* );

    virtual void Bound( ELEM*, PROJECT*, double**, double* );
    virtual void Region( ELEM*, PROJECT*, double**, double* );

    virtual void Bound_dt( ELEM*, PROJECT*, double**, double* );
    virtual void Region_dt( ELEM*, PROJECT*, double**, double* );
};
#endif
