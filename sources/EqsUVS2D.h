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

#ifndef EQS_UVS2D_INCL
#define EQS_UVS2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_UVS2D : public EQS
{
  protected:
    double  relaxThdt_UV;
    double  relaxThdt_H;

    // dispersion terms
    double *Duu;
    double *Dvv;
    double *Duv;

  public:
    EQS_UVS2D();
    ~EQS_UVS2D();

    virtual void Execute( PROJECT*, int );
    virtual void Predict( PROJECT*, int, double, double );
    virtual void Timegrad( PROJECT*, double, double );

  protected:
    virtual int  Coefs( ELEM*, PROJECT*, double**, double* );
    virtual void Bound( ELEM*, PROJECT*, double**, double* );
    virtual void Region( ELEM*, PROJECT*, double**, double* );

    virtual void Bound_pinc( ELEM*, PROJECT*, double**, double* );
    virtual void Region_pinc( ELEM*, PROJECT*, double**, double* );
};
#endif
