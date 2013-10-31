// ======================================================================================
//                              E Q S _ U V S 2 D _ L V
// ======================================================================================
// This class implements a differential equation systems for shallow water flow.
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
// 01.01.2000     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef EQS_UVS2D_LV_INCL
#define EQS_UVS2D_LV_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_UVS2D_LV : public EQS
{
  protected:
    double** UElimEq;
    double** VElimEq;
    double** PElimEq;
    double   relaxThdt_UV;
    double   relaxThdt_H;

  public:
    EQS_UVS2D_LV();
    ~EQS_UVS2D_LV();

    void Execute( PROJECT*, int );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Bound( ELEM*, PROJECT*, double**, double*, int, int, int, int, int );
    void Region( ELEM*, PROJECT*, double**, double* );
    void RegionBT( ELEM*, PROJECT*, double**, double* );
    void RegionAI( ELEM*, PROJECT*, double**, double* );
    int  PToNode( PROJECT* );
};

#endif
