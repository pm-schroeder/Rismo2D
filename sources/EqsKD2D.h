// ======================================================================================
//                                   E Q S _ K D 2 D
// ======================================================================================
// This class implements a differential equation system for the
// k-epsilon (turbulence) equations.
// ======================================================================================
//
// Copyright (C) 1992-2010  by  P.M. SCHROEDER
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

#ifndef EQS_KD2D_INCL
#define EQS_KD2D_INCL

#include "Eqs.h"


class  DRYREW;
class  ELEM;
class  MODEL;
class  NODE;
class  PROJECT;
class  TIMEINT;


class EQS_KD2D : public EQS
{
  public:
    int     linearShape;
    int     quarterShape;

  protected:
    double  relaxThdt_KD;
    NODE*   cbuf;
    NODE**  cent;

  public:
    EQS_KD2D();
    ~EQS_KD2D();

    void Execute( PROJECT*, int, int );
    int  Coefs( ELEM*, PROJECT*, double**, double* );

  protected:
    void Region( ELEM*, PROJECT*, double**, double* );
    void RegionAI( ELEM*, PROJECT*, double**, double* );

    void RegionL( ELEM*, PROJECT*, double**, double* );
    void RegionLAI( ELEM*, PROJECT*, double**, double* );

    void RegionQ( ELEM*, PROJECT*, double**, double*, int*, int* );
    void RegionQAI( ELEM*, PROJECT*, double**, double* );

    void Validate( PROJECT*, int, NODE**, int =0, NODE** =NULL, ELEM** =NULL );
};

#endif
