// ======================================================================================
//                                      F R O N T
// ======================================================================================
// This class implements the frontal solver.
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
// 01.01.1992     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once ...

#ifndef FRONT_INCL
#define FRONT_INCL

#include "Defs.h"
#include "Solver.h"


class FRONT : public SOLVER
{
  public:
    // -------------------------- scratch file (frontal method) -------------------
    FILE*  id;                 // file pointers
    char   scratch[160];       // scratch file name

    int    count;              // equation counter during writing
    int    recno;              // number of record actually in buffer

    int    recLen;             // length of one record in matrix coefficients


  public:
    FRONT();
    virtual ~FRONT();

    virtual void Direct( PROJECT* prj, MODEL* m, EQS* eqs, double* vec );
    virtual void Direct( PROJECT* prj, CRSMAT* M, double* rhs, double* x )
    { };

    void ReadEq( int, size_t, long*, int*, REALPR* );
    void WriteEq( long, size_t, int*, REALPR* );

    void OpenEquation();
    void CloseEquation();
};
#endif
