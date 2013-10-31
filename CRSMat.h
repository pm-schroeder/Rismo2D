// ======================================================================================
//                                    C R S M A T
// ======================================================================================
// This class implements a matrix object with compressed row storage (CRS).
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
// 01.01.1994     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef CRSMAT_INCL
#define CRSMAT_INCL

#include "Defs.h"

class EQS;
class SUBDOM;
class MODEL;
class PROJECT;
class FROMAT;


class CRSMAT
{
  public:
    int      m_neq;
    int      m_ceq;

    int      m_entries;

    int      m_neq_up;       // used for MPI; see class EQS
    int      m_neq_dn;

    int*     m_width;
    int**    m_index;
    REALPR** m_A;

  protected:
    int      m_user;
    int      m_buffer;
    int      m_shadow;

    enum     { k_nbuf = 100 };
    int      m_nbuf;
    long     m_bufsz;
    long*    m_size;
    int**    m_Ibuf;
    REALPR** m_Abuf;

  public:
    // CRSMat.cpp ------------------------------------------------------------------------
    CRSMAT();
    CRSMAT( int n, int m=0, int entries=0 );
    CRSMAT( CRSMAT* A );

    virtual ~CRSMAT();

    void    Alloc_Buf( int nb=1 );
    void    Alloc_I();
    void    Alloc_A();

    void    Init();

    void    Append( int eq, double pivot, FROMAT* fromat );

    void    Addrow( int from, int to, REALPR factor= (REALPR)1.0 );

    double  ScaleL2Norm( double* B, SUBDOM* subdom );
    void    ScaleDiag( double* B );

    int     Getneq();
    int     Getwidth( int r );
    int     Getindex( int r, int i );

    void    Setwidth( int* w );
    void    Setwidth( int r, int w );
    void    Setindex( int r, int i, int c );

    int     Find( int r, int c );

    REALPR  Get( int r, int i );
    void    Set( int r, int i, REALPR v );

    double* MulVec( double* x, double* r, PROJECT* project, EQS* eqs );

    // Assemble.cpp ----------------------------------------------------------------------
    void    AssembleEstifm_im( EQS* eqs, MODEL* m, PROJECT* p );
    void    AssembleEqs_im( EQS* eqs, double* rhs, MODEL* m, PROJECT* p );
    void    AssembleForce( EQS* eqs, double* rhs, MODEL* m, PROJECT* p );
};

#endif
