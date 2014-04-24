// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// C R S M A T
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// CRSMat.h     : definition file of the class.
// CRSMat.cpp   : implementation file of the class.
//
// Assemble.cpp : methods CRSMAT::AssembleEstifm_im()
//                        CRSMAT::AssembleEqs_im()
//                        CRSMAT::AssembleForce()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a matrix object with compressed row storage (CRS).
//
// -------------------------------------------------------------------------------------------------
//
// COPYRIGHT (C) 2011 - 2014  by  P.M. SCHROEDER  (sc)
//
// This program is free software; you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with this program; if
// not, write to the
//
// Free Software Foundation, Inc.
// 59 Temple Place
// Suite 330
// Boston
// MA 02111-1307 USA
//
// -------------------------------------------------------------------------------------------------
//
// P.M. Schroeder
// Walzbachtal / Germany
// michael.schroeder@hnware.de
//
// -------------------------------------------------------------------------------------------------
//
// HISTORY
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.1994     sc     first implementation
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

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
