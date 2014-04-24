// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E Q S
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Eqs.h        : definition file of the class.
// Eqs.cpp      : implementation file of the class.
//
// IndexMat.cpp : methods EQS::KillCrsm()
//                        EQS::SetIndexMat()
//                        EQS::SortIndex()
// Rotate.cpp   : method  EQS::Rotate2D()
// SetEqno.cpp  : methods EQS::SetEqno()
//                        EQS::LastEquation()
//                        EQS::ResetEqOrder()
// Solve.cpp    : method  EQS::Solve()
// Update.cpp   : method  EQS::EQS::Update()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the base class for differential equation systems.
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
//    date              changes
// ------------  ----  -----------------------------------------------------------------------------
//  01.01.200x    sc     first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EQS_INCL
#define EQS_INCL

#include "Defs.h"
#include "Solver.h"

class NODE;
class ELEM;
class GRID;
class MODEL;
class PROJECT;
class SUBDOM;
class CRSMAT;


// ---------------------------------------------------------------------------------------
// equation system

class EQS
{
  protected:
    int**           nodeEqno;           // array of node equation numbers
    int**           elemEqno;           // array of element equation numbers

    int             iterCountCG;        // total counter for CG iterations

  public:
    int             dfcn;               // degree of freedom at corner nodes
    int             dfmn;               // degree of freedom at midside nodes
    int             dfel;               // degree of freedom in elements

    int             maxEleq;            // maximum number of element equations
                                        // maxEleq = df * kMaxNodes2D
    double*         force;              // element force vector
    double**        estifm;             // element stiffness matrix

    int             modelInit;          // set to the counter "model->init" and
                                        // used to initialize the model structure
    EQS*            next;               // used in createEquation

    int             initStructure;      // flag to initialize the index matrix

    NODE**          eqnoNode;           // list of node pointers to determine the
                                        // corresponding node of an equation number
    int*            eqid;               // index of equation at node (0,1,2,...)

    int             neq;                // total number of equations
                                        // number of subdomain equations...
    int             neq_up;             // ...without any interface nodes
    int             neq_dn;             // ...without nodes on downstream interface

    int*            order;              // order of equations
    int*            end;                // end of equations
    int*            record;             // record numbers for equations

    int*            eqnoBuf;            // buffer for equation numbers
    REALPR*         eqBuf_L;            // buffer for elimination factors
    REALPR*         eqBuf_U;            // buffer for upper diagonal matrix

    double*         pivot;

    // -------------------------------- index matrices -----------------------------------
    CRSMAT*         crsm;               // CRS matrix

  public:
    // Eqs.cpp ---------------------------------------------------------------------------
    EQS( int dfcn, int dfmn, int dfel );
    virtual ~EQS();

    int          GetEqno( NODE*, int );
    int          GetEqno( ELEM*, int );

    void         ExportIM( const char* filename, CRSMAT* crsm, PROJECT* project );
    void         ExportEQS( const char* filename, CRSMAT* crsm, PROJECT* project );
    void         ExportVEC( const char* filename, double* vec, PROJECT* project );

    void         ScaleDiag( double* B, PROJECT* project );

    void         Mpi_assemble( double* vec, PROJECT* project );

    void         DataOut( char* name, int step, char* time, int release,
                          GRID* region, char* label[], ... );

    // IndexMat.cpp ----------------------------------------------------------------------
    void         KillCrsm();
    void         SetIndexMat( MODEL* model, int mceq );
    void         SortIndex( int, int*, int** );

    // Rotate.cpp ------------------------------------------------------------------------
    void         Rotate2D( int, NODE**, int, double**, double* );

    // SetEqno.cpp -----------------------------------------------------------------------
    void         SetEqno( MODEL*, int, int, int, unsigned int*, int );
    void         LastEquation( MODEL*, long );
    void         ResetEqOrder( MODEL* );

    // Solve.cpp -------------------------------------------------------------------------
    int          Solve( MODEL* model, int neq, double* rhs, double* x, PROJECT* project,
                        SOLVER* solver=NULL, PRECON** precon=NULL, int assemble=true );

    // Update.cpp ------------------------------------------------------------------------
    void         Update( MODEL*,SUBDOM*,double*,int,int,double*,double*,double*,double*,int*,int* );

    // -----------------------------------------------------------------------------------
    virtual int  Coefs( ELEM*, PROJECT*, double**, double* );
};

#endif
