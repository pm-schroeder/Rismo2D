// ======================================================================================
//                                        E Q S
// ======================================================================================
// This class implements the base class for differential equation systems.
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
