// ======================================================================================
//                                     S O L V E R
// ======================================================================================
// This class implements the base class for solvers.
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

// read this header file only once ...

#ifndef SOLVER_INCL
#define SOLVER_INCL

#include "Defs.h"
#include "Report.h"

// ---------------------------------------------------------------------------------------
// equation solver

#define kFront                   1      // direct solvers
#define kFrontm                  2

#define kBicgstab                5

#define kParmsBcgstabd           6      // iterators from PARMS
#define kParmsFgmresd            7

#define kPreco_null              0      // preconditioner
#define kPreco_ilu0              1
#define kPreco_ilut              2


class EQS;
class CRSMAT;
class PRECON;
class MODEL;
class PROJECT;


class SOLVER
{
  ////////////////////////////////////////////////////////////////////////////////////////
  // members
  ////////////////////////////////////////////////////////////////////////////////////////

  public:
    static int      m_neqs;             // number of specified equation solvers
    static SOLVER** m_solver;

    // ----------------------------------------------------------------------------------
    // statical method Get() returns the number of user-defined solvers

  public:
    static int Get()
    {
      return m_neqs;
    }

    // ----------------------------------------------------------------------------------
    // statical method Getno() returns a pointer to the solver with number <no>

    static SOLVER* Getno( int no )
    {
      if( !m_solver )  return NULL;

      for( int i=0; i<m_neqs; i++ )
      {
        if( m_solver[i]->no == no )  return m_solver[i];
      }

      if( no == 0 )  return m_solver[0];

      REPORT::rpt.Error( kParameterFault,
                         "undefined solver type %d %s (%s)",
                         no, "used in time step file", "SOLVER::Getno #1" );
      return NULL;
    }

  public:
    int     no;

    int     solverType;            // solver type

    int     reinit;

    EQS*    eqs;                   // pointer to equation system

                                   // frontal solver -------------------------------------
    int     mfw;                   // maximum front width
    char    path[400];             // name for scratch file
    int     size;                  // size of equation buffer

                                   // cg and gmres solver --------------------------------
    int     preconType;            // preconditioner type

    int     mceq;                  // maximum number of connected equations

    int     proceed;               // how to procced if solver has failed

    int     maxIter;               // maximum of iterations in cg-solvers
    double  maxDiff;               // convergence criterion in cg-solvers

    int     mkyrl;                 // dimension of krylov subspace

    double  accuracy;              // accuracy of CG-Solver
    int     iterCountCG;           // total counter for CG iterations


  ////////////////////////////////////////////////////////////////////////////////////////
  // methods
  ////////////////////////////////////////////////////////////////////////////////////////

  public:
    SOLVER();
    virtual ~SOLVER();

    virtual int  Iterate( PROJECT* prj, CRSMAT* M, double* B, double* X, PRECON* P )
    { return false; };

    virtual void Direct( PROJECT* prj, MODEL* m, EQS* eqs, double* vec )
    { };

    virtual void Direct( PROJECT* prj, CRSMAT* M, double* rhs, double* x )
    { };
};

#endif
