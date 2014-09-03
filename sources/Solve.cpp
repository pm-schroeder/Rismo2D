// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
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
// /////////////////////////////////////////////////////////////////////////////////////////////////

#include "Defs.h"
#include "Report.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"
#include "Memory.h"
#include "CRSMat.h"
#include "Solver.h"
#include "Front.h"
#include "Frontm.h"
#include "Preco_ilu0.h"
#include "Preco_ilut.h"

#include "Eqs.h"

//#define kDebug


int EQS::Solve( MODEL*   model,
                int      neq,
                double*  B,
                double*  X,
                PROJECT* project,
                SOLVER*  slv,
                PRECON** precon,
                int      assemble )
{
  if( !slv )  slv = project->actualSolver;
  if( !slv )  return false;

  slv->eqs = this;

  double  scale = 0.0;
  double* rhs   = NULL;

  PRECON* pre = NULL;                   // a pointer to the preconditioner
  if( !precon )  precon = &pre;         // precon is the handle of the preconditioner

  int err = false;

  switch( slv->solverType )
  {
    // -----------------------------------------------------------------------------------

    case kFront:
#     ifdef _MPI_
      if( project->subdom.npr > 1 )
      {
        REPORT::rpt.Error( "frontal solver not supported in parallel- EQS::Solve(1)" );
      }
#     endif
      {
        FRONT* front = (FRONT*) slv;

        KillCrsm();

        front->OpenEquation();

        front->Direct( project, model, this, B );

        if( X )  for( int i=0; i<neq; i++ )  X[i] = B[i];

        front->CloseEquation();
      }
      break;

    // -----------------------------------------------------------------------------------

    case kFrontm:
      if( !X )  REPORT::rpt.Error( "solver not supported - EQS::Solve(2)" );

      if( this->initStructure )
      {
        KillCrsm();

        REPORT::rpt.Screen( 3, "\n ... setting index matrix\n" );

        SetIndexMat( model, slv->mceq );          // initialize index matrix EQS::crsm

        crsm->Alloc_A();                          // allocate memory for equation matrix

        this->initStructure = false;
      }


      // ---------------------------------------------------------------------------------
      // assemble equation system

      if( assemble )
      {
        crsm->Init();                             // initialize the matrix

        crsm->AssembleEqs_im( this, B, model, project );

        //////////////////////////////////////////////////////////////////////////////////
#       ifdef kDebug
        {
          char  filename[80];

          if( project->subdom.npr == 1 )
            sprintf( filename, "B.dbg" );
          else
            sprintf( filename, "B_%02d.dbg", project->subdom.pid+1 );

          ExportVEC( filename, B, project );

          if( project->subdom.npr == 1 )
            sprintf( filename, "D.dbg" );
          else
            sprintf( filename, "D_%02d.dbg", project->subdom.pid+1 );

          ExportVEC( filename, crsm->m_diag, project );

/*
          if( project->subdom.npr == 1 )
            sprintf( filename, "im.dbg" );
          else
            sprintf( filename, "im_%02d.dbg", project->subdom.pid );

          ExportIM( filename, crsm, project );


          if( project->subdom.npr == 1 )
            sprintf( filename, "eqs.dbg" );
          else
            sprintf( filename, "eqs_%02d.dbg", project->subdom.pid );

          ExportEQS( filename, crsm, project );
*/
        }
#       endif
        //////////////////////////////////////////////////////////////////////////////////
      }

      {
        FRONTM* frontm = (FRONTM*) slv;

        frontm->OpenEquation();
        frontm->Direct( project, crsm, B, X );
        frontm->CloseEquation();
      }
      break;

    // -----------------------------------------------------------------------------------

    case kBicgstab:
    case kParmsBcgstabd:
    case kParmsFgmresd:
      if( !X )  REPORT::rpt.Error( "solver not supported - EQS::Solve(3)" );

      if( this->initStructure )
      {
        KillCrsm();

        REPORT::rpt.Screen( 3, "\n ... setting index matrix\n" );

        SetIndexMat( model, slv->mceq );          // initialize index matrix EQS::crsm

        crsm->Alloc_A();                          // allocate memory for equation matrix

        this->initStructure = false;
      }


      // ---------------------------------------------------------------------------------
      // assemble equation system

      if( assemble )
      {
        crsm->Init();                             // initialize the matrix
        crsm->AssembleEqs_im( this, B, model, project );
      }

      ////////////////////////////////////////////////////////////////////////////////////
      // assemble local vectors from all adjacent subdomains
#     ifdef _MPI_
      rhs = (double*) MEMORY::memo.Array_eq( neq );
      memcpy( rhs, B, neq );
      Mpi_assemble( B, project );
#     endif
      ////////////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////////////
#     ifdef kDebug
      {
        char  filename[80];

        if( project->subdom.npr == 1 )
          sprintf( filename, "B.dbg" );
        else
          sprintf( filename, "B_%02d.dbg", project->subdom.pid+1 );

        ExportVEC( filename, B, project );

        if( project->subdom.npr == 1 )
          sprintf( filename, "D.dbg" );
        else
          sprintf( filename, "D_%02d.dbg", project->subdom.pid+1 );

        ExportVEC( filename, crsm->m_diag, project );
/*
        if( project->subdom.npr == 1 )
          sprintf( filename, "im.dbg" );
        else
          sprintf( filename, "im_%02d.dbg", project->subdom.pid );

        ExportIM( filename, crsm, project );

        if( project->subdom.npr == 1 )
          sprintf( filename, "eqs.dbg" );
        else
          sprintf( filename, "eqs_%02d.dbg", project->subdom.pid );

        ExportEQS( filename, crsm, project );
*/
      }
#     endif
      ////////////////////////////////////////////////////////////////////////////////////

      scale = crsm->ScaleL2Norm( B, &project->subdom );

      if( !(*precon) )
      {
        switch( slv->preconType )
        {
          default:
            *precon = NULL;
            break;

          case kPreco_ilu0:
            *precon = new PRECO_ILU0();

            // allocate memory for incomplete LU matrix ----------------------------------

            (*precon)->crsi = new CRSMAT( crsm );
            if( !(*precon)->crsi )
              REPORT::rpt.Error( kMemoryFault, "can not allocate memory - EQS::Solve(4)" );

            (*precon)->Factor( project, this, crsm );
            break;

          case kPreco_ilut:
            *precon = new PRECO_ILUT();

            // allocate memory for incomplete LU matrix ----------------------------------

            (*precon)->crsi = new CRSMAT( crsm );
            if( !(*precon)->crsi )
              REPORT::rpt.Error( kMemoryFault, "can not allocate memory - EQS::Solve(5)" );

            (*precon)->Factor( project, this, crsm );
            break;

          //case kPreco_ilut:
          //  *precon = new PRECO_ILUT();

          //  // allocate memory for incomplete LU matrix ----------------------------------

          //  (*precon)->crsi = new CRSMAT( crsm->m_neq, crsm->m_ceq, crsm->m_entries );
          //  if( !(*precon)->crsi )
          //    REPORT::rpt.Error( kMemoryFault, "can not allocate memory - EQS::Solve(6)" );

          //  (*precon)->Factor( project, this, crsm );
          //  break;
        }
      }

      // ---------------------------------------------------------------------------------

      if( !slv->Iterate( project, crsm, B, X, *precon ) )
      {
        err = true;

        // reset the right hand side to a not assembled vector ---------------------------
#       ifdef _MPI_
        memcpy( B, rhs, neq );
        for( int i=0; i<neq; i++ )  B[i] /= scale;
#       endif
      }

      if( rhs ) MEMORY::memo.Detach( rhs );

      iterCountCG += slv->iterCountCG;
      break;


    // -----------------------------------------------------------------------------------

    default:
      REPORT::rpt.Error( kParameterFault, "no valid solver specification - EQS::Solve(7)" );
  }

  if( err )
  {
    if( slv->proceed > 0 )
    {
      if( slv->accuracy >= 1.0 )  for( int i=0; i<neq; i++ ) X[i] = 0.0;

      SOLVER* next_slv = SOLVER::Getno(slv->proceed);
      next_slv->eqs = this;

      if( next_slv->maxDiff < slv->accuracy )
        err = Solve( model, neq, B, X, project, next_slv, precon, true );
      else
        err = false;
    }

    else if( slv->proceed == -1 ) // ||  slv->accuracy >= 1.0 )
    {
      REPORT::rpt.Message( 1, "\n" );
      REPORT::rpt.Warning( kSolverFault,
                  "Solver %d has not converged - iteration cancelled",
                  slv->solverType );

      for( int i=0; i<neq; i++ )  X[i] = 0.0;
      err = kErr_interrupt;
    }

    else if( slv->proceed == 0 )
    {
      REPORT::rpt.Message( 1, "\n" );
      REPORT::rpt.Warning( kSolverFault,
                  "Solver %d has not converged - iteration continued",
                  slv->solverType );
      err = kErr_some_errors;
    }
  }


  // release memory for preconditioner ---------------------------------------------------

  if( *precon )
  {
    delete *precon;
    *precon = NULL;
  }


  return err;
}
