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
#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsSL2D.h"


EQS_SL2D::EQS_SL2D() : EQS( 1, 1, 0 )
{
  neq  = 0;
}


EQS_SL2D::~EQS_SL2D()
{
}

// ---------------------------------------------------------------------------------------

void EQS_SL2D::Execute( PROJECT* project )
{
  MODEL* model       = project->M2D;
  GRID*  rg          = project->M2D->region;
  NODE** node        = model->node;
  int    np          = model->np;
  int    diverged_cg = 0;

  int    i;
  char   text[200];

  model->Incinit();

  // print information on actual iteration -----------------------------------------------
  project->PrintTheCycle( 1 );
  REPORT::rpt.PrintTime( 1 );

  // compute friction coefficients -------------------------------------------------------
  model->DoFriction( project );

  // initialize Reynolds stresses and eddy viscosity -------------------------------------
  rg->Turbulence( project );

  // initialize turbulent diffusivity/dispersion -----------------------------------------
  rg->EddyDisp();

  // set up equation numbers -------------------------------------------------------------
  project->fix[0]   = BCON::kSetC;
  project->elemKind = ELEM::kRegion | ELEM::kBound;

  SetEqno( model, 1, 1, 0, project->fix, project->elemKind );

  double* B = (double*) MEMORY::memo.Array_eq( neq );
  double* X = (double*) MEMORY::memo.Array_eq( neq );

  if( model->Getinit() != modelInit )
  {
    initStructure = true;
    modelInit = model->Getinit();
  }


  // -------------------------------------------------------------------------------------

  for( int it=0; it<project->actualCycit; it++ )
  {
    // print information on actual iteration -----------------------------------------------
    project->PrintTheCycle( it+1 );
    REPORT::rpt.PrintTime( 2 );

    for( i=0; i<neq; i++ )  X[i] = 0.0;


    // solve equations -------------------------------------------------------------------
    diverged_cg = Solve( model, neq, B, X, project );


    // update node variables with equation solution and determine convergence parameters -

    // statistics of correction vector ---------------------------------------------------
    int    noAbs, noPer;
    double maxAbs, maxPer, avAbs, avPer;

    sprintf( text, "\n\n%-25s%s\n\n %s\n\n",
                   " (EQS_SL2D::Execute)", "convergence parameters ...",
                   " variable   node          average        maximum" );
    REPORT::rpt.Message( 2, text );

    Update( model, &project->subdom,
            X, 0, kVarC, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer );

    sprintf( text, "      %2s    %5d       %12.5le   %12.5le %s\n",
                   "qs", noAbs, avAbs, maxAbs, "     (abs)" );
    REPORT::rpt.Message( 2, text );

    sprintf( text, "      %2s    %5d       %12.5lf   %12.5lf %s\n\n",
                   "  ", noPer, avPer, maxPer, "     ( % )" );
    REPORT::rpt.Message( 2, text );

    // apply correction ------------------------------------------------------------------
    for( int n=0; n<np; n++ )
    {
      int eqno = GetEqno( node[n], 0 );

      if( eqno >= 0 )  node[n]->v.C += X[eqno];

      if( node[n]->v.C < 0.0 )  node[n]->v.C = 0.0;
    }

    model->ContinuityBSL();
  }

  rg->ReportCuPe( project->timeint.incTime.Getsec(), project->vk );

  MEMORY::memo.Detach( X );
  MEMORY::memo.Detach( B );

  if( diverged_cg )  project->errLevel |= diverged_cg | kErr_no_conv_cg;
}
