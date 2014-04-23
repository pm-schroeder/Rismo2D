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

#include "EqsPPE2D.h"


EQS_PPE2D::EQS_PPE2D() : EQS( 3, 2, 0 )
{
  neq = 0;
}


EQS_PPE2D::~EQS_PPE2D()
{
}


// ---------------------------------------------------------------------------------------

void EQS_PPE2D::Execute( PROJECT* project )
{
  MODEL*  model  = project->M2D;
  GRID*   rg     = project->M2D->region;
  DRYREW* dryRew = &rg->dryRew;
  NODE**  node   = model->node;
  int     np     = model->np;

  int diverged_cg = 0;

  model->Incinit();


  // print information on actual iteration -----------------------------------------------
  project->PrintTheCycle( 1 );
  REPORT::rpt.PrintTime( 1 );

  // check for minimum flow depth --------------------------------------------------------
  for( int i=0; i<np; i++ )
  {
    double H;

    if( !isFS(node[i]->flag, NODE::kDry) )
    {
      H = node[i]->v.S - node[i]->z;

      if( H < dryRew->dryLimit )
      {
        if( dryRew->method == 2 )
        {
          node[i]->z = node[i]->v.S - dryRew->dryLimit;
        }

        else if( H <= 0.0 )
        {
          node[i]->z   = node[i]->zor;
          node[i]->v.S = node[i]->z + project->hmin;
        }
      }
    }
  }


  // compute friction coefficients -------------------------------------------------------

  model->DoFriction( project );


  // set inflow and outflow condition ----------------------------------------------------

  TIME bcTime = project->timeint.actualTime;

  project->timeint.actualBcSet->InitBcon( project, &bcTime );

  model->SetLocation();
  model->SetNormal();
  model->SetRotation();


  // -------------------------------------------------------------------------------------
  // set up equation numbers

  // set Dirichlet boundary conditions ---------------------------------------------------

  for( int i=0; i<np; i++ )
  {
    BCON *bc = &node[i]->bc;

    CF( bc->kind, BCON::kFix_1 | BCON::kFix_2 );

    if( isFS(bc->kind, BCON::kFixU) )  SF( bc->kind, BCON::kFix_1 );
    if( isFS(bc->kind, BCON::kFixV) )  SF( bc->kind, BCON::kFix_2 );

	if( isFS(bc->kind, BCON::kQInlet) )
    {
      SF( bc->kind, BCON::kFix_1 | BCON::kFix_2 );

      double H = node[i]->v.S - node[i]->z;

      if( H > 0.0 )
      {
        node[i]->v.U = bc->niox * bc->val->U / H;
        node[i]->v.V = bc->nioy * bc->val->U / H;
	  }
      else
      {
        node[i]->v.U = 0.0;
        node[i]->v.V = 0.0;
      }
    }

	else if( isFS(bc->kind, BCON::kInlet) )
    {
      SF( bc->kind, BCON::kFix_1 | BCON::kFix_2 );

      double H = node[i]->v.S - node[i]->z;

      if( H > 0.0 )
      {
        node[i]->v.U = -bc->val->U / H;
        node[i]->v.V = -bc->val->V / H;
      }
      else
      {
        node[i]->v.U = 0.0;
        node[i]->v.V = 0.0;
      }
    }
/*
    if( isFS(bc->kind, BCON::kOutlet) )
    {
      SF( node[i]->flag, NODE::kRotat );

      SF( bc->kind, BCON::kFix_2 );

      bc->rot[0][0] =  bc->niox;
      bc->rot[1][0] =  bc->nioy;

      bc->rot[0][1] =  bc->rot[1][0];   // tangential vector
      bc->rot[1][1] = -bc->rot[0][0];
    }
*/
  }


  project->fix[0] = BCON::kFix_1;
  project->fix[1] = BCON::kFix_2;
  project->fix[2] = 0;

  project->elemKind = ELEM::kRegion | ELEM::kBound;

  SetEqno( model, 3, 2, 0, project->fix, project->elemKind );

  double* B = (double*) MEMORY::memo.Array_eq( neq );
  double* X = (double*) MEMORY::memo.Array_eq( neq );


  for( int it=0; it<project->actualCycit; it++ )
  {
    // print information on actual iteration ---------------------------------------------

    if( it )
    {
      project->PrintTheCycle( it+1 );
      REPORT::rpt.PrintTime( 1 );
    }


    // -----------------------------------------------------------------------------------
    // solve for divergence free flow field

    for( int i=0; i<neq; i++ ) X[i] = 0.0;

    diverged_cg = Solve( model, neq, B, X, project );


    // statistics ------------------------------------------------------------------------

    int    noAbs, noPer;
    double maxAbs, maxPer, avAbs, avPer;

    REPORT::rpt.Message( 2, "\n\n %s\n\n %s\n\n",
                            "(EQS_PPE2D)     convergence parameters ...",
                            " variable   node          average          maximum" );


    Update( model, &project->subdom,
            X, 0, kVarU, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer);

    REPORT::rpt.Message( 2, "     %2s    %6d    %16.4le  %16.4le %s\n",
                            " U", noAbs, avAbs, maxAbs, "     (abs)" );

    REPORT::rpt.Message( 2, "           %6d    %16.4le  %16.4le %s\n\n",
                            noPer, avPer, maxPer, "     ( % )" );


    Update( model, &project->subdom,
            X, 1, kVarV, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer );

    REPORT::rpt.Message( 2, "     %2s    %6d    %16.4le  %16.4le %s\n",
                            " V", noAbs, avAbs, maxAbs, "     (abs)" );

    REPORT::rpt.Message( 2, "           %6d    %16.4le  %16.4le %s\n\n",
                            noPer, avPer, maxPer, "     ( % )" );


    // apply correction to velocities ----------------------------------------------------

    for( int i=0; i<np; i++ )
    {
      int eqno = GetEqno( node[i], 0 );

      if( eqno >= 0 )
      {
        BCON* bc = &node[i]->bc;

        if( isFS(node[i]->flag, NODE::kRotat) )
        {
          node[i]->v.U += node[i]->bc.Getrot(0,0) * X[eqno];
          node[i]->v.V += node[i]->bc.Getrot(1,0) * X[eqno];
        }

        else
        {
          node[i]->v.U += X[eqno];
        }
      }

      eqno = GetEqno( node[i], 1 );

      if( eqno >= 0 )
      {
        BCON* bc = &node[i]->bc;

        if( isFS(node[i]->flag, NODE::kRotat) )
        {
          node[i]->v.U += node[i]->bc.Getrot(0,1) * X[eqno];
          node[i]->v.V += node[i]->bc.Getrot(1,1) * X[eqno];
        }

        else
        {
          node[i]->v.V += X[eqno];
        }
      }
    }


    // -----------------------------------------------------------------------------------
    // determine discharge through continuity lines

    model->Continuity();


    // check for interrupt ---------------------------------------------------------------

    if( (diverged_cg & kErr_interrupt)  ||  it == project->actualCycit-1 )
    {
      REPORT::rpt.Message( 1, "\n\n%-25s%s: %d\n\n",
                              " (EQS_PPE2D::Execute)", "finished in iteration step", it+1 );

      REPORT::rpt.Message( 1, "\n\n%-25s%s\n\n %s\n\n",
                              " (EQS_PPE2D::Execute)", "convergence parameters ...",
                              " variable   node          average          maximum" );


      Update( model, &project->subdom,
              X, 0, kVarU, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer);

      REPORT::rpt.Message( 1, "     %2s    %6d    %16.4le  %16.4le %s\n",
                              " U", noAbs, avAbs, maxAbs, "     (abs)" );

      REPORT::rpt.Message( 1, "           %6d    %16.4le  %16.4le %s\n\n",
                              noPer, avPer, maxPer, "     ( % )" );


      Update( model, &project->subdom,
              X, 1, kVarV, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer );

      REPORT::rpt.Message( 1, "     %2s    %6d    %16.4le  %16.4le %s\n",
                              " V", noAbs, avAbs, maxAbs, "     (abs)" );

      REPORT::rpt.Message( 1, "           %6d    %16.4le  %16.4le %s\n\n",
                              noPer, avPer, maxPer, "     ( % )" );
      break;
    }
  }

  // compute friction coefficients -------------------------------------------------------

  model->DoFriction( project );


  // initialize Reynolds stresses and eddy viscosity -------------------------------------

  model->region->Turbulence( project );


  MEMORY::memo.Detach( B );
  MEMORY::memo.Detach( X );


  if( diverged_cg )  project->errLevel |= diverged_cg | kErr_no_conv_cg;
}
