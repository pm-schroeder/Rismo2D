// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_TM
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
#include "Vars.h"
#include "Shape.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsUVS2D_TM.h"


//#define kDebug


EQS_UVS2D_TM::EQS_UVS2D_TM() : EQS_UVS2D()
{
  neq      = 0;
  timegrad = 0;
}


EQS_UVS2D_TM::~EQS_UVS2D_TM()
{
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_TM::Timegrad( PROJECT* project, double dt, double th )
{
  // call the base class ---------------------------------------------------------------------------
  EQS_UVS2D::Timegrad( project, dt, th );
  return;

  MODEL*  model  = project->M2D;
  GRID*   rg     = model->region;
  DRYREW* dryRew = &rg->dryRew;
  int     np     = model->np;
  double  Qout   = -1.0;

  // set computation of time gradients on
  timegrad = 1;


  for( int n=0; n<rg->Getnp(); n++ )
  {
    NODE *nd = rg->Getnode(n);

    nd->v.dUdt = 0.0;
    nd->v.dVdt = 0.0;
    nd->v.dSdt = 0.0;
  }

  // convergence criterion for drying and rewetting ------------------------------------
  if( dryRew->method )
  {
    int wet;
    model->DoDryRewet( project, NULL, &wet );
  }

  // compute friction coefficients for each element ------------------------------------
  model->DoFriction( project );

  double* X = NULL;      // solution of equation system
  double* B = NULL;      // right hand side of equation system

  // initialize Reynolds stresses and eddy viscosity -----------------------------------
  rg->Turbulence( project );

  // set inflow and outflow condition --------------------------------------------------
  TIME bcTime = project->timeint.actualTime + project->timeint.incTime;

  project->timeint.actualBcSet->InitBcon( project, &bcTime, &Qout );

  model->SetLocation();
  model->SetNormal();
  model->SetRotation();

  // set equation numbers --------------------------------------------------------------
  project->fix[0] = BCON::kFixU;
  project->fix[1] = BCON::kFixV;
  project->fix[2] = BCON::kSetS;

  project->elemKind = ELEM::kRegion | ELEM::kBound;

  SetEqno( model, 3, 2, 0, project->fix, project->elemKind );

  initStructure = true;

  if( X ) MEMORY::memo.Detach( X );
  if( B ) MEMORY::memo.Detach( B );

  X = (double*) MEMORY::memo.Array_eq( neq );
  B = (double*) MEMORY::memo.Array_eq( neq );

  // solve equations -------------------------------------------------------------------
  for( int i=0; i<neq; i++ ) X[i] = 0.0;

  Solve( model, neq, B, X, project );

  // statistic of correction vector ----------------------------------------------------
  REPORT::rpt.Message( 2, "\n\n %s\n\n %s\n\n",
                          "(EQS_UVS2D_TM::Predict)    time gradients ...",
                          " variable   node           average          maximum" );

  for( int i=0; i<3; i++ )
  {
    int    varIndex = 0;
    char   varName  = ' ';
    int    noAbs,  noPer;
    double maxPer, avPer;
    double maxAbs, avAbs;

    switch( i )
    {
      case 0:  varName = 'U';  varIndex = kVarU;  break;
      case 1:  varName = 'V';  varIndex = kVarV;  break;
      case 2:  varName = 'h';  varIndex = kVarH;  break;
    }

    Update( model, &project->subdom,
            X, i, varIndex, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer );

    REPORT::rpt.Message( 2, "      %1c    %7d      %14.5le   %14.5le %s\n",
                            varName, noAbs, avAbs, maxAbs, "     (abs)" );

    REPORT::rpt.Message( 2, "      %1c    %7d      %14.5lf   %14.5lf %s\n\n",
                            ' ',     noPer, avPer, maxPer, "     ( % )" );
  }

  // apply correction ------------------------------------------------------------------
  for( int n=0; n<np; n++ )
  {
    NODE *nd = model->node[n];

    double dUdt, dVdt, dSdt;

    dUdt = dVdt = dSdt = 0.0;

    int eqnoU = GetEqno( nd, 0 );
    int eqnoV = GetEqno( nd, 1 );
    int eqnoS = GetEqno( nd, 2 );

    if( eqnoU >= 0 ) dUdt = X[eqnoU];
    if( eqnoV >= 0 ) dVdt = X[eqnoV];
    if( eqnoS >= 0 ) dSdt = X[eqnoS];

    if( isFS(nd->flag, NODE::kRotat) )
    {
      nd->v.dUdt  = nd->bc.Getrot(0,0) * dUdt;
      nd->v.dVdt  = nd->bc.Getrot(1,0) * dUdt;

      nd->v.dUdt += nd->bc.Getrot(0,1) * dVdt;
      nd->v.dVdt += nd->bc.Getrot(1,1) * dVdt;
    }

    else
    {
      nd->v.dUdt  = dUdt;
      nd->v.dVdt  = dVdt;
    }

    if( eqnoS >= 0 )
    {
      nd->v.dSdt  = dSdt;
    }
  }

  // compute S and time gradient at midside nodes (linear interpolation) -------------------------
  for( int e=0; e<rg->Getne(); e++ )
  {
    ELEM* el = rg->Getelem(e);

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    for( int i=ncn; i<nnd; i++ )
    {
      // get left and right corner node to midside node i ----------------------------------------
      int    il, ir;
      double left, rght;

      el->GetQShape()->getCornerNodes( i, &il, &ir );

      left = el->nd[il]->v.S;
      rght = el->nd[ir]->v.S;
      el->nd[i]->v.S = 0.5 * (left + rght);

      left = el->nd[il]->v.dSdt;
      rght = el->nd[ir]->v.dSdt;
      el->nd[i]->v.dSdt = 0.5 * (left + rght);
    }
  }

  MEMORY::memo.Detach( X );
  MEMORY::memo.Detach( B );


  // set computation of time gradients off
  timegrad = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_TM::Predict( PROJECT* project, int init, double dt, double th )
{
  // call the base class ---------------------------------------------------------------------------
  EQS_UVS2D::Predict( project, init, dt, th );
  return;

  MODEL*  model  = project->M2D;
  GRID*   rg     = model->region;
  DRYREW* dryRew = &rg->dryRew;
  int     np     = model->np;
  double  Qout   = -1.0;

  // set computation of time gradients on
  timegrad = 1;


  if( init )
  {
    for( int n=0; n<rg->Getnp(); n++ )
    {
      NODE *nd = rg->Getnode(n);

      nd->vo.dUdt = 0.0;
      nd->vo.dVdt = 0.0;
      nd->vo.dSdt = 0.0;
    }
  }

  else
  {
    // convergence criterion for drying and rewetting ------------------------------------
    if( dryRew->method )
    {
      int wet;
      model->DoDryRewet( project, NULL, &wet );
    }

    // compute friction coefficients for each element ------------------------------------
    model->DoFriction( project );

    double* X = NULL;      // solution of equation system
    double* B = NULL;      // right hand side of equation system

    // initialize Reynolds stresses and eddy viscosity -----------------------------------
    rg->Turbulence( project );

    // set inflow and outflow condition --------------------------------------------------
    TIME bcTime = project->timeint.actualTime + project->timeint.incTime;

    project->timeint.actualBcSet->InitBcon( project, &bcTime, &Qout );

    model->SetLocation();
    model->SetNormal();
    model->SetRotation();

    // set equation numbers --------------------------------------------------------------
    project->fix[0] = BCON::kFixU;
    project->fix[1] = BCON::kFixV;
    project->fix[2] = BCON::kSetS;

    project->elemKind = ELEM::kRegion | ELEM::kBound;

    SetEqno( model, 3, 2, 0, project->fix, project->elemKind );

    initStructure = true;

    if( X ) MEMORY::memo.Detach( X );
    if( B ) MEMORY::memo.Detach( B );

    X = (double*) MEMORY::memo.Array_eq( neq );
    B = (double*) MEMORY::memo.Array_eq( neq );

    // solve equations -------------------------------------------------------------------
    for( int i=0; i<neq; i++ )  X[i] = 0.0;

    Solve( model, neq, B, X, project );

    // statistic of correction vector ----------------------------------------------------
    REPORT::rpt.Message( 2, "\n\n %s\n\n %s\n\n",
                            "(EQS_UVS2D_TM::Predict)    time gradients ...",
                            " variable   node           average          maximum" );

    for( int i=0; i<3; i++ )
    {
      int    varIndex = 0;
      char   varName  = ' ';
      int    noAbs,  noPer;
      double maxPer, avPer;
      double maxAbs, avAbs;

      switch( i )
      {
        case 0:  varName = 'U';  varIndex = kVarU;  break;
        case 1:  varName = 'V';  varIndex = kVarV;  break;
        case 2:  varName = 'h';  varIndex = kVarH;  break;
      }

      Update( model, &project->subdom,
              X, i, varIndex, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer );

      REPORT::rpt.Message( 2, "      %1c    %7d      %14.5le   %14.5le %s\n",
                              varName, noAbs, avAbs, maxAbs, "     (abs)" );

      REPORT::rpt.Message( 2, "      %1c    %7d      %14.5lf   %14.5lf %s\n\n",
                              ' ',     noPer, avPer, maxPer, "     ( % )" );
    }


    // apply correction ------------------------------------------------------------------
    for( int n=0; n<np; n++ )
    {
      NODE *nd = model->node[n];

      double dUdt, dVdt, dSdt;

      dUdt = dVdt = dSdt = 0.0;

      int eqnoU = GetEqno( nd, 0 );
      int eqnoV = GetEqno( nd, 1 );
      int eqnoS = GetEqno( nd, 2 );

      if( eqnoU >= 0 ) dUdt = X[eqnoU];
      if( eqnoV >= 0 ) dVdt = X[eqnoV];
      if( eqnoS >= 0 ) dSdt = X[eqnoS];

      if( isFS(nd->flag, NODE::kRotat) )
      {
        nd->vo.dUdt  = nd->bc.Getrot(0,0) * dUdt;
        nd->vo.dVdt  = nd->bc.Getrot(1,0) * dUdt;

        nd->vo.dUdt += nd->bc.Getrot(0,1) * dVdt;
        nd->vo.dVdt += nd->bc.Getrot(1,1) * dVdt;
      }

      else
      {
        nd->vo.dUdt  = dUdt;
        nd->vo.dVdt  = dVdt;
      }

      if( eqnoS >= 0 )
      {
        nd->vo.dSdt  = dSdt;
      }
    }

    // ---------------------------------------------------------------------------------------------
    // initialize U,V,S and time gradients on previously dry nodes to their current values
    // NOTE: Especially the water elevation of the previously dry node has been initialized
    // with the bottom elevation and therefore the difference v.S - vo.S might be very large.

    double dryLimit = model->region->dryRew.dryLimit;

    for( int n=0; n<np; n++ )
    {
      NODE *nd = model->node[n];

      // check for previously dry nodes
      // if( isFS(nd->flag, NODE::kDryPrev)
      if( nd->vo.S <= nd->zor + dryLimit )
      {
        nd->vo.U = nd->v.U;
        nd->vo.V = nd->v.V;
        nd->vo.S = nd->v.S;

        nd->vo.dUdt = 0.0;
        nd->vo.dVdt = 0.0;
        nd->vo.dSdt = 0.0;
      }
    }

    // compute S and time gradient at midside nodes (linear interpolation) -------------------------
    for( int e=0; e<rg->Getne(); e++ )
    {
      ELEM* el = rg->Getelem(e);

      int ncn = el->Getncn();
      int nnd = el->Getnnd();

      for( int i=ncn; i<nnd; i++ )
      {
        // get left and right corner node to midside node i ----------------------------------------
        int    il, ir;
        double left, rght;

        el->GetQShape()->getCornerNodes( i, &il, &ir );

        left = el->nd[il]->vo.S;
        rght = el->nd[ir]->vo.S;
        el->nd[i]->vo.S = 0.5 * (left + rght);

        left = el->nd[il]->vo.dSdt;
        rght = el->nd[ir]->vo.dSdt;
        el->nd[i]->vo.dSdt = 0.5 * (left + rght);
      }
    }

    MEMORY::memo.Detach( X );
    MEMORY::memo.Detach( B );
  }

  // set computation of time gradients off
  timegrad = 0;
}
