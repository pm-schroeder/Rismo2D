// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_ME
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

#include "EqsUVS2D_ME.h"

////////////////////////////////////////////////////////////////////////////////////////////////////

EQS_UVS2D_ME::EQS_UVS2D_ME() : EQS_UVS2D( 3, 0, 0 )
{
  UElimEq = NULL;
  VElimEq = NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

EQS_UVS2D_ME::~EQS_UVS2D_ME()
{
  if( UElimEq )
  {
    delete[] UElimEq[0];
    delete[] UElimEq;
  }
  if( VElimEq )
  {
    delete[] VElimEq[0];
    delete[] VElimEq;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_ME::Execute( PROJECT *project, int steadyFlow )
{
  MODEL*  model       = project->M2D;
  GRID*   rg          = model->region;
  DRYREW* dryRew      = &rg->dryRew;

  int     diverged_cg = 0;

  // allocate memory for elimination equations U,V of center node (bubble) -------------------------

  if( !modelInit )
  {
    int ne = rg->Getne();

    UElimEq = new double*[ne];
    VElimEq = new double*[ne];

    if( !UElimEq || !VElimEq )
      REPORT::rpt.Error( "can not allocate memory - EQS_UVS2D_ME::Execute(1)" );

    // maximum number of columns per element equation:
    //   4 corner nodes with 3 equations (U,V,S) on LHS = 12
    // + 1 center node  with 2 equations (U,V)   on LHS =  2
    // + RHS                                            =  1
    //                                              sum = 15
    UElimEq[0] = new double[15*ne];
    VElimEq[0] = new double[15*ne];

    if( !UElimEq[0] || !VElimEq[0] )
      REPORT::rpt.Error( "can not allocate memory - EQS_UVS2D_ME::Execute(2)" );

    for( int i=1; i<ne; i++ )
    {
      UElimEq[i] = UElimEq[i-1] + 15;
      VElimEq[i] = VElimEq[i-1] + 15;
    }
  }

  model->Incinit();

  double* X = NULL;      // solution of equation system
  double* B = NULL;      // right hand side of equation system

  // initialize the counter for cg-iterations --------------------------------------------
  iterCountCG = 0;

  // parameters for the time integration -------------------------------------------------
  double th = project->timeint.thetaFlow;
  double dt = project->timeint.incTime.Getsec();

  if( fabs(th) < 1.0e-10  ||  fabs(dt) < 1.0e-10 )
  {
    REPORT::rpt.Error( "theta or timeInterval too small (EQS_UVS2D::execute - 1)" );
  }

  // setup relyxytion parameter ----------------------------------------------------------
  double relaxDt_H  = project->timeint.relaxTimeFlow.Getsec();
  double relaxDt_UV = project->timeint.relaxTimeFlow.Getsec();

  double dt_H  = relaxDt_H;
  double dt_UV = relaxDt_UV;

  int relaxMethod = project->relaxMethod;

  if( steadyFlow )
  {
    if( relaxMethod >= 3 ) // stationary flow time relaxed computation
    {
      REPORT::rpt.Message( 5, "\n%-25s%s\n", " (EQS_UVS2D::Execute)",
                              "stationary flow / time relaxed iteration" );
      dt_H  = relaxDt_H;
      dt_UV = relaxDt_UV;

      relaxThdt_H  = 1.0 / dt_H;
      relaxThdt_UV = 1.0 / dt_UV;
    }
    else
    {
      REPORT::rpt.Message( 5, "\n%-25s%s\n", " (EQS_UVS2D::Execute)",
                              "stationary flow" );
      dt_H  = dt;
      dt_UV = dt;

      relaxThdt_H  = 0.0;
      relaxThdt_UV = 0.0;
    }
  }

  else // instationary flow
  {
    if( relaxMethod >= 3 ) // instationary flow time relaxed computation
    {
      REPORT::rpt.Message( 5, "\n%-25s%s\n", " (EQS_UVS2D::Execute)",
                              "instationary flow / time relaxed iteration" );
      dt_H  = relaxDt_H;
      dt_UV = relaxDt_UV;
    }
    else
    {
      REPORT::rpt.Message( 5, "\n%-25s%s\n", " (EQS_UVS2D::Execute)",
                              "instationary flow" );
      dt_H  = dt;
      dt_UV = dt;
    }

    relaxThdt_H  = 1.0 / dt_H  / th;
    relaxThdt_UV = 1.0 / dt_UV / th;
  }

  // write a note, if dispersion model is active -----------------------------------------
  if( project->actualDisp > 0 )
  {
    REPORT::rpt.Message( 4, "\n%-25s%s\n\n", " (EQS_UVS2D::Execute)",
                            "using dispersion coefficients from secondary flow model" );
  }

  // initialize countDown for dry-rewet-algorithm ----------------------------------------
  //for( int n=0; n<rg->Getnp(); n++ )  rg->Getnode(n)->countDown = 0;

  // print information on actual iteration -----------------------------------------------
  //project->PrintTheCycle( 0 );
  //REPORT::rpt.PrintTime( 1 );

  // do the predition step; determine time gradients -------------------------------------
  Predict( project, steadyFlow, dt, th );

  // -------------------------------------------------------------------------------------

  int     NRconv   = true;
  int     DWconv   = true;

  int     NRcontin = false;

  int     DWfreq   = dryRew->dryRewFreq;
  int     nextDW   = 0;

  int     eddy     = true;

  double  theNorm  = -1.0;
  double  relax    =  1.0;

  double  Qout     = -1.0;

  int     maxit    = project->actualCycit;

  double  last_dt  = dt_H;


  for( int it=0; it<maxit; it++ )
  {
    if( it == 0 || isFS(project->actualTurb, BCONSET::kVtIterat) )  eddy = true;

    // print information on actual iteration ---------------------------------------------
    if( REPORT::rpt.level > 1 )
    {
      project->PrintTheCycle( it+1 );
      REPORT::rpt.PrintTime( 2 );
    }

    // convergence criterion for drying and rewetting ------------------------------------
    if( dryRew->method  &&  DWfreq < project->actualCycit )
    {
      if( nextDW >= 0  &&  (it == nextDW || NRconv) )
      {
        nextDW = it + DWfreq;

        int wet;
        model->DoDryRewet( project, NULL, &wet );

        // prevent dry-rewet-convergence before time step is convergent
        if( wet  ||  dt_H < dt || dt_UV < dt )  DWconv = false;
        else                                    DWconv = true;

        eddy = true;

        if( nextDW > maxit )  // prevent further drying and rewetting
        {
          maxit  = nextDW;
          nextDW = -1;
          DWconv = true;
        }
      }
    }
/*
    else if( dryRew->method == 0 )
    {
      int any = false;

      // mark nodes in partially wetted elements with H <= 0 (marsh-nodes)
      for( int n=0; n<np; n++ )
      {
        NODE* nd = node[n];

        CF( nd->flag, NODE::kMarsh );
        CF( nd->flag, NODE::kDry );
        CF( nd->flag, NODE::kNoMoment );

        if( nd->v.S - nd->zor < dryRew->dryLimit )
        {
          any = true;
          SF( nd->flag, NODE::kDry );
        }
      }

      if( any )
      {
        for( int n=0; n<np; n++ )
        {
          CF( node[n]->flag, NODE::kNoMoment );
        }

        for( int e=0; e<rg->Getne(); e++ )
        {
          ELEM* el = rg->Getelem(e);

          // count number of dry nodes in element el
          int cnt = 0;
          int ncn = el->getncn();

          for( int i=0; i<ncn; i++ )
          {
            if( isFS(el->Getnode(i)->flag, NODE::kDry) )  cnt++;
          }

          // if not all nodes of the element are dry, set the kMarsh-flag on dry nodes
          if( cnt < ncn )
          {
            for( int i=0; i<ncn; i++ )
            {
              NODE* nd = el->Getnode(i);
              if( isFS(nd->flag, NODE::kDry) )  SF( nd->flag, NODE::kMarsh );
            }
          }
        }

        // set kNoMoment-condition for nodes with flag kDry set and kMarsh unset
        for( int n=0; n<np; n++ )
        {
          NODE* nd = node[n];

          if( isFS(nd->flag, NODE::kDry) )
          {
            SF( nd->flag, NODE::kNoMoment );
          }
        }

        // clear the kDry-flag for each node
        for( int n=0; n<np; n++ )
        {
          CF( node[n]->flag, NODE::kDry );
        }

        // reset boundary conditions for nodes with kNoMoment-flag
        model->SetRotation();
      }
    }
*/

    // exchange integer value DWconv
#   ifdef _MPI_
    DWconv = project->subdom.Mpi_max( DWconv );
#   endif

    // convergence criterion for Newton-Raphson ------------------------------------------
    NRconv = true;

    // compute friction coefficients for each element ------------------------------------
    model->DoFriction( project );

    // initialize Reynolds stresses and eddy viscosity -----------------------------------
    if( eddy )
    {
      eddy = false;
      rg->Turbulence( project );
    }

    // set inflow and outflow condition --------------------------------------------------
    TIME bcTime = project->timeint.actualTime + project->timeint.incTime;

    project->timeint.actualBcSet->InitBcon( project, &bcTime, &Qout, 1 );

    // update the iterated difference of the water level for gauge controlled nodes ------
    BCONSET *abcs = project->timeint.actualBcSet;

    if( project->ngct != abcs->ngct )
    {
      if( project->ngct > 0 ) delete[] project->gct;
      project->gct  = new BCVAL::GAUGECT[abcs->ngct];
      project->ngct = abcs->ngct;
    }

    for( int i=0; i<project->ngct; i++ ) project->gct[i] = abcs->gct[i];

    model->SetLocation();
    model->SetNormal();
    model->SetRotation();

    // set equation numbers --------------------------------------------------------------
    if( model->Getinit() != modelInit )
    {
      modelInit = model->Getinit();

      project->fix[0] = BCON::kFixU;
      project->fix[1] = BCON::kFixV;
      project->fix[2] = BCON::kSetS;

      project->elemKind = ELEM::kRegion | ELEM::kBound;

      SetEqno( model, 3, 0, 0, project->fix, project->elemKind );

      initStructure = true;

      if( X ) MEMORY::memo.Detach( X );
      if( B ) MEMORY::memo.Detach( B );

      X = (double*) MEMORY::memo.Array_eq( neq, "EQS_UVS2D_ME::Execute(1)" );
      B = (double*) MEMORY::memo.Array_eq( neq, "EQS_UVS2D_ME::Execute(2)" );
    }


    // solve equations -------------------------------------------------------------------

    for( int i=0; i<neq; i++ )  X[i] = 0.0;

    diverged_cg = Solve( model, neq, B, X, project );

    //////////////////////////////////////////////////////////////////////////////////////
#   ifdef kDebug
    {
      char  name[80];

      sprintf( name, "Eqs_UVS2D_%02d_%02d.report", it, project->subdom.pid+1 );
      FILE* id = fopen( name, "w" );

      fprintf( id, "Solution of equation system (U,V,S):\n" );
      fprintf( id, "%d\n", rg->Getnp() );

      for( int n=0; n<rg->Getnp(); n++ )
      {
        NODE* nd = rg->Getnode(n);

        int eqnoU = GetEqno( nd, 0 );
        int eqnoV = GetEqno( nd, 1 );
        int eqnoS = GetEqno( nd, 2 );

        if( eqnoU >= 0 )  fprintf( id, "%6d  %12.6lf | ", nd->Getname(), X[eqnoU] );
        else              fprintf( id, "%6d  %12s | ",    nd->Getname(), "************" );

        if( eqnoV >= 0 )  fprintf( id, "%12.6lf | ", X[eqnoV] );
        else              fprintf( id, "%12s | ",    "************" );

        if( eqnoS >= 0 )  fprintf( id, "%12.6lf\n", X[eqnoS] );
        else              fprintf( id, "%12s\n",    "************" );
      }

      fclose( id );
    }
#   endif
    //////////////////////////////////////////////////////////////////////////////////////


    // statistic of correction vector ----------------------------------------------------

    REPORT::rpt.Message( 2, "\n\n%-25s%s\n\n %s\n\n",
                            " (EQS_UVS2D::Execute)", "convergence parameters ...",
                            " variable   node           average          maximum" );

    int    noAbs[3], noPer[3];
    double mxPer[3], avPer[3];
    double mxAbs[3], avAbs[3];

    for( int i=0; i<3; i++ )
    {
      int    varIndex = 0;
      char   varName  = ' ';

      switch( i )
      {
        case 0:  varName = 'U';  varIndex = kVarU;  break;
        case 1:  varName = 'V';  varIndex = kVarV;  break;
        case 2:  varName = 'h';  varIndex = kVarH;  break;
      }

      Update( model, &project->subdom,
              X, i, varIndex, &mxAbs[i], &mxPer[i], &avAbs[i], &avPer[i], &noAbs[i], &noPer[i] );

      REPORT::rpt.Message( 2, "      %1c    %7d      %14.5le   %14.5le %s\n",
                              varName, noAbs[i], avAbs[i], mxAbs[i], "     (abs)" );

      REPORT::rpt.Message( 2, "      %1c    %7d      %14.5lf   %14.5lf %s\n\n",
                              ' ',     noPer[i], avPer[i], mxPer[i], "     ( % )" );


      switch( i )
      {
        case 0:  if( fabs(mxAbs[i]) > project->convUV )  NRconv = false;  break;
        case 1:  if( fabs(mxAbs[i]) > project->convUV )  NRconv = false;  break;
        case 2:  if( fabs(mxAbs[i]) > project->convS )   NRconv = false;  break;
      }
    }

    // exchange integer value NRconv (secure)
#   ifdef _MPI_
    NRconv = project->subdom.Mpi_max( NRconv );
#   endif

    // determine relaxation parameter for NEWTON-RAPHSON ---------------------------------

    double maxUs;
    double maxS;

    double relaxUs = 1.0;
    double relaxS  = 1.0;

    switch( relaxMethod )
    {
      default:
        REPORT::rpt.Warning( kParameterFault,
                             "relaxation method %d not supported", relaxMethod );

      case 0:
        theNorm = 0.0;

        for( int n=0; n<model->np; n++ )
        {
          NODE* nd = model->node[n];

          double dU = 0.0;
          double dV = 0.0;
          double dS = 0.0;

          int eqnoU = GetEqno( nd, 0 );
          int eqnoV = GetEqno( nd, 1 );
          int eqnoS = GetEqno( nd, 2 );

          if( eqnoU >= 0 )  dU = fabs( X[eqnoU] );
          if( eqnoV >= 0 )  dV = fabs( X[eqnoV] );
          if( eqnoS >= 0 )  dS = fabs( X[eqnoS] );

          theNorm += dU*dU + dV*dV + dS*dS;
        }

        theNorm = sqrt( theNorm );

        relax   = 1.0;
        relaxUs = relax;
        relaxS  = relax;

        REPORT::rpt.Message( 2, "\n%25s%s %12.4le\n%25s%s %12.4le\n",
                                " ", "relaxation:    norm    =", theNorm,
                                " ", "               relax   =", relax );
        break;

      case 2:
        theNorm = 0.0;
        maxUs   = 0.0;
        maxS    = 0.0;

        for( int n=0; n<model->np; n++ )
        {
          NODE* nd = model->node[n];

          double dU = 0.0;
          double dV = 0.0;
          double dS = 0.0;

          int eqnoU = GetEqno( nd, 0 );
          int eqnoV = GetEqno( nd, 1 );
          int eqnoS = GetEqno( nd, 2 );

          if( eqnoU >= 0 )  dU = fabs( X[eqnoU] );
          if( eqnoV >= 0 )  dV = fabs( X[eqnoV] );
          if( eqnoS >= 0 )  dS = fabs( X[eqnoS] );

          theNorm += dU*dU + dV*dV + dS*dS;

          double dUs = sqrt( dU*dU + dV*dV );

          if( dUs > maxUs )  maxUs = dUs;
          if( dS  > maxS  )  maxS  = dS;
        }

        //////////////////////////////////////////////////////////////////////////////////
        // MPI: broadcast max changes
#       ifdef _MPI_
        maxUs   = project->subdom.Mpi_max( maxUs );
        maxS    = project->subdom.Mpi_max( maxS );
        theNorm = project->subdom.Mpi_sum( theNorm );
#       endif
        //////////////////////////////////////////////////////////////////////////////////

        theNorm = sqrt( theNorm );

        relaxUs = project->maxDeltaUV / maxUs;
        relaxS  = project->maxDeltaS  / maxS;

        if( relaxUs > 1.0 )  relaxUs = 1.0;
        if( relaxS  > 1.0 )  relaxS  = 1.0;

        REPORT::rpt.Message( 2, "\n%25s%s %12.4le\n%25s%s %12.4le\n%25s%s %12.4le\n",
                                " ", "relaxation:    norm    =", theNorm,
                                " ", "               relaxUs =", relaxUs,
                                " ", "               relaxS  =", relaxS );
        break;

      case 3:
      case 4:
      {
        if( dt_H < dt  ||  dt_UV < dt )  NRconv = false;

        theNorm = 0.0;
        maxUs   = 0.0;
        maxS    = 0.0;

        for( int n=0; n<model->np; n++ )
        {
          NODE* nd = model->node[n];

          double dU = 0.0;
          double dV = 0.0;
          double dS = 0.0;

          int eqnoU = GetEqno( nd, 0 );
          int eqnoV = GetEqno( nd, 1 );
          int eqnoS = GetEqno( nd, 2 );

          if( eqnoU >= 0 )  dU = fabs( X[eqnoU] );
          if( eqnoV >= 0 )  dV = fabs( X[eqnoV] );
          if( eqnoS >= 0 )  dS = fabs( X[eqnoS] );

          theNorm += dU*dU + dV*dV + dS*dS;

          double dUs = sqrt( dU*dU + dV*dV );

          if( dUs > maxUs )  maxUs = dUs;
          if( dS  > maxS  )  maxS  = dS;
        }

        //////////////////////////////////////////////////////////////////////////////////
        // MPI: broadcast max changes
#       ifdef _MPI_
        maxUs   = project->subdom.Mpi_max( maxUs );
        maxS    = project->subdom.Mpi_max( maxS );
        theNorm = project->subdom.Mpi_sum( theNorm );
#       endif
        //////////////////////////////////////////////////////////////////////////////////

        theNorm = sqrt( theNorm );

        relaxUs = project->maxDeltaUV / maxUs;
        relaxS  = project->maxDeltaS  / maxS;

        REPORT::rpt.Message( 2, "\n%25s%s %12.4le\n%25s%s %12.4le\n",
                                " ", "relaxed time:  dt_Us   =", dt_UV,
                                " ", "               dt_S    =", dt_H  );

        last_dt = dt_H;

        dt_UV *= relaxUs;
        dt_H  *= relaxS;

        if( dt_UV > dt )          dt_UV = dt;
        if( dt_UV < relaxDt_UV )  dt_UV = relaxDt_UV;

        if( dt_H  > dt )          dt_H  = dt;
        if( dt_H  < relaxDt_H )   dt_H  = relaxDt_H;

//      if( !steadyFlow )
        {
          (dt_H < dt_UV)? (dt_UV = dt_H):(dt_H = dt_UV);
        }

        if( steadyFlow )
        {
          relaxThdt_UV = 1.0 / dt_UV;
          relaxThdt_H  = 1.0 / dt_H;
        }
        else
        {
          relaxThdt_UV = 1.0 / dt_UV / th;
          relaxThdt_H  = 1.0 / dt_H  / th;
        }

        if( relaxUs > 1.0 )  relaxUs = 1.0;
        if( relaxS  > 1.0 )  relaxS  = 1.0;

        REPORT::rpt.Message( 2, "\n%25s%s %12.4le\n%25s%s %12.4le\n%25s%s %12.4le\n",
                                " ", "relaxation:    norm    =", theNorm,
                                " ", "               relaxUs =", relaxUs,
                                " ", "               relaxS  =", relaxS );
        break;
      }
    }


    // apply correction ------------------------------------------------------------------

    int exUV = 0;
    int exS  = 0;

    for( int n=0; n<model->np; n++ )
    {
      NODE* nd = model->node[n];

      double U, V, Us, H, dU, dV, dS;

      dU = dV = dS = 0.0;

      int eqnoU = GetEqno( nd, 0 );
      int eqnoV = GetEqno( nd, 1 );
      int eqnoS = GetEqno( nd, 2 );

      if( eqnoU >= 0 ) dU = X[eqnoU] * relaxUs;
      if( eqnoV >= 0 ) dV = X[eqnoV] * relaxUs;
      if( eqnoS >= 0 ) dS = X[eqnoS] * relaxS;

      if( isFS(nd->flag, NODE::kRotat) )
      {
        nd->v.U += nd->bc.Getrot(0,0) * dU;
        nd->v.V += nd->bc.Getrot(1,0) * dU;

        nd->v.U += nd->bc.Getrot(0,1) * dV;
        nd->v.V += nd->bc.Getrot(1,1) * dV;
      }

      else
      {
        nd->v.U += dU;
        nd->v.V += dV;
      }

      U  = nd->v.U;
      V  = nd->v.V;
      Us = sqrt( U*U + V*V );

      if( Us > project->maxUs )
      {
        exUV++;

        if( exUV == 1 )
        {
          REPORT::rpt.Message( 4, "\n\n----------------------------------------" );
          REPORT::rpt.Message( 4, "------------------------------\n" );
          REPORT::rpt.Message( 4, "\n ... maximum velocity of %.2lf [m/s] exceeded at nodes:\n",
                               project->maxUs );
        }

        REPORT::rpt.Message( 4, "  %5d", nd->Getname() );

        if( exUV % 10 == 0 ) REPORT::rpt.Message( 4, "\n" );

        nd->v.U = U * project->maxUs / Us;
        nd->v.V = V * project->maxUs / Us;
      }

      if( eqnoS >= 0 )
      {
        nd->v.S += dS;

        switch( dryRew->method )
        {
          case 0:
          case 2:
          case 3:
            H = nd->v.S - nd->zor;

            if( H < dryRew->dryLimit )
            {
              exS++;
              nd->z = nd->v.S - dryRew->dryLimit;
            }
            else
            {
              nd->z = nd->zor;
            }
            break;

          default:
            H = nd->v.S - nd->z;

            if( H <= 0.0 )
            {
              exS++;
              nd->v.S = nd->z + project->hmin;
            }
            break;
        }
      }
    }

    if( exUV ) REPORT::rpt.Message( 4, "\n" );


    // solve velocities (U,V) for the center node (bubble) -----------------------------------------

    for( int e=0; e<rg->Getne(); e++ )
    {
      ELEM *el = rg->Getelem(e);

      if( !isFS(el->flag, ELEM::kDry) )
      {
        double* UElimPtr = UElimEq[el->Getno()];
        double* VElimPtr = VElimEq[el->Getno()];

        int ncn = el->Getncn();
        int nbn = el->GetBShape()->nnd;

        double fu = UElimPtr[2*nbn + ncn];        // RHS
        double fv = VElimPtr[2*nbn + ncn];

        for( int i=0; i<ncn; i++ )                // LHS
        {
          NODE *nd = el->nd[i];

          int eqnoU = GetEqno( nd, 0 );
          int eqnoV = GetEqno( nd, 1 );
          int eqnoS = GetEqno( nd, 2 );

          double Xu = (eqnoU >= 0)? (X[eqnoU]) : (0.0);
          double Xv = (eqnoV >= 0)? (X[eqnoV]) : (0.0);
          double Xs = (eqnoS >= 0)? (X[eqnoS]) : (0.0);

          fu -= UElimPtr[        i] * Xu;
          fu -= UElimPtr[nbn   + i] * Xv;
          fu -= UElimPtr[2*nbn + i] * Xs;

          fv -= VElimPtr[        i] * Xu;
          fv -= VElimPtr[nbn   + i] * Xv;
          fv -= VElimPtr[2*nbn + i] * Xs;
        }

        double du, dv;

        dv  = fv / VElimPtr[nbn + ncn];

        fu -= UElimPtr[nbn + ncn] * dv;
        du  = fu / UElimPtr[ncn];

        el->U += relaxUs * du;
        el->V += relaxUs * dv;
      }
    }


    if( exUV )
    {
      REPORT::rpt.Message( 2, "\n%-25s%s %d-times exceeded\n",
                              " (EQS_UVS2D::Execute)", "limit for maximum velocity", exUV );
    }

    if( exS )
    {
      REPORT::rpt.Message( 5, "\n%-25s%s %d nodes\n",
                              " (EQS_UVS2D::Execute)", "negative flow depth at", exS );
    }


    // compute velocities and flow depth at midside nodes (linear interpolation) -------------------

    for( int e=0; e<rg->Getne(); e++ )
    {
      ELEM *el = rg->Getelem(e);

      int ncn = el->Getncn();
      int nnd = el->Getnnd();

      for( int i=ncn; i<nnd; i++ )
      {
        // get left and right corner node to midside node i
        int lnd, rnd;
        el->GetQShape()->getCornerNodes( i, &lnd, &rnd );

        // interpolate water elevation
        double lS = el->nd[lnd]->v.S;
        double rS = el->nd[rnd]->v.S;
        el->nd[i]->v.S = 0.5 * (lS + rS);

//        // interpolate bottom elevation
//        double lz = el->nd[lnd]->z;
//        double rz = el->nd[rnd]->z;
//        el->nd[i]->z = 0.5 * (lz + rz);

//        double lH = lS - lz;
//        double rH = rS - rz;
//        double mH = el->nd[i]->v.S - el->nd[i]->z;

        // interpolate water elevation and compute flow depth
        double lU = el->nd[lnd]->v.U;
        double rU = el->nd[rnd]->v.U;
        el->nd[i]->v.U = 0.5 * (lU + rU);

        double lV = el->nd[lnd]->v.V;
        double rV = el->nd[rnd]->v.V;
        el->nd[i]->v.V = 0.5 * (lV + rV);
      }
    }


    // experimental: force flow through outlet boundary > 0 ----------------------------------------

    if( false && isFS(project->actualTurb, BCONSET::kVtLES) )
    {
      int    adapt = false;
      double maxUn = 0.001 * project->maxUs;

      for( int n=0; n<model->np; n++ )
      {
        NODE* nd = model->node[n];

        if( isFS(nd->bc.kind, BCON::kOutlet) )
        {
          // compute normal component of flow velocity on outlet node
          double Un = nd->v.U * nd->bc.niox  +  nd->v.V * nd->bc.nioy;

          if( Un < -maxUn )
          {
            adapt = true;

            double dUn = nd->v.U * (fabs(maxUn/Un) - 1.0);
            double dVn = nd->v.V * (fabs(maxUn/Un) - 1.0);

            nd->v.U += dUn;
            nd->v.V += dVn;

            nd->vo.U += dUn;
            nd->vo.V += dVn;

            REPORT::rpt.Message( 5, "\n%-25s%s %d %s (%le,%le)\n",
                                 " (EQS_UVS2D::Execute)",
                                 "recirculating flow on outlet node", nd->Getname(),
                                 "reduced by", dUn, dVn );
          }
        }
      }

      if( adapt )
      {
        REPORT::rpt.Message( 3, "\n%-25s%s\n",
                             " (EQS_UVS2D::Execute)",
                             "recirculating flow on outlet nodes reduced" );
      }
    }

    // report of Courant and Peclet number ---------------------------------------------------------
    if( REPORT::rpt.level > 1 ) rg->ReportCuPe( dt, project->vk );

    // compute time derivatives --------------------------------------------------------------------
    if( !steadyFlow ) Timegrad( project, dt, th );

    // determine discharge through continuity lines ------------------------------------------------
    if( REPORT::rpt.level > 1 ) model->Continuity();

    // check target water elevation at node (if set) -----------------------------------------------
    for( int i=0; i<project->ngauge; i++ )
    {
      if( project->gauge[i] > 0 )
      {
        NODE *ndg = NULL;

        if( project->subdom.npr > 1 ) ndg = project->subdom.node[project->gauge[i] - 1];
        else                          ndg = rg->Getnode(project->gauge[i] - 1);

        int    name = 0;
        double S    = 0.0;
        double So   = project->gaugeSo[i];

        if( ndg )
        {
          name = ndg->Getname();
          S    = ndg->v.S;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
#       ifdef _MPI_
        if( project->subdom.npr > 1 )
        {
//        MPI_Bcast( &name,  1, MPI_INT,    project->subdom.pid, MPI_COMM_WORLD );
//        MPI_Bcast( &S,     1, MPI_DOUBLE, project->subdom.pid, MPI_COMM_WORLD );

          if( project->subdom.pid == 0 ) // the master process (pid == 0) is receiving
          {
            MPI_Status status;

            for( int s=1; s<project->subdom.npr; s++ )
            {
              int    mpi_name;
              double mpi_S;

              MPI_Recv( &mpi_name, 1, MPI_INT,    s, 1, MPI_COMM_WORLD, &status );
              MPI_Recv( &mpi_S,    1, MPI_DOUBLE, s, 2, MPI_COMM_WORLD, &status );

              if( mpi_name > 0 )
              {
                name = mpi_name;
                S    = mpi_S;
              }
            }
          }
          else
          {
            MPI_Send( &name, 1, MPI_INT,    0, 1, MPI_COMM_WORLD );
            MPI_Send( &S,    1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD );
          }
        }

        if( project->subdom.pid == 0 ) // the master process (pid == 0) is sending
        {
          for( int s=1; s<project->subdom.npr; s++ )
          {
            MPI_Send( &name, 1, MPI_INT,    s, 1, MPI_COMM_WORLD );
            MPI_Send( &S,    1, MPI_DOUBLE, s, 2, MPI_COMM_WORLD );
          }
        }
        else
        {
          MPI_Status status;
          MPI_Recv( &name, 1, MPI_INT,    0, 1, MPI_COMM_WORLD, &status );
          MPI_Recv( &S,    1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status );
        }
#       endif
        ////////////////////////////////////////////////////////////////////////////////////////////

        REPORT::rpt.Message( 1, "\n\n%-25s%s%d\n%-25s%s%.5lf%s%.5lf\n\n",
                                " (EQS_UVS2D::Execute)", "water elevation at node ",
                                name,
                                " ", "S = ", S, "  |  target S = ", So );
      }
    }

    // check convergence ---------------------------------------------------------------------------
    if( (NRconv && DWconv)  ||  (diverged_cg & kErr_interrupt)  ||  it == project->actualCycit-1 )
    {
      theNorm = -1.0;

      if( REPORT::rpt.level == 1 )
      {
        project->PrintTheCycle( it+1 );
        REPORT::rpt.PrintTime( 1 );

        REPORT::rpt.Message( 1, "\n\n%-25s%s\n\n %s\n\n",
                                " (EQS_UVS2D::Execute)", "convergence parameters ...",
                                " variable   node           average          maximum" );
        for( int i=0; i<3; i++ )
        {
          char   varName  = ' ';

          switch( i )
          {
            case 0:  varName = 'U';  break;
            case 1:  varName = 'V';  break;
            case 2:  varName = 'h';  break;
          }

          REPORT::rpt.Message( 1, "      %1c    %7d      %14.5le   %14.5le %s\n",
                                  varName, noAbs[i], avAbs[i], mxAbs[i], "     (abs)" );

          REPORT::rpt.Message( 1, "      %1c    %7d      %14.5lf   %14.5lf %s\n\n",
                                  ' ',     noPer[i], avPer[i], mxPer[i], "     ( % )" );
        }

        rg->ReportCuPe( dt, project->vk );
        model->Continuity();
      }
      break;
    }

    // or reset time increment to length of last relaxed time step -----------------------
    else if( it == maxit-1  &&  !steadyFlow  &&  !NRconv  &&  relaxMethod >= 3 )
    {
      REPORT::rpt.Line1( 2 );
      REPORT::rpt.PrintTime( 2 );

      if( 1.0-last_dt/dt > 1.0e-6 )     // length of specified time not gained -----------
      {
//        if( last_dt < dt_H + 1.0e-6 )
//        {
//          REPORT::rpt.Message( 2, "\n%-25s%s  %.4le\n",
//                                  " (EQS_UVS2D::Execute)",
//                                  "continuing cycle to gain specified time", dt );
//          NRcontin = true;
//          maxit++;
//        }
//        else
        {
          REPORT::rpt.Message( 2, "\n%-25s%s  %.4le\n",
                                  " (EQS_UVS2D::Execute)",
                                  "restarting cycle with reduced time", last_dt );

          project->timeint.incTime.Setsec( last_dt );

          dt_H = dt_UV = dt = project->timeint.incTime.Getsec();

          relaxThdt_H  = 1.0 / dt_H  / th;
          relaxThdt_UV = 1.0 / dt_UV / th;

          maxit += project->actualCycit;

          NRcontin = false;
        }
      }

      else if( NRcontin )     // not yet convergent: try some further iterations ---------
      {
        REPORT::rpt.Message( 2, "\n%-25s%s\n",
                                " (EQS_UVS2D::Execute)",
                                "continuing cycle to gain convergence" );
        NRcontin = false;
        maxit += project->actualCycit/4;
      }
    }
  }

  MEMORY::memo.Detach( X );
  MEMORY::memo.Detach( B );


  // experimental: force flow through outlet boundary > 0 --------------------------------

  if( isFS(project->actualTurb, BCONSET::kVtLES) )
  {
    int    adapt = false;
    double maxUn = 0.001 * project->maxUs;

    for( int n=0; n<model->np; n++ )
    {
      NODE* nd = model->node[n];

      if( isFS(nd->bc.kind, BCON::kOutlet) )
      {
        // compute normal component of flow velocity on outlet node
        double Un = nd->v.U * nd->bc.niox  +  nd->v.V * nd->bc.nioy;

        if( Un < -maxUn )
        {
          adapt = true;

          double dUn = nd->v.U * (fabs(maxUn/Un) - 1.0);
          double dVn = nd->v.V * (fabs(maxUn/Un) - 1.0);

          nd->v.U += dUn;
          nd->v.V += dVn;

          nd->vo.U += dUn;
          nd->vo.V += dVn;

          REPORT::rpt.Message( 5, "\n%-25s%s %d %s (%le,%le)\n",
                               " (EQS_UVS2D::Execute)",
                               "recirculating flow on outlet node", nd->Getname(),
                               "reduced by", dUn, dVn );
        }
      }
    }

    if( adapt )
    {
      REPORT::rpt.Message( 3, "\n%-25s%s\n",
                           " (EQS_UVS2D::Execute)",
                           "recirculating flow on outlet nodes reduced" );
    }
  }

  // report time gradients ---------------------------------------------------------------
  if( !steadyFlow )
  {
    double maxU, maxV, maxS;
    double aveU, aveV, aveS;

    maxU = aveU = 0.0;
    maxV = aveV = 0.0;
    maxS = aveS = 0.0;

    int j = 0;

    for( int n=0; n<model->np; n++ )
    {
      NODE* nd = model->node[n];

      if( !isFS(nd->flag, NODE::kDry) )
      {
        if( fabs(nd->v.dUdt) > fabs(maxU) )  maxU = nd->v.dUdt;
        if( fabs(nd->v.dVdt) > fabs(maxV) )  maxV = nd->v.dVdt;
        if( fabs(nd->v.dSdt) > fabs(maxS) )  maxS = nd->v.dSdt;

        aveU += fabs(nd->v.dUdt);
        aveV += fabs(nd->v.dVdt);
        aveS += fabs(nd->v.dSdt);

        j++;
      }
    }

    if( j )
    {
      aveU /= j;
      aveV /= j;
      aveS /= j;


      REPORT::rpt.Message( 1, "\n%-25s%s\n\n%25s%s\n\n",
                              " (EQS_UVS2D::Execute)", "time gradients ...",
                              " ", "variable      average        maximum" );

      REPORT::rpt.Message( 1, "%25s    %1c     %12.5le   %12.5le\n",
                              " ", 'U', aveU, maxU );

      REPORT::rpt.Message( 1, "%25s    %1c     %12.5le   %12.5le\n",
                              " ", 'V', aveV, maxV );

      REPORT::rpt.Message( 1, "%25s    %1c     %12.5le   %12.5le\n",
                              " ", 'S', aveS, maxS );
    }
  }


  rg->SetSlipFlow();


  if( iterCountCG )
  {
    REPORT::rpt.Message( 5, "\n%-25s%s %d\n",
                            " (EQS_UVS2D::Execute)", "total number of cg-iterations", iterCountCG );
  }

  if( !NRconv )
  {
    REPORT::rpt.Message( 0, "\n" );
    REPORT::rpt.Warning( kSolverFault, "time step %d: %s\n",
                         project->GetTimeStep(), "Newton Raphson has not converged" );
  }

  if( !NRconv )      project->errLevel |= kErr_some_errors | kErr_no_conv_nr;
  if( diverged_cg )  project->errLevel |= diverged_cg | kErr_no_conv_cg;
}
