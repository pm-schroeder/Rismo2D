// ======================================================================================
//
// Copyright (C) 1992-2010  by  P.M. SCHROEDER
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
#include "Type.h"
#include "Vars.h"
#include "Shape.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsD2D.h"


EQS_D2D::EQS_D2D() : EQS( 1, 1, 0 )
{
  neq  = 0;
}


EQS_D2D::~EQS_D2D()
{
}


// ---------------------------------------------------------------------------------------

void EQS_D2D::Execute( PROJECT* project, int steadyFlow, int linearShape )
{
  MODEL* model = project->M2D;
  GRID*  rg    = project->M2D->region;

  NODE** node  = model->node;
  ELEM** elem  = model->elem;

  int np       = model->np;
  int ne       = model->ne;

  int div_cg   = 0;

  double* B    = NULL;
  double* X    = NULL;

  model->Incinit();

  // print information on actual iteration -----------------------------------------------
  project->PrintTheCycle( 1 );
  REPORT::rpt.PrintTime( 1 );

  // set parameters according to time integration and relaxation -------------------------
  double th = project->timeint.thetaTurb;
  double dt = project->timeint.incTime.Getsec();

  if( fabs(th) < 1.0e-10  &&  fabs(dt) < 1.0e-10 )
  {
    REPORT::rpt.Error( kParameterFault, "theta and timeInterval too small (EQS_D2D::execute - 1)" );
  }

  double thdt = 1.0 / dt / th;

  double dt_KD;
  double relaxDt_KD = project->timeint.relaxTimeTurb.Getsec();

  int relaxMethod = project->relaxMethod;

  if( steadyFlow )
  {
    if( relaxMethod >= 3 )
    {
      dt_KD = relaxDt_KD;
      relaxThdt_KD = 1.0 / dt_KD;
    }
    else
    {
      dt_KD = dt;
      relaxThdt_KD = 0.0;
    }

    for( int i=0; i<np; i++ )
    {
      node[i]->v.dDdt = 0.0;
    }
  }
  else
  {
    if( relaxMethod >= 3 )
    {
      dt_KD = relaxDt_KD;
    }
    else
    {
      dt_KD = dt;
    }

    relaxThdt_KD = 1.0 / dt_KD / th;

    // time prediction -------------------------------------------------------------------
    for( int i=0; i<np; i++ )
    {
      VARS* v = &node[i]->v;

      //VARS* vo = &(node[i]->vo);

      //v->D    = vo->D  +  vo->dDdt * dt;
      //v->dDdt = (1.0 - 1.0/th)*vo->dDdt + thdt*(v->D - vo->D);

      v->dDdt = 0.0;
    }
  }

  // check KD-values for validity (>= 0) -------------------------------------------------
  Validate( np, node, project );

  // determine friction coefficients -----------------------------------------------------
  model->DoFriction( project );

  // initialize Reynolds stresses and eddy viscosity -------------------------------------
  rg->Turbulence( project );

  // set KD boundary conditions ----------------------------------------------------------
  model->SetBoundKD( project );


  // -------------------------------------------------------------------------------------

  int conv = true;

  for( int it=0; it<project->actualCycit; it++ )
  {
    conv = true;

    // print information on actual iteration ---------------------------------------------
    if( it > 0 )
    {
      project->PrintTheCycle( it+1 );
      REPORT::rpt.PrintTime( 1 );
    }

    // initialize Reynolds stresses and eddy viscosity -----------------------------------
    if( it > 0  &&  isFS(project->actualTurb, BCONSET::kVtIterat)
                && !isFS(project->actualTurb, BCONSET::kVtPrandtlKol) )
    {
      rg->Turbulence( project );
    }

    // set up equation numbers -----------------------------------------------------------
    if( model->Getinit() != modelInit )
    {
      initStructure = true;
      modelInit = model->Getinit();

      project->fix[0] = BCON::kFixD;

      project->elemKind = ELEM::kRegion;

      
      SetEqno( model, 1, 1, 0, project->fix, project->elemKind );

      if( B )  MEMORY::memo.Detach( B );
      if( X )  MEMORY::memo.Detach( X );

      B = (double*) MEMORY::memo.Array_eq( neq );
      X = (double*) MEMORY::memo.Array_eq( neq );
    }

    // solve equations with frontal solving algorithm ------------------------------------
    for( int i=0; i<neq; i++ ) X[i] = 0.0;

    div_cg = Solve( model, neq, B, X, project );


    // -----------------------------------------------------------------------------------
    // statistics

    REPORT::rpt.Message( 1, "\n\n%-25s%s\n\n %s\n\n",
                            " (EQS_D2D::Execute)","convergence parameters ...",
                            " variable   node           average          maximum" );


    // compute averaged changes and standard deviation of D ------------------------------
    double stdevD  = 0.0;
    double aveAbsD = 0.0;

    for( int i=0; i<np; i++ )
    {
      int eqno;
      double dD = 0.0;

      eqno = GetEqno( node[i], 0 );
      if( eqno >= 0 )  dD = X[eqno];

      aveAbsD += dD;
      stdevD  += dD*dD;
    }

    int nptot = 0;

    for( int i=0; i<np; i++ )
    {
      NODE* nd = rg->Getnode(i);
      if( !isFS(nd->flag, NODE::kInface_DN) )  nptot++;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: broadcast statistic
#   ifdef _MPI_
    aveAbsD = project->subdom.Mpi_sum( aveAbsD );
    stdevD  = project->subdom.Mpi_sum( stdevD );

    nptot   = project->subdom.Mpi_sum( nptot );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    aveAbsD /= nptot;
    stdevD   = sqrt( stdevD / nptot );

    double norm   = stdevD;
    double fractD = 2.0 * sqrt( stdevD );


    // compute maximum changes of K and D limited to ~95% fractile -----------------------

    int    cntD    = 0;

    int    noPerD  = 0;
    double avePerD = 0.0;
    double maxPerD = 0.0;

    int    noAbsD  = 0;
    double maxAbsD = 0.0;

    for( int i=0; i<np; i++ )
    {
      int eqno;
      double dK = 0.0;
      double dD = 0.0;

      eqno = GetEqno( node[i], 0 );
      if( eqno >= 0 )
      {
        dD = X[eqno];
        if( fabs(dD) > project->maxDeltaKD  &&  fabs(dD) > fractD )
        {
          X[eqno] = dD/fabs(dD) * fractD;
        }
      }

      if( node[i]->v.D + dD > 0.0 )
      {
        if( fabs(dD) > fabs(maxAbsD) )
        {
          maxAbsD = dD;
          noAbsD = node[i]->Getname();
        }

        if( node[i]->v.D > 0.0 )
        {
          double per = dD / node[i]->v.D;

          avePerD += fabs(per);
          cntD++;

          if( fabs(per) > fabs(maxPerD) )
          {
            maxPerD = per;
            noPerD = node[i]->Getname();
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: broadcast statistic
#   ifdef _MPI_
    cntD    = project->subdom.Mpi_sum( cntD );
    maxAbsD = project->subdom.Mpi_max( maxAbsD );
    avePerD = project->subdom.Mpi_sum( avePerD );
    maxPerD = project->subdom.Mpi_max( maxPerD );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    avePerD /= cntD;

    REPORT::rpt.Message( 1, "      %1c     %5d       %12.5le   %12.5le %s\n",
                            'D', noAbsD+1, aveAbsD, maxAbsD, "     (abs)" );

    REPORT::rpt.Message( 1, "      %1c     %5d       %12.5lf   %12.5lf %s\n\n",
                            ' ', noPerD+1, avePerD, maxPerD, "     ( % )" );

    if( fabs(maxAbsD) > project->convKD )  conv = false;


    // determine relaxation parameter for NEWTON-RAPHSON ---------------------------------

    double relax;

    double maxKD = fabs(maxAbsD);

    switch( relaxMethod )
    {
      default:
        REPORT::rpt.Warning( kParameterFault,
                             "relaxation method %d not supported", relaxMethod );

      case 0:
        relax = 1.0;

        REPORT::rpt.Message( 1, "\n%-25s%s %12.4le\n %s %12.4le\n\n",
                                " ", "relaxation:   norm  =", norm,
                                " ", "              relax =", relax );
        break;

      case 2:
        relax = project->maxDeltaKD / maxKD;

        if( relax > 1.0 )  relax = 1.0;

        REPORT::rpt.Message( 1, "\n%-25s%s %12.4le\n%-25s%s %12.4le\n\n",
                                " ", "relaxation:   norm  =", norm,
                                " ", "              relax =", relax );
        break;

      case 3:
      case 4:
        REPORT::rpt.Message( 1, "\n%-25s%s %12.4le\n",
                                " ", "relaxed time: dt_KD =", dt_KD );

        if( dt_KD < dt )  conv = false;

        relax = project->maxDeltaKD / maxKD;

        dt_KD *= relax;

        if( dt_KD > dt )          dt_KD = dt;
        if( dt_KD < relaxDt_KD )  dt_KD = relaxDt_KD;

        if( steadyFlow ) relaxThdt_KD = 1.0 / dt_KD;
        else             relaxThdt_KD = 1.0 / dt_KD / th;

        if( relax > 1.0 )  relax = 1.0;

        REPORT::rpt.Message( 1, "\n%-25s%s %12.4le\n%-25s%s %12.4le\n\n",
                                " ", "relaxation:   norm  =", norm,
                                " ", "              relax =", relax );
        break;
    }


    // update ----------------------------------------------------------------------------

    for( int i=0; i<np; i++ )
    {
      int n = GetEqno( node[i], 0 );
      if( n >= 0 )  node[i]->v.D += relax * X[n];
    }


    // check for range of values ---------------------------------------------------------

    double minD = 0.0;
    double maxD = 0.0;

    int first = true;

    int jd = 0;

    for( int i=0; i<np; i++ )
    {
      if( first )
      {
        first = false;
        minD  =
        maxD  = node[i]->v.D;
      }
      else
      {
        if( node[i]->v.D < minD )  minD = node[i]->v.D;
        if( node[i]->v.D > maxD )  maxD = node[i]->v.D;
      }

      if( node[i]->v.D <= 0.0 )
      {
//        node[i]->v.D = project->minD;
        if( GetEqno(node[i],1) >= 0 )  jd++;
      }
    }


    // compute useful values for K and D where they are negative -------------------------
    if( jd )
    {
      REPORT::rpt.Message( 1, "%-25s%d %s\n\n", " ", jd, "nodes with D out of range" );
      Validate( np, node, project );
    }

    REPORT::rpt.Message( 1, "%-25sminimum of D: %le\n", " ", minD );
    REPORT::rpt.Message( 1, "%-25smaximum of D: %le\n", " ", maxD );


    // compute time derivatives ----------------------------------------------------------

    rg->ReportCuPe( dt, project->vk );


    if( !steadyFlow )
    {
      double iTheta;

      iTheta  = 1.0  -  1.0 / project->timeint.thetaTurb;

      for( int i=0; i<np; i++ )
      {
        if(     isFS(node[i]->flag, NODE::kDry)
            ||  isFS(node[i]->flag, NODE::kMarsh) )
        {
          node[i]->v.dDdt = 0.0;
        }

        else
        {
          double D, pD, pdDdt;
          VARS *v, *vo;

          v  = &(node[i]->v);
          vo = &(node[i]->vo);

          D     = v->D;
          pD    = vo->D;
          pdDdt = vo->dDdt;

          // compute derivatives at actual time step
          node[i]->v.dDdt = iTheta*pdDdt + thdt*(D - pD);
        }
      }
    }


    if( conv )  break;
  }


  // finally: compute eddy viscosity from revised turbulence parameters K,D --------------

  rg->Turbulence( project );

  // -------------------------------------------------------------------------------------


  MEMORY::memo.Detach( B );
  MEMORY::memo.Detach( X );


  if( !conv )   project->errLevel |= kErr_some_errors | kErr_no_conv_nr;
  if( div_cg )  project->errLevel |= div_cg | kErr_no_conv_cg;
}


//////////////////////////////////////////////////////////////////////////////////////////
// check D-values for validity (>= 0)
//////////////////////////////////////////////////////////////////////////////////////////

void EQS_D2D::Validate( int np, NODE** node, PROJECT* project )
{
//  for( int iter=0; iter<1000; iter++ )
//  {
//    int further = false;
//
//    for( int i=0; i<np; i++ )
//    {
//      NODE* nd = node[i];
//
//      if( isFS(nd->flag, NODE::kDry) )  continue;
//
//      double U  = nd->v.U;
//      double V  = nd->v.V;
//      double Us = sqrt( nd->cf * (U*U + V*V) );
//
//      double H  = nd->v.S - nd->z;
//      if( H <= 0.0 ) H = project->hmin;
//
//      if( nd->v.D <= 0.0 )
//      {
//        further = true;
//
//        int noel = nd->noel;
//
//        double Dave = 0.0;
//        int    cntD = 0;
//
//        for( int j=0; j<noel; j++ )
//        {
//          ELEM* el = nd->el[j];
//
//          for( int k=0; k<el->getnnd(); k++ )
//          {
//            double D = el->nd[k]->v.D;
//            if( D > 0.0 )
//            {
//              Dave += D;
//              cntD++;
//            }
//          }
//        }
//
//        // ***PENDING*** MPI broadcasting of averaged values ...
//
//        if( cntD > 0 )
//        {
//          nd->v.D = Dave / cntD;
//        }
//      }
//    }
//
//    if( !further )  break;
//  }


  // -------------------------------------------------------------------------------------

  for( int i=0; i<np; i++ )
  {
    NODE* nd = node[i];
    int   no = nd->Getno();

    if( isFS(nd->flag, NODE::kDry) )  continue;

    if( nd->v.D < project->minD )
    {
      nd->v.K = project->minK;
      nd->v.D = project->minD;
    }
  }
}
