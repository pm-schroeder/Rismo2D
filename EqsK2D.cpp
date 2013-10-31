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

#include "EqsK2D.h"


EQS_K2D::EQS_K2D() : EQS( 1, 1, 0 )
{
  neq  = 0;
}


EQS_K2D::~EQS_K2D()
{
}


// ---------------------------------------------------------------------------------------

void EQS_K2D::Execute( PROJECT* project, int steadyFlow, int linearShape )
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
    REPORT::rpt.Error( kParameterFault,
                       "theta and timeInterval too small (EQS_K2D::execute - 1)" );
  }

  double thdt = 1.0 / dt / th;

  double dt_KD;
  double relaxDt_KD = project->timeint.relaxTimeTurb.Getsec();

  int relaxMethod = project->relaxMethod;

  if( steadyFlow )
  {
    if( relaxMethod >= 3 ) // stationary flow time relaxed computation
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
      node[i]->v.dKdt = 0.0;
    }
  }
  else // instationary flow
  {
    if( relaxMethod >= 3 ) // instationary flow time relaxed computation
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
      //v->K = vo->K  +  vo->dKdt * dt;
      //v->dKdt = (1.0 - 1.0/th)*vo->dKdt + thdt*(v->K - vo->K);

      v->dKdt = 0.0;
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
      REPORT::rpt.PrintTime( 2 );
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

      project->fix[0] = BCON::kFixK;

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

    REPORT::rpt.Message( 2, "\n\n%-25s%s\n\n %s\n\n",
                            " (EQS_K2D::Execute)","convergence parameters ...",
                            " variable   node           average          maximum" );


    // compute averaged changes and standard deviation of K ------------------------------
    double stdevK  = 0.0;
    double aveAbsK = 0.0;

    for( int i=0; i<np; i++ )
    {
      int eqno;
      double dK = 0.0;

      eqno = GetEqno( node[i], 0 );
      if( eqno >= 0 )  dK = X[eqno];

      aveAbsK += dK;
      stdevK  += dK*dK;
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
    aveAbsK = project->subdom.Mpi_sum( aveAbsK );
    stdevK  = project->subdom.Mpi_sum( stdevK );

    nptot   = project->subdom.Mpi_sum( nptot );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    aveAbsK /= nptot;
    stdevK   = sqrt( stdevK / nptot );

    double norm = stdevK;

    double fractK = 2.0 * sqrt( stdevK );


    // compute maximum changes of K and D limited to ~95% fractile -----------------------

    int    cntK    = 0;

    int    noPerK  = 0;
    double avePerK = 0.0;
    double maxPerK = 0.0;

    int    noAbsK  = 0;
    double maxAbsK = 0.0;

    for( int i=0; i<np; i++ )
    {
      int eqno;
      double dK = 0.0;

      eqno = GetEqno( node[i], 0 );
      if( eqno >= 0 )
      {
        dK = X[eqno];
        if( fabs(dK) > project->maxDeltaKD  &&  fabs(dK) > fractK )
        {
          X[eqno] = dK/fabs(dK) * fractK;
        }
      }

      // maximum changes and percentage
      if( node[i]->v.K + dK > 0.0 )
      {
        if( fabs(dK) > fabs(maxAbsK) )
        {
          maxAbsK = dK;
          noAbsK = node[i]->Getname();
        }

        if( node[i]->v.K > 0.0 )
        {
          double per = dK / node[i]->v.K;

          avePerK += fabs(per);
          cntK++;

          if( fabs(per) > fabs(maxPerK) )
          {
            maxPerK = per;
            noPerK = node[i]->Getname();
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: broadcast statistic
#   ifdef _MPI_
    cntK    = project->subdom.Mpi_sum( cntK );
    maxAbsK = project->subdom.Mpi_max( maxAbsK );
    avePerK = project->subdom.Mpi_sum( avePerK );
    maxPerK = project->subdom.Mpi_max( maxPerK );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    avePerK /= cntK;

    REPORT::rpt.Message( 2, "      %1c     %5d       %12.5le   %12.5le %s\n",
                            'K', noAbsK+1, aveAbsK, maxAbsK, "     (abs)" );

    REPORT::rpt.Message( 2, "      %1c     %5d       %12.5lf   %12.5lf %s\n\n",
                            ' ', noPerK+1, avePerK, maxPerK, "     ( % )" );


    if( fabs(maxAbsK) > project->convKD )  conv = false;


    // determine relaxation parameter for NEWTON-RAPHSON ---------------------------------

    double relax;
    double maxKD = fabs(maxAbsK);

    switch( relaxMethod )
    {
      default:
        REPORT::rpt.Warning( kParameterFault, "relaxation method %d not supported", relaxMethod );

      case 0:
        relax = 1.0;

        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n %s %12.4le\n\n",
                                " ", "relaxation:   norm  =", norm,
                                " ", "               relax =", relax );
        break;

      case 2:
        relax = project->maxDeltaKD / maxKD;

        if( relax > 1.0 )  relax = 1.0;

        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n%-25s%s %12.4le\n\n",
                                " ", "relaxation:   norm  =", norm,
                                " ", "               relax =", relax );
        break;

      case 3:
      case 4:
        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n",
                                " ", "relaxed time: dt_KD =", dt_KD );

        if( dt_KD < dt )  conv = false;

        relax = project->maxDeltaKD / maxKD;

        dt_KD *= relax;

        if( dt_KD > dt )          dt_KD = dt;
        if( dt_KD < relaxDt_KD )  dt_KD = relaxDt_KD;

        if( steadyFlow ) relaxThdt_KD = 1.0 / dt_KD;
        else             relaxThdt_KD = 1.0 / dt_KD / th;

        if( relax > 1.0 )  relax = 1.0;

        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n%-25s%s %12.4le\n\n",
                                " ", "relaxation:   norm  =", norm,
                                " ", "              relax =", relax );
        break;
    }


    // update ----------------------------------------------------------------------------

    for( int i=0; i<np; i++ )
    {
      int n = GetEqno( node[i], 0 );
      if( n >= 0 )  node[i]->v.K += relax * X[n];

    }


    // check for range of values ---------------------------------------------------------

    double minK = 0.0;
    double maxK = 0.0;

    int first = true;

    int jk = 0;

    for( int i=0; i<np; i++ )
    {
      if( first )
      {
        first = false;
        minK  =
        maxK  = node[i]->v.K;
      }
      else
      {
        if( node[i]->v.K < minK )  minK = node[i]->v.K;
        if( node[i]->v.K > maxK )  maxK = node[i]->v.K;

      }

      if( node[i]->v.K <= 0.0 )
      {
//        node[i]->v.K = project->minK;
        if( GetEqno(node[i],0) >= 0 )  jk++;
      }
    }



    // compute useful values for K and D where they are negative -------------------------
    if( jk )
    {
      REPORT::rpt.Message( 2, "%-25s%d %s\n\n", " ", jk, "nodes with K out of range" );
      Validate( np, node, project );
    }

    REPORT::rpt.Message( 2, "%-25sminimum of K: %le\n", " ", minK );
    REPORT::rpt.Message( 2, "%-25smaximum of K: %le\n", " ", maxK );


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
          node[i]->v.dKdt = 0.0;
        }

        else
        {
          double K, pK, pdKdt;
          VARS *v, *vo;

          v  = &(node[i]->v);
          vo = &(node[i]->vo);

          K     = v->K;
          pK    = vo->K;
          pdKdt = vo->dKdt;

          // compute derivatives at actual time step
          node[i]->v.dKdt = iTheta*pdKdt + thdt*(K - pK);
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
// check K-values for validity (>= 0)
//////////////////////////////////////////////////////////////////////////////////////////

void EQS_K2D::Validate( int np, NODE** node, PROJECT* project )
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
//      if( nd->v.K <= 0.0 )
//      {
//        further = true;
//
//        double Kave = 0.0;
//        int    cntK = 0;
//        int    noel = nd->noel;
//
//        for( int j=0; j<noel; j++ )
//        {
//          ELEM* el = nd->el[j];
//
//          for( int k=0; k<el->getnnd(); k++ )
//          {
//            double K = el->nd[k]->v.K;
//
//            if( K > 0.0 )
//            {
//              Kave += K;
//              cntK++;
//            }
//          }
//        }
//
//        // ***PENDING*** MPI broadcasting of averaged values ...
//
//        if( cntK > 0 )
//        {
//          nd->v.K = Kave / cntK;
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

    if( nd->v.K < project->minK )
    {
      nd->v.K = project->minK;
    }
  }
}
