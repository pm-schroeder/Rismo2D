// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_KL2D
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
#include "Type.h"
#include "Vars.h"
#include "Shape.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsKL2D.h"


EQS_KL2D::EQS_KL2D() : EQS( 1, 0, 0 )
{
  neq  = 0;
}


EQS_KL2D::~EQS_KL2D()
{
}


// ---------------------------------------------------------------------------------------

void EQS_KL2D::Execute( PROJECT* project, int steadyFlow )
{
  MODEL* model    = project->M2D;
  GRID*  rg       = project->M2D->region;

  NODE** node     = model->node;
  ELEM** elem     = model->elem;

  int np          = model->np;
  int ne          = model->ne;

  int diverged_cg = 0;

  double  theNorm = -1.0;

  double cm = project->KD.cm * project->KD.cd;

  model->Incinit();

  // print information on actual iteration -----------------------------------------------
  project->PrintTheCycle( 1 );
  REPORT::rpt.PrintTime( 1 );

  // set parameters according to time integration and relaxation -------------------------
  double th = project->timeint.thetaTurb;
  double dt = project->timeint.incTime.Getsec();

  if( fabs(th) < 1.0e-10  ||  fabs(dt) < 1.0e-10 )
  {
    REPORT::rpt.Error( kParameterFault, "theta or timeInterval too small (EQS_KL2D::execute - 1)" );
  }

  double thdt = 1.0 / dt / th;

  double dt_K;
  double relaxDt_KD = project->timeint.relaxTimeTurb.Getsec();

  int relaxMethod = project->relaxMethod;

  if( steadyFlow )
  {
    if( relaxMethod >= 3 ) // stationary flow time relaxed computation
    {
      dt_K = relaxDt_KD;
      relaxThdt_KD = 1.0 / dt_K;
    }
    else
    {
      dt_K = dt;
      relaxThdt_KD = 0.0;
    }

    for( int i=0; i<np; i++ )  node[i]->v.dKdt = 0.0;
  }

  else // instationary flow
  {
    if( relaxMethod >= 3 ) // instationary flow time relaxed computation
    {
      dt_K = relaxDt_KD;
    }
    else
    {
      dt_K = dt;
    }

    relaxThdt_KD = 1.0 / dt_K / th;

    // time prediction -------------------------------------------------------------------

    for( int i=0; i<np; i++ )
    {
      VARS* v = &(node[i]->v);
/*
      VARS* vo = &(node[i]->vo);

      v->K = vo->K  +  vo->dKdt * dt;

      v->dKdt = (1.0 - 1.0/th)*vo->dKdt + thdt*(v->K - vo->K);
*/
      v->dKdt = 0.0;
    }
  }


  if( model->Getinit() != modelInit )
  {
    initStructure = true;
    modelInit = model->Getinit();
  }


  // determine friction coefficients -----------------------------------------------------

  model->DoFriction( project );


  // set KD boundary conditions ----------------------------------------------------------

  model->SetBoundKD( project );


  // set up equation numbers -------------------------------------------------------------

  project->fix[0] = BCON::kFixK;

  project->elemKind = ELEM::kRegion;

  SetEqno( model, 1, 0, 0, project->fix, project->elemKind );

  double* B = (double*) MEMORY::memo.Array_eq( neq );
  double* X = (double*) MEMORY::memo.Array_eq( neq );

  if( !B || !X )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (EQS_KL2D::execute - 2)" );


  // -------------------------------------------------------------------------------------

  for( int it=0; it<project->actualCycit; it++ )
  {
    double relax;

    int conv = true;


    // print information on actual iteration ---------------------------------------------

    if( it )
    {
      project->PrintTheCycle( it+1 );
      REPORT::rpt.PrintTime( 2 );
    }


    // compute dissipation ---------------------------------------------------------------

    Dissipation( project );


    // initialize Reynolds stresses and eddy viscosity -----------------------------------

    rg->Turbulence( project );


    // solve equations with frontal solving algorithm ------------------------------------

    for( int i=0; i<neq; i++ ) X[i] = 0.0;

    diverged_cg = Solve( model, neq, B, X, project );


    // statistics ------------------------------------------------------------------------

    int    noAbs, noPer;
    double maxAbs, maxPer, avAbs, avPer;

    REPORT::rpt.Message( 1, "\n\n %s\n\n %s\n\n",
                            "(KLCycle)       convergence parameters ...",
                            " variable   node          average        maximum" );


    Update( model, &project->subdom,
            X, 0, kVarK, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer );

    REPORT::rpt.Message( 1, "      %1c     %5d       %12.5le   %12.5le %s\n",
                            'K', noAbs, avAbs, maxAbs, "     (abs)" );

    REPORT::rpt.Message( 1, "      %1c     %5d       %12.5lf   %12.5lf %s\n\n",
                            ' ', noPer, avPer, maxPer, "     ( % )" );

    if( fabs(maxAbs) > project->convKD )  conv = false;


    // determine relaxation parameter for NEWTON-RAPHSON ---------------------------------

    double minK = 0.0;
    double maxK = 0.0;

    switch( relaxMethod )
    {
      case 2:
        theNorm = 0.0;
        maxK    = 0.0;

        for( int i=0; i<np; i++ )
        {
          int    eqno;
          double dK = 0.0;

          eqno = GetEqno( node[i], 0 );
          if( eqno >= 0 )  dK = fabs( X[eqno] );

          theNorm += dK*dK;

          if( dK > maxK )  maxK = dK;
        }

        theNorm = sqrt( theNorm );

        relax = project->maxDeltaKD / maxK;

        if( relax > 1.0 )  relax = 1.0;

        REPORT::rpt.Message( 1, "\n %s %12.4le\n %s %12.4le\n\n",
                                "                relaxation:   norm  =", theNorm,
                                "                              relax =", relax );
        break;

      case 3:
      case 4:
        REPORT::rpt.Message( 1, "\n %s %12.4le\n",
                                "(KLCycle)       relaxed time: dt_K =", dt_K );

        if( dt_K < dt )  conv = false;

        theNorm = 0.0;
        maxK    = 0.0;

        for( int i=0; i<np; i++ )
        {
          int    eqno;
          double dK = 0.0;

          eqno = GetEqno( node[i], 0 );
          if( eqno >= 0 )  dK = fabs( X[eqno] );

          theNorm += dK*dK;

          if( dK > maxK )  maxK = dK;
        }

        theNorm = sqrt( theNorm );

        relax = project->maxDeltaKD / maxK;

        dt_K *= relax;

        if( dt_K > dt )          dt_K = dt;
        if( dt_K < relaxDt_KD )  dt_K = relaxDt_KD;

        if( steadyFlow ) relaxThdt_KD = 1.0 / dt_K;
        else             relaxThdt_KD = 1.0 / dt_K / th;

        if( relax > 1.0 )  relax = 1.0;

        REPORT::rpt.Message( 1, "\n %s %12.4le\n %s %12.4le\n\n",
                                "                relaxation:   norm  =", theNorm,
                                "                              relax =", relax );
        break;

      default:
        relax = 1.0;

        theNorm = 0.0;

        for( int i=0; i<np; i++ )
        {
          int    eqno;
          double dK = 0.0;

          eqno = GetEqno( node[i], 0 );
          if( eqno >= 0 )  dK = fabs( X[eqno] );

          theNorm += dK*dK;
        }

        theNorm = sqrt( theNorm );

        REPORT::rpt.Message( 1, "\n %s %12.4le\n %s %12.4le\n\n",
                                "                relaxation:   norm  =", theNorm,
                                "                              relax =", relax );
        break;
    }


    // update ----------------------------------------------------------------------------

    for( int i=0; i<np; i++ )
    {
      int n = GetEqno( node[i], 0 );
      if( n >= 0 )  node[i]->v.K += relax * X[n];
    }


    // check for range of values ---------------------------------------------------------

    int first = true;
    int jk = 0;

    for( int i=0; i<np; i++ )
    {
      if( first )
      {
        first = false;
        minK =
        maxK = node[i]->v.K;
      }
      else
      {
        if( node[i]->v.K < minK )  minK = node[i]->v.K;
        if( node[i]->v.K > maxK )  maxK = node[i]->v.K;
      }

      if( node[i]->v.K <= 0.0 )
      {
//        node[i]->v.K = project->minK;
        node[i]->v.K = sqrt( node[i]->vt * node[i]->v.D / cm );

        if( GetEqno(node[i],0) >= 0 )  jk++;
      }
    }


    if( jk )
    {
      REPORT::rpt.Message( 1, " (KLCycle)       %d %s\n\n", jk, "nodes with K out of range");
    }


    REPORT::rpt.Message( 1, " (KLCycle)       minimum of K: %le\n", minK );
    REPORT::rpt.Message( 1, "                 maximum of K: %le\n", maxK );


    // compute midside values (linear interpolation) -------------------------------------

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = elem[e];

      int ncn = el->Getncn();
      int nnd = el->Getnnd();

      for( int i=ncn; i<nnd; i++ )
      {
        // get left and right corner node to midside node i

        int    il, ir;
        double left, rght;

        el->GetQShape()->getCornerNodes( i, &il, &ir );

        left = el->nd[il]->v.K;
        rght = el->nd[ir]->v.K;
        el->nd[i]->v.K = 0.5 * (left + rght);
      }
    }


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


  MEMORY::memo.Detach( B );
  MEMORY::memo.Detach( X );


  if( diverged_cg )  project->errLevel |= diverged_cg | kErr_no_conv_cg;
}


//////////////////////////////////////////////////////////////////////////////////////////
// compute dissipation of turbulent kinetic energy
//////////////////////////////////////////////////////////////////////////////////////////

void EQS_KL2D::Dissipation( PROJECT* project )
{
  double cm = project->KD.cm;
  double cd = project->KD.cd;

  MODEL* model = project->M2D;
  GRID*  rg    = model->region;


  // initialization ----------------------------------------------------------------------

  for( int n=0; n<rg->Getnp(); n++ )
  {
    NODE* nd = rg->Getnode(n);
    nd->v.D = 0.0;
  }

  double* ndarea = (double*) MEMORY::memo.Array_nd( rg->Getnp() );

  for( int n=0; n<rg->Getnp(); n++ )  ndarea[n] = 0.0;


  // -------------------------------------------------------------------------------------
  // loop on elements
  // -------------------------------------------------------------------------------------

  for( int e=0; e<rg->Getne(); e++ )
  {
    ELEM* el = rg->Getelem(e);

    if( !isFS(el->flag, ELEM::kDry) )
    {
      TYPE* type = TYPE::Getid( el->type );

      int nnd = el->Getnnd();

      NODE** nd = el->nd;


      // ---------------------------------------------------------------------------------
      // compute coordinates relative to first node
      // ---------------------------------------------------------------------------------

      double xmin, xmax, ymin, ymax;
      double x[kMaxNodes2D], y[kMaxNodes2D];

      xmin = xmax = x[0] = nd[0]->x;
      ymin = ymax = y[0] = nd[0]->y;

      for( int i=1; i<nnd; i++ )
      {
        if( nd[i]->x < xmin )  xmin = nd[i]->x;
        if( nd[i]->x > xmax )  xmax = nd[i]->x;

        if( nd[i]->y < ymin )  ymin = nd[i]->y;
        if( nd[i]->y > ymax )  ymax = nd[i]->y;

        x[i] = nd[i]->x - x[0];
        y[i] = nd[i]->y - y[0];
      }

      x[0] = y[0] = 0.0;


      double lm = type->lm;


      // grid depending lm ---------------------------------------------------------------

      if( isFS(project->actualTurb, BCONSET::kVtLES) )
      {
        lm *= sqrt( (xmax-xmin)*(ymax-ymin) );
      }


      double area = el->area();

      for( int i=0; i<nnd; i++ )
      {
        double K  = nd[i]->v.K;
        if( K <= 0.0 )  K = project->minK;

        double ls = lm / sqrt( sqrt(cm*cm*cm/cd) );

        nd[i]->v.D += area * cd * sqrt(K*K*K) / ls;

        ndarea[nd[i]->Getno()] += area;
      }
    }
  }


  for( int n=0; n<rg->Getnp(); n++ )
  {
    NODE* nd = rg->Getnode(n);

    if( ndarea[n] > 0.0 )  nd->v.D /= ndarea[n];
  }

  MEMORY::memo.Detach( ndarea );
}
