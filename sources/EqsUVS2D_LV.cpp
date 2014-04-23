// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_LV
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

#include "EqsUVS2D_LV.h"


EQS_UVS2D_LV::EQS_UVS2D_LV() : EQS( 2, 0, 1 )
{
  neq   = 0;
  crsm  = NULL;

  UElimEq = NULL;
  VElimEq = NULL;
  PElimEq = NULL;
}


EQS_UVS2D_LV::~EQS_UVS2D_LV()
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
  if( PElimEq )
  {
    delete[] PElimEq[0];
    delete[] PElimEq;
  }
}


// ---------------------------------------------------------------------------------------

void EQS_UVS2D_LV::Execute( PROJECT* project, int steadyFlow )
{
  MODEL*  model   = project->M2D;
  NODE**  node    = model->node;
  GRID*   rg      = model->region;

  int     np      = model->np;
  int     ne      = rg->Getne();

  int diverged_cg = 0;

  if( steadyFlow  &&  project->relaxMethod != 3 )
  {
    REPORT::rpt.Error( kUserFault, "steady flow computation (cycle 4) not supported" );
  }

  model->Incinit();

  // print information on actual iteration -----------------------------------------------
  project->PrintTheCycle ( 1 );
  REPORT::rpt.PrintTime( 1 );

  // allocate memory for partial elimination equations S ---------------------------------
  if( !modelInit )
  {
    UElimEq = new double* [ne];
    VElimEq = new double* [ne];
    PElimEq = new double* [ne];
    if( !UElimEq || !VElimEq || !PElimEq )
      REPORT::rpt.Error( "can not allocate memory - EQS_UVS2D_LV::execute(1)" );

    // maximum number of columns is 10:
    //   4 corner nodes with 2 equations on LHS
    // + 1 element equation on LHS
    // + RHS
    UElimEq[0] = new double[10*ne];
    VElimEq[0] = new double[10*ne];
    PElimEq[0] = new double[10*ne];

    if( !UElimEq[0] || !VElimEq[0] || !PElimEq[0] )
      REPORT::rpt.Error( "can not allocate memory - EQS_UVS2D_LV::execute(2)" );

    for( int i=1; i<ne; i++ )
    {
      UElimEq[i] = UElimEq[i-1] + 10;
      VElimEq[i] = VElimEq[i-1] + 10;
      PElimEq[i] = PElimEq[i-1] + 10;
    }
  }

  // -------------------------------------------------------------------------------------
  // compute element average of water elevation from nodes

  if( modelInit > 1  &&  modelInit != model->Getinit() )
  {
    for( int e=0; e<rg->Getne(); e++ )
    {
      ELEM* el = rg->Getelem(e);

      el->U = 0.0;
      el->V = 0.0;
      el->P = 0.0;

      if( !isFS(el->flag, ELEM::kDry) )
      {
        int ncn = el->Getncn();

        for( int i=0; i<ncn; i++ )
        {
          el->U += el->nd[i]->v.U;
          el->V += el->nd[i]->v.V;
          el->P += el->nd[i]->v.S;
        }

        el->U /= ncn;
        el->V /= ncn;
        el->P /= ncn;
      }
    }
  }

  if( modelInit != model->Getinit() )
  {
    initStructure = true;
    modelInit = model->Getinit();
  }


  // set variables for time integration and relaxation -----------------------------------

  double th = project->timeint.thetaFlow;
  double dt = project->timeint.incTime.Getsec();

  if( fabs(th) < 1.0e-10  ||  fabs(dt) < 1.0e-10 )
  {
    REPORT::rpt.Error( "theta or timeInterval too small - EQS_UVS2D_LV::execute(2)" );
  }


  double relaxDt_H  = project->timeint.relaxTimeFlow.Getsec();
  double relaxDt_UV = project->timeint.relaxTimeFlow.Getsec();

  double dt_H;
  double dt_UV;

  int relaxMethod = project->relaxMethod;

  if( steadyFlow )
  {
    if( relaxMethod >= 3 ) // stationary flow time relaxed computation
    {
      dt_H  = relaxDt_H;
      dt_UV = relaxDt_UV;

      relaxThdt_H  = 1.0 / dt_H;
      relaxThdt_UV = 1.0 / dt_UV;
    }

//    not supported combination
//    else
//    {
//      dt_H  = dt;
//      dt_UV = dt;
//
//      relaxThdt_H  = 0.0;
//      relaxThdt_UV = 0.0;
//    }

    for( int i=0; i<np; i++ )
    {
      node[i]->v.dUdt = 0.0;
      node[i]->v.dVdt = 0.0;
      node[i]->v.dSdt = 0.0;
    }
  }
  else// instationary flow
  {
    if( relaxMethod >= 3 ) // instationary flow time relaxed computation
    {
      dt_H  = relaxDt_H;
      dt_UV = relaxDt_UV;
    }
    else
    {
      dt_H  = dt;
      dt_UV = dt;
    }

    relaxThdt_H  = 1.0 / dt_H  / th;
    relaxThdt_UV = 1.0 / dt_UV / th;

    // time prediction -------------------------------------------------------------------
    for( int i=0; i<np; i++ )
    {
      VARS* v     = &(node[i]->v);
/*
      VARS* vo = &(node[i]->vo);
      double thdt = 1.0 / dt / th;

      v->U = vo->U  +  vo->dUdt * dt;
      v->V = vo->V  +  vo->dVdt * dt;
      v->S = vo->S  +  vo->dSdt * dt;

      v->dUdt = (1.0 - 1.0/th)*vo->dUdt + thdt*(v->U - vo->U);
      v->dVdt = (1.0 - 1.0/th)*vo->dVdt + thdt*(v->V - vo->V);
      v->dSdt = (1.0 - 1.0/th)*vo->dSdt + thdt*(v->S - vo->S);
*/
      v->dUdt =
      v->dVdt = 0.0;
      v->dSdt = 0.0;
    }
  }


  // set Dirichlet boundary conditions ---------------------------------------------------

  model->SetNormal();
  model->SetRotation();

  for( int i=0; i<np; i++ )
  {
    BCON* bc = &node[i]->bc;

    if( bc )
    {
      CF( bc->kind, BCON::kFix_1 );
      CF( bc->kind, BCON::kFix_2 );

      if( isFS(bc->kind, BCON::kFixU) )   SF( bc->kind, BCON::kFix_1 );
      if( isFS(bc->kind, BCON::kFixV) )   SF( bc->kind, BCON::kFix_2 );
    }
  }


  // set up equation numbers -------------------------------------------------------------

  project->fix[0] = BCON::kFix_1;
  project->fix[1] = BCON::kFix_2;
  project->fix[2] = 0;

  project->elemKind = ELEM::kRegion;

  SetEqno( model, 2, 0, 0, project->fix, project->elemKind );


  // allocate memory ---------------------------------------------------------------------

  double* X  = (double*) MEMORY::memo.Array_eq( neq );
  double* B  = (double*) MEMORY::memo.Array_eq( neq );

  double* dU = (double*) MEMORY::memo.Array_nd( np );
  double* dV = (double*) MEMORY::memo.Array_nd( np );

  double* dUe = (double*) MEMORY::memo.Array_el( ne );
  double* dVe = (double*) MEMORY::memo.Array_el( ne );
  double* dPe = (double*) MEMORY::memo.Array_el( ne );


  // -------------------------------------------------------------------------------------

  double Qout = -1.0;

  int cycleIter = 0;

  int conv = true;

  for( int it=0; it<project->actualCycit; it++ )
  {
    conv = true;

    // increment iteration counter (per cycle) -------------------------------------------
    cycleIter++;

    // print information on actual iteration ---------------------------------------------
    if( it )
    {
      project->PrintTheCycle ( cycleIter );
      REPORT::rpt.PrintTime( 2 );
    }

    // compute friction coefficients for each element ------------------------------------
    model->DoFriction( project );

    // initialize Reynolds stresses and eddy viscosity -----------------------------------
    if( it == 0 || isFS(project->actualTurb, BCONSET::kVtIterat) )
    {
      rg->Turbulence( project );
    }

    // set inflow and outflow condition --------------------------------------------------

    TIME bcTime = project->timeint.actualTime + project->timeint.incTime;

    project->timeint.actualBcSet->InitBcon( project, &bcTime, &Qout );

    model->SetLocation();
    model->SetNormal();
    model->SetRotation();


    // output current time intervals -----------------------------------------------------
    if( relaxMethod >= 3 )
    {
      REPORT::rpt.Message( 2, "\n %s %12.4le\n %s %12.4le\n",
                              "(EQS_UVS2D_LV)  relaxed time: dt_H  =", dt_H,
                              "                              dt_UV =", dt_UV );
    }


    // solve equations (UV) --------------------------------------------------------------

    for( int i=0; i<neq; i++ )  X[i] = 0.0;

    diverged_cg = Solve( model, neq, B, X, project );


    // rotate flow velocity (where necessary) and statistic ------------------------------

    double aveU = 0.0;
    double maxU = 0.0;
    double aveV = 0.0;
    double maxV = 0.0;
    int    noU  = 0;
    int    noV  = 0;
    int    cntU = 0;
    int    cntV = 0;

    double maxUs = 0.0;

    for( int i=0; i<np; i++ )
    {
      NODE* nd = node[i];

      int eqnoU = GetEqno( nd, 0 );
      int eqnoV = GetEqno( nd, 1 );

      double du = 0.0;
      double dv = 0.0;

      if( eqnoU >= 0 )  du = X[eqnoU];
      if( eqnoV >= 0 )  dv = X[eqnoV];

      if( isFS(nd->flag, NODE::kRotat) )
      {
        dU[i]  = nd->bc.Getrot(0,0) * du;
        dV[i]  = nd->bc.Getrot(1,0) * du;

        //dU[i] += nd->bc.Getrot(0,1) * dv;        // dv is allways zero on
        //dV[i] += nd->bc.Getrot(1,1) * dv;        // rotational boundaries
      }

      else
      {
        dU[i] = du;
        dV[i] = dv;
      }


      // maximum change in flow velocity -------------------------------------------------

      double dUs = sqrt( dU[i]*dU[i] + dV[i]*dV[i] );

      if( dUs > maxUs )  maxUs = dUs;


      // statistical values --------------------------------------------------------------

      if( eqnoU >= 0  &&  !noU )
      {
        noU  = i + 1;
        maxU = dU[i];
        aveU = fabs( dU[i] );
        cntU = 1;
      }

      else if( eqnoU >= 0 )
      {
        if( fabs(dU[i]) > fabs(maxU) )
        {
          maxU = dU[i];
          noU  = i + 1;
        }

        aveU += fabs( dU[i] );
        cntU++;
      }

      if( eqnoV >= 0  &&  !noV )
      {
        noV  = i + 1;
        maxV = dV[i];
        aveV = fabs( dV[i] );
        cntV = 1;
      }

      else if( eqnoV >= 0 )
      {
        if( fabs(dV[i]) > fabs(maxV) )
        {
          maxV = dV[i];
          noV  = i + 1;
        }

        aveV += fabs( dV[i] );
        cntV++;
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: broadcast statistic
#   ifdef _MPI_
    maxU = project->subdom.Mpi_max( maxU );
    cntU = project->subdom.Mpi_sum( cntU );
    aveU = project->subdom.Mpi_sum( aveU );

    maxV = project->subdom.Mpi_max( maxV );
    cntV = project->subdom.Mpi_sum( cntV );
    aveV = project->subdom.Mpi_sum( aveV );

    maxUs = project->subdom.Mpi_max( maxUs );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    if( cntU ) aveU /= cntU;
    if( cntV ) aveV /= cntV;


    // determine relaxation parameter ----------------------------------------------------

    double aspectUs = project->maxDeltaUV / maxUs;
    double relaxUV  = aspectUs;

    if( relaxUV > 1.0 )  relaxUV = 1.0;

    relaxUV = 1.0;   // seems to work better


    // solve for element equations and statistics ----------------------------------------

    double aveS = 0.0;
    double maxS = 0.0;
    int    noS  = 0;
    int    cntS = 0;

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = rg->Getelem(e);

      if( !isFS(el->flag, ELEM::kDry) )
      {
        double* UElimPtr = UElimEq[el->Getno()];
        double* VElimPtr = VElimEq[el->Getno()];
        double* PElimPtr = PElimEq[el->Getno()];

        // bubble triangular element

        if( el->shape == kTri )
        {
          int ncn = el->Getncn();
          int nbn = el->GetBShape()->nnd;

          double du = UElimPtr[2*nbn + 1];         // RHS
          double dv = VElimPtr[2*nbn + 1];
          double dp = PElimPtr[2*nbn + 1];

          for( int i=0; i<ncn; i++ )               // LHS
          {
            int eqnoU = GetEqno( el->nd[i], 0 );
            int eqnoV = GetEqno( el->nd[i], 1 );

            double u = 0.0;
            double v = 0.0;

            if( eqnoU >= 0 )  u = X[eqnoU];
            if( eqnoV >= 0 )  v = X[eqnoV];

            du -= UElimPtr[i]     * u;
            du -= UElimPtr[i+nbn] * v;

            dv -= VElimPtr[i]     * u;
            dv -= VElimPtr[i+nbn] * v;

            dp -= PElimPtr[i]     * u;
            dp -= PElimPtr[i+nbn] * v;
          }

          dUe[el->Getno()] = du / UElimPtr[2*nbn];

          dv -= VElimPtr[nbn-1] * dUe[el->Getno()];
          dVe[el->Getno()] = dv / VElimPtr[2*nbn];

          dp -= PElimPtr[nbn-1]   * dUe[el->Getno()];
          dp -= PElimPtr[2*nbn-1] * dVe[el->Getno()];
          dPe[el->Getno()] = dp / PElimPtr[2*nbn];
        }

        // quadrilateral element

        else
        {
          int ncn = el->Getncn();

          double dp = PElimPtr[2*ncn + 1];         // RHS

          for( int i=0; i<ncn; i++ )               // LHS
          {
            int eqnoU = GetEqno( el->nd[i], 0 );
            int eqnoV = GetEqno( el->nd[i], 1 );

            double u = 0.0;
            double v = 0.0;

            if( eqnoU >= 0 )  u = X[eqnoU];   //relaxUV * X[eqnoU];
            if( eqnoV >= 0 )  v = X[eqnoV];   //relaxUV * X[eqnoV];

            dp -= PElimPtr[i] * u;
            dp -= PElimPtr[i+ncn] * v;
          }

          dPe[el->Getno()] = dp / PElimPtr[2*ncn];
        }


        // statistical values ------------------------------------------------------------

        if( !noS )
        {
          noS  = el->Getno() + 1;
          maxS = dPe[el->Getno()];
          aveS = fabs( dPe[el->Getno()] );
          cntS = 1;
        }

        else
        {
          if( fabs(dPe[el->Getno()]) > fabs(maxS) )
          {
            maxS = dPe[el->Getno()];
            noS  = el->Getno() + 1;
          }

          aveS += fabs( dPe[el->Getno()] );
          cntS++;
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: broadcast statistic
#   ifdef _MPI_
    maxS = project->subdom.Mpi_max( maxS );
    cntS = project->subdom.Mpi_sum( cntS );
    aveS = project->subdom.Mpi_sum( aveS );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    if( cntS ) aveS /= cntS;


    // determine relaxation parameter for flow depth -------------------------------------

    double relaxH = 1.0;

    switch( relaxMethod )
    {
      default:
        REPORT::rpt.Warning( kParameterFault, "relaxation method %d not supported", relaxMethod );
        break;

      case 2:
        relaxH  = project->maxDeltaS / fabs(maxS);
//      relaxH  = project->maxDeltaS / fabs(aveS);

        if( relaxH > 1.0 )  relaxH = 1.0;
        break;

      case 3:
      case 4:
        {
          double aspectH = project->maxDeltaS / fabs(maxS);
//        double aspectH = project->maxDeltaS / fabs(aveS);

          relaxH = aspectH;
          if( relaxH > 1.0 ) relaxH = 1.0;

          if( dt_H < dt )  conv = false;

               if( aspectH < 0.5 )  aspectH = 0.5;
          else if( aspectH < 1.0 )  aspectH = 1.0;
          else if( aspectH < 1.5 )  aspectH = 1.5;

          dt_H *= aspectH;

          if( dt_H > dt )         dt_H = dt;
          if( dt_H < relaxDt_H )  dt_H = relaxDt_H;

          // -----------------------------------------------------------------------------

          if( dt_UV < dt )  conv = false;

               if( aspectUs < 0.5 )  aspectUs = 0.5;
          else if( aspectUs < 1.0 )  aspectUs = 1.0;
          else if( aspectUs < 1.5 )  aspectUs = 1.5;

          dt_UV *= aspectUs;

          if( dt_UV > dt )          dt_UV = dt;
          if( dt_UV < relaxDt_UV )  dt_UV = relaxDt_UV;

          dt_UV = dt_H;

          // -----------------------------------------------------------------------------

          if( steadyFlow )
          {
            relaxThdt_H  = 1.0 / dt_H;
            relaxThdt_UV = 1.0 / dt_UV;
          }
          else
          {
            relaxThdt_H  = 1.0 / dt_H  / th;
            relaxThdt_UV = 1.0 / dt_UV / th;  break;
          }
        }
        break;
    }


    // apply correction to flow velocity -------------------------------------------------

    int exUV = 0;

    for( int i=0; i<np; i++ )
    {
      NODE* nd = node[i];

      nd->v.U += relaxUV * dU[i];
      nd->v.V += relaxUV * dV[i];

      double U  = nd->v.U;
      double V  = nd->v.V;
      double Us = sqrt( U*U + V*V );

      if( Us > project->maxUs )
      {
        exUV++;

        nd->v.U = U * project->maxUs / Us;
        nd->v.V = V * project->maxUs / Us;
      }
    }


    // apply correction to element center values -----------------------------------------

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = rg->Getelem(e);

      if( !isFS(el->flag, ELEM::kDry) )
      {
        el->U += relaxUV * dUe[el->Getno()];
        el->V += relaxUV * dVe[el->Getno()];
        el->P += relaxH  * dPe[el->Getno()];
      }
    }


    // compute change of water elevation at nodes ----------------------------------------

    int exS = PToNode( project );


    // output statistics -----------------------------------------------------------------

    if( fabs(maxU) > project->convUV )  conv = false;
    if( fabs(maxV) > project->convUV )  conv = false;
    if( fabs(maxS) > project->convS )   conv = false;

    REPORT::rpt.Message( 2, "\n\n %s\n\n %s\n\n",
                            "(EQS_UVS2D_LV)  convergence parameters ...",
                            " variable   node           average           maximum" );

    REPORT::rpt.Message( 2, "     %2s    %6d    %16.4le  %16.4le %s\n",
                            " U", noU, aveU, maxU, "     (abs)" );

    REPORT::rpt.Message( 2, "     %2s    %6d    %16.4le  %16.4le %s\n",
                            " V", noV, aveV, maxV, "     (abs)" );


    REPORT::rpt.Message( 2, "\n\n %s\n\n",
                            " variable   elem           average           maximum" );

    REPORT::rpt.Message( 2, "     %2s    %6d    %16.4le  %16.4le %s\n",
                            " S", noS, aveS, maxS, "     (abs)" );


    if( exUV )
    {
      REPORT::rpt.Message( 2, "\n %s %d-times exceeded\n",
                              "                limit for maximum velocity", exUV );
    }

    if( exS )
    {
      REPORT::rpt.Message( 2, "\n %s %d nodes\n",
                              "                negative flow depth at", exS );
    }

//  REPORT::rpt.Message( 2, "\n %s %12.4le\n %s %12.4le\n",
//                          "                relaxation:  relaxH  =", relaxH,
//                          "                             relaxUV =", relaxUV );
    REPORT::rpt.Message( 2, "\n %s %12.4le\n",
                            "                relaxation:  relaxH  =", relaxH );


    // compute values at midside nodes ---------------------------------------------------

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = rg->Getelem(e);

      int ncn = el->Getncn();
      int nnd = el->Getnnd();

      for( int i=ncn; i<nnd; i++ )
      {
        int    lnd, rnd;
        double left, rght;


        // --- get left and right corner node to midside node i --------------------------

        el->GetQShape()->getCornerNodes( i, &lnd, &rnd );

        left = el->nd[lnd]->v.S;
        rght = el->nd[rnd]->v.S;
        el->nd[i]->v.S = 0.5 * (left + rght);

        left = el->nd[lnd]->z;
        rght = el->nd[rnd]->z;
        el->nd[i]->z = 0.5 * (left + rght);

        left = el->nd[lnd]->v.U;
        rght = el->nd[rnd]->v.U;
        el->nd[i]->v.U = 0.5 * (left + rght);

        left = el->nd[lnd]->v.V;
        rght = el->nd[rnd]->v.V;
        el->nd[i]->v.V = 0.5 * (left + rght);
      }
    }


    // compute time derivatives ----------------------------------------------------------

    project->M2D->region->ReportCuPe( dt, project->vk );


    if( !steadyFlow )
    {
      double iTheta = 1.0  -  1.0 / project->timeint.thetaFlow;

      for( int i=0; i<np; i++ )
      {
        if(     isFS(node[i]->flag, NODE::kDry)
            ||  isFS(node[i]->flag, NODE::kMarsh) )
        {
          node[i]->v.dUdt =
          node[i]->v.dVdt =
          node[i]->v.dSdt = 0.0;
        }

        else
        {
          VARS* v  = &(node[i]->v);
          VARS* vo = &(node[i]->vo);

          double U     = v->U;
          double pU    = vo->U;
          double pdUdt = vo->dUdt;

          double V     = v->V;
          double pV    = vo->V;
          double pdVdt = vo->dVdt;

          double S     = v->S;
          double pS    = vo->S;
          double pdSdt = vo->dSdt;

          // compute derivatives at actual time step
          node[i]->v.dUdt = iTheta*pdUdt + relaxThdt_UV*(U - pU);
          node[i]->v.dVdt = iTheta*pdVdt + relaxThdt_UV*(V - pV);
          node[i]->v.dSdt = iTheta*pdSdt + relaxThdt_H *(S - pS);
        }
      }
    }


    // determine discharge through continuity lines --------------------------------------

    model->Continuity();


    // check convergence -----------------------------------------------------------------

    if( conv )  break;
  }


  MEMORY::memo.Detach( B );
  MEMORY::memo.Detach( X );

  MEMORY::memo.Detach( dU );
  MEMORY::memo.Detach( dV );

  MEMORY::memo.Detach( dUe );
  MEMORY::memo.Detach( dVe );
  MEMORY::memo.Detach( dPe );


  // report time gradients ---------------------------------------------------------------

  if( !steadyFlow )
  {
    int    noU = 0;
    int    noV = 0;
    int    noS = 0;

    double maxU = 0.0;
    double aveU = 0.0;
    double maxV = 0.0;
    double aveV = 0.0;
    double maxS = 0.0;
    double aveS = 0.0;

    int cnt = 0;

    for( int i=0; i<np; i++ )
    {
      NODE* nd = node[i];

      if( !isFS(nd->flag, NODE::kDry)  &&  isFS(nd->flag, NODE::kCornNode) )
      {
        if( fabs(nd->v.dUdt) > fabs(maxU) ) { maxU = nd->v.dUdt; noU = i; }
        if( fabs(nd->v.dVdt) > fabs(maxV) ) { maxV = nd->v.dVdt; noV = i; }
        if( fabs(nd->v.dSdt) > fabs(maxS) ) { maxS = nd->v.dSdt; noS = i; }

        aveU += fabs(nd->v.dUdt);
        aveV += fabs(nd->v.dVdt);
        aveS += fabs(nd->v.dSdt);

        cnt++;
      }
    }

    if( cnt )
    {
      aveU /= cnt;
      aveV /= cnt;
      aveS /= cnt;

      REPORT::rpt.Message( 1, "\n\n %s\n\n %s\n\n",
                              "(EQS_UVS2D_LV)  time gradients ...",
                              " variable      average        maximum  at  node" );

      REPORT::rpt.Message( 1, "      %1c    %14.4le  %14.4le  %6d\n",
                              'U', aveU, maxU, noU );

      REPORT::rpt.Message( 1, "      %1c    %14.4le  %14.4le  %6d\n",
                              'V', aveV, maxV, noV );

      REPORT::rpt.Message( 1, "      %1c    %14.4le  %14.4le  %6d\n",
                              'S', aveS, maxS, noS );
    }
  }


  rg->SetSlipFlow();


  if( !conv )        project->errLevel |= kErr_some_errors | kErr_no_conv_nr;
  if( diverged_cg )  project->errLevel |= diverged_cg | kErr_no_conv_cg;
}



// =======================================================================================
// Solve for water elevation at nodes.
// =======================================================================================

int EQS_UVS2D_LV::PToNode( PROJECT* project )
{
  MODEL*  model  = project->M2D;
  GRID*   rg     = model->region;
  int     np     = rg->Getnp();
  int     ne     = rg->Getne();
  DRYREW* dryRew = &rg->dryRew;

  double* HNode = (double*) MEMORY::memo.Array_nd( np );
  double* ANode = (double*) MEMORY::memo.Array_nd( np );

  for( int i=0; i<np; i++ )
  {
    HNode[i] =
    ANode[i] = 0.0;
  }


  // -------------------------------------------------------------------------------------
  // loop on elements
  // -------------------------------------------------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = rg->Getelem(e);

    if( !isFS(el->flag, ELEM::kDry) )
    {
      SHAPE* lShape = el->GetLShape();
      SHAPE* qShape = el->GetQShape();

      int ngp = lShape->ngp;              // number of GAUSS points
      int ncn = lShape->nnd;              // number of corner nodes
      int nnd = qShape->nnd;              // number of nodes

      double Se = el->P;                  // element water elevation

      NODE** nd = el->nd;


      // ---------------------------------------------------------------------------------
      // compute coordinates relative to first node

      double x[kMaxNodes2D], y[kMaxNodes2D];

      x[0] = nd[0]->x;
      y[0] = nd[0]->y;

      for( int i=1; i<ncn; i++ )
      {
        x[i] = nd[i]->x - x[0];
        y[i] = nd[i]->y - y[0];
      }

      x[0] = y[0] = 0.0;


      // ---------------------------------------------------------------------------------
      // GAUSS point integration

      for( int i=0; i<ngp; i++ )
      {
        // -------------------------------------------------------------------------------
        // form JACOBIAN transformation matrix with linear shape functions

        double* dfdxPtr = lShape->dfdx[i];
        double* dfdyPtr = lShape->dfdy[i];

        double trafo[2][2];

        double detj = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );

        double weight = detj * lShape->weight[i];

        for( int j=0; j<nnd; j++ )
        {
          HNode[nd[j]->Getno()] += weight * ( Se - nd[j]->z );
          ANode[nd[j]->Getno()] += weight;
        }
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // compute water elevation at nodes
  // -------------------------------------------------------------------------------------

  int exS  = 0;

  for( int i=0; i<np; i++ )
  {
    NODE* nd = rg->Getnode(i);

    if( !isFS(nd->flag, NODE::kDry) )
    {
      double S = HNode[i] / ANode[i]  +  nd->z;

      double H = S - nd->zor;

      if( dryRew->method == 1 )
      {
        if( H <= 0.0 )
        {
          exS++;
          S = nd->z + project->hmin;
        }
      }

      else if( dryRew->method == 2 )
      {
        if( H < dryRew->dryLimit )
        {
          exS++;
          nd->z = S - dryRew->dryLimit;
        }

        else
        {
          nd->z = nd->zor;
        }
      }

      nd->v.S = S;
    }
  }


  MEMORY::memo.Detach( HNode );
  MEMORY::memo.Detach( ANode );


  return exS;
}
