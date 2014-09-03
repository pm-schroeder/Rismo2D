// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_KD2D
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

#include "EqsKD2D.h"


EQS_KD2D::EQS_KD2D() : EQS( 2, 2, 2 )
{
  neq  = 0;
  cbuf = NULL;
  cent = NULL;
}


EQS_KD2D::~EQS_KD2D()
{
}


// ---------------------------------------------------------------------------------------

void EQS_KD2D::Execute( PROJECT* project, int steadyFlow, int shape )
{
  switch( shape )
  {
    case 0:
      linearShape  = false;
      quarterShape = false;
      break;
    case 1:
      linearShape  = true;
      quarterShape = false;
      break;
    case 2:
      linearShape  = false;
      quarterShape = true;
      break;
  }

  MODEL* model    = project->M2D;
  GRID*  rg       = project->M2D->region;

  NODE** node     = model->node;
  ELEM** elem     = model->elem;

  int np          = model->np;
  int ne          = model->ne;

  int diverged_cg = 0;

  double* B    = NULL;
  double* xKD  = NULL;
  double* xKDo = NULL;
  double* cxKo = NULL;
  double* cxDo = NULL;

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
                       "theta and timeInterval too small (EQS_KD2D::execute - 1)" );
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
      node[i]->v.dKdt = 0.0;
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
      VARS* v = &(node[i]->v);

      //VARS* vo = &(node[i]->vo);

      //v->K = vo->K  +  vo->dKdt * dt;
      //v->D = vo->D  +  vo->dDdt * dt;

      //v->dKdt = (1.0 - 1.0/th)*vo->dKdt + thdt*(v->K - vo->K);
      //v->dDdt = (1.0 - 1.0/th)*vo->dDdt + thdt*(v->D - vo->D);

      v->dKdt = 0.0;
      v->dDdt = 0.0;
    }
  }

  // check KD-values for validity (>= 0) -------------------------------------------------
  Validate( project, np, node );

  // determine friction coefficients -----------------------------------------------------
  model->DoFriction( project );

  // initialize Reynolds stresses and eddy viscosity -------------------------------------
  rg->Turbulence( project );

  // set KD boundary conditions ----------------------------------------------------------
  model->SetBoundKD( project );

  // -------------------------------------------------------------------------------------
  // in case of quartered elements:
  // compute averaged values of U, V, S, K and D for the virtual center node in quads

  int     nq   = 0;
  double* cxK  = NULL;
  double* cxD  = NULL;

  if( quarterShape )
  {
    // count number of quadrilaterals
    for( int e=0; e<ne; e++ )
    {
      if( elem[e]->Getncn() == 4 )  nq++;
    }

    // allocate memory for nq center nodes
    if( !cbuf )  cbuf = new NODE [nq];
    if( !cbuf )
    {
      REPORT::rpt.Error( kMemoryFault, "cannot allocate memory (EQS_KD2D::execute - 2)" );
    }

    cent = (NODE**)  MEMORY::memo.Array_el( ne );
    cxK  = (double*) MEMORY::memo.Array_el( ne );
    cxD  = (double*) MEMORY::memo.Array_el( ne );

    nq = 0;
    for( int e=0; e<ne; e++ )
    {
      ELEM* el  = elem[e];

      if( isFS(el->flag, ELEM::kRegion) )
      {
        int  no  = el->Getno();
        int  ncn = el->Getncn();

        cent[no] = NULL;

        if( ncn == 4 )
        {
          cent[no] = &cbuf[nq++];

          cent[no]->Setno( no );

          cent[no]->x      = 0.0;
          cent[no]->y      = 0.0;
          cent[no]->z      = 0.0;
          cent[no]->cf     = 0.0;
          cent[no]->v.U    = 0.0;
          cent[no]->v.V    = 0.0;
          cent[no]->v.S    = 0.0;
          cent[no]->v.K    = 0.0;
          cent[no]->v.D    = 0.0;
          cent[no]->v.dKdt = 0.0;
          cent[no]->v.dDdt = 0.0;

          for( int i=0; i<ncn; i++ )
          {
            cent[no]->x      += el->nd[i]->x;
            cent[no]->y      += el->nd[i]->y;
            cent[no]->z      += el->nd[i]->z;
            cent[no]->cf     += el->nd[i]->cf;
            cent[no]->v.U    += el->nd[i]->v.U;
            cent[no]->v.V    += el->nd[i]->v.V;
            cent[no]->v.S    += el->nd[i]->v.S;
            cent[no]->v.K    += el->nd[i]->v.K;
            cent[no]->v.D    += el->nd[i]->v.D;
            cent[no]->v.dKdt += el->nd[i]->v.dKdt;
            cent[no]->v.dDdt += el->nd[i]->v.dDdt;
          }

          cent[no]->x      /= ncn;
          cent[no]->y      /= ncn;
          cent[no]->z      /= ncn;
          cent[no]->cf     /= ncn;
          cent[no]->v.U    /= ncn;
          cent[no]->v.V    /= ncn;
          cent[no]->v.S    /= ncn;
          cent[no]->v.K    /= ncn;
          cent[no]->v.D    /= ncn;
          cent[no]->v.dKdt /= ncn;
          cent[no]->v.dDdt /= ncn;
        }
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // iteration loop

  double dt_KDo =  0.0;
  double relaxo =  1.0;
  double maxKDo = -1.0;

  int conv;

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

      project->fix[0] = BCON::kFixK;
      project->fix[1] = BCON::kFixD;

      project->elemKind = ELEM::kRegion;

      if( linearShape )
      {
        SetEqno( model, 2, 0, 0, project->fix, project->elemKind );
      }
      else
      {
        SetEqno( model, 2, 2, 0, project->fix, project->elemKind );
      }

      if( B )   MEMORY::memo.Detach( B );
      if( xKD ) MEMORY::memo.Detach( xKD );

      B   = (double*) MEMORY::memo.Array_eq( neq );
      xKD = (double*) MEMORY::memo.Array_eq( neq );

      // allocate memory for relaxed Newton-Rahpson
      if( relaxMethod >= 3 )
      {
        xKDo = (double*) MEMORY::memo.Array_eq( neq );

        if( quarterShape )
        {
          cxKo = (double*) MEMORY::memo.Array_el( ne );
          cxDo = (double*) MEMORY::memo.Array_el( ne );
        }
      }
    }

    // solve equations with frontal solving algorithm ------------------------------------
    for( int i=0; i<neq; i++ )  xKD[i] = 0.0;

    diverged_cg = Solve( model, neq, B, xKD, project );


    // -----------------------------------------------------------------------------------
    // determine new values for K and D at virtual center nodes

    if( quarterShape )
    {
      for( int e=0; e<ne; e++ )
      {
        ELEM* el = elem[e];

        if( isFS(el->flag, ELEM::kRegion) )
        {
          int  no  = el->Getno();
          int  ncn = el->Getncn();

          if( ncn == 4 )
          {
            Coefs( el, project, estifm, force );

            cxK[no] = force[16];
            cxD[no] = force[17];

            for( int i=0; i<8; i++ )
            {
              int eqK = GetEqno( el->nd[i], 0 );
              int eqD = GetEqno( el->nd[i], 1 );

              cxK[no] -= estifm[16][i] * xKD[eqK] + estifm[16][i+8] * xKD[eqD];
              cxD[no] -= estifm[17][i] * xKD[eqK] + estifm[17][i+8] * xKD[eqD];
            }

            cxK[no] /= estifm[16][16];

            cxD[no] -= estifm[17][16] * cxK[no];
            cxD[no] /= estifm[17][17];
          }
        }
      }
    }


    // -----------------------------------------------------------------------------------
    // statistics

    REPORT::rpt.Message( 2, "\n\n%-25s%s\n\n %s\n\n",
                            " (EQS_KD2D::Execute)", "convergence parameters ...",
                            " variable   node          average        maximum" );


    // compute averaged changes and standard deviation of K and D ------------------------

    double stdevK  = 0.0;
    double stdevD  = 0.0;

    double aveAbsK = 0.0;
    double aveAbsD = 0.0;

    for( int i=0; i<np; i++ )
    {
      int eqno;

      double dK = 0.0;
      double dD = 0.0;

      eqno = GetEqno( node[i], 0 );
      if( eqno >= 0 )  dK = xKD[eqno];

      eqno = GetEqno( node[i], 1 );
      if( eqno >= 0 )  dD = xKD[eqno];

      aveAbsK += dK;
      stdevK  += dK*dK;

      aveAbsD += dD;
      stdevD  += dD*dD;
    }

    if( quarterShape )
    {
      for( int e=0; e<ne; e++ )
      {
        ELEM* el = elem[e];

        if( isFS(el->flag, ELEM::kRegion) )
        {
          int  no  = el->Getno();
          int  ncn = el->Getncn();

          if( ncn == 4 )
          {
            double dK = cxK[no];
            double dD = cxD[no];

            aveAbsK += dK;
            stdevK  += dK*dK;

            aveAbsD += dD;
            stdevD  += dD*dD;
          }
        }
      }
    }


    // -----------------------------------------------------------------------------------

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

    aveAbsD = project->subdom.Mpi_sum( aveAbsD );
    stdevD  = project->subdom.Mpi_sum( stdevD );

    nptot   = project->subdom.Mpi_sum( nptot );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    aveAbsK /= nptot;
    stdevK   = sqrt( stdevK / nptot );

    aveAbsD /= nptot;
    stdevD   = sqrt( stdevD / nptot );

    double norm = stdevK + stdevD;

    double fractK = 2.0 * sqrt( stdevK );
    double fractD = 2.0 * sqrt( stdevD );


    // compute maximum changes of K and D limited to ~95% fractile -----------------------

    int    cntK    = 0;
    int    cntD    = 0;

    int    noPerK  = 0;
    double avePerK = 0.0;
    double maxPerK = 0.0;

    int    noAbsK  = 0;
    double maxAbsK = 0.0;

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
        dK = xKD[eqno];
        if( fabs(dK) > project->maxDeltaKD  &&  fabs(dK) > fractK )
        {
          xKD[eqno] = dK/fabs(dK) * fractK;
        }
      }

      eqno = GetEqno( node[i], 1 );
      if( eqno >= 0 )
      {
        dD = xKD[eqno];
        if( fabs(dD) > project->maxDeltaKD  &&  fabs(dD) > fractD )
        {
          xKD[eqno] = dD/fabs(dD) * fractD;
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

    if( quarterShape )
    {
      for( int e=0; e<ne; e++ )
      {
        ELEM* el = elem[e];

        if( isFS(el->flag, ELEM::kRegion) )
        {
          int  no  = el->Getno();
          int  ncn = el->Getncn();

          if( ncn == 4 )
          {
            double dK = cxK[no];
            double dD = cxD[no];

            if( fabs(dK) > project->maxDeltaKD  &&  fabs(dK) > fractK )
            {
              cxK[no] = dK/fabs(dK) * fractK;
            }

            if( fabs(dD) > project->maxDeltaKD  &&  fabs(dD) > fractD )
            {
              cxD[no] = dD/fabs(dD) * fractD;
            }

            // maximum changes and percentage
            if( cent[no]->v.K + dK > 0.0 )
            {
              if( fabs(dK) > fabs(maxAbsK) )
              {
                maxAbsK = dK;
                noAbsK  = -(no+1);
              }

              if( cent[no]->v.K > 0.0 )
              {
                double per = dK / cent[no]->v.K;

                avePerK += fabs(per);
                cntK++;

                if( fabs(per) > fabs(maxPerK) )
                {
                  maxPerK = per;
                  noPerK  = -(no+1);
                }
              }
            }

            if( cent[no]->v.D + dD > 0.0 )
            {
              if( fabs(dD) > fabs(maxAbsD) )
              {
                maxAbsD = dD;
                noAbsD  = -(no+1);
              }

              if( cent[no]->v.D > 0.0 )
              {
                double per = dD / cent[no]->v.D;

                avePerD += fabs(per);
                cntD++;

                if( fabs(per) > fabs(maxPerD) )
                {
                  maxPerD = per;
                  noPerD  = -(no+1);
                }
              }
            }
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: broadcast statistic
#   ifdef _MPI_
    cntK    = project->subdom.Mpi_sum( cntK );
    maxAbsK = project->subdom.Mpi_maxabs( maxAbsK );
    avePerK = project->subdom.Mpi_sum( avePerK );
    maxPerK = project->subdom.Mpi_maxabs( maxPerK );

    cntD    = project->subdom.Mpi_sum( cntD );
    maxAbsD = project->subdom.Mpi_maxabs( maxAbsD );
    avePerD = project->subdom.Mpi_sum( avePerD );
    maxPerD = project->subdom.Mpi_maxabs( maxPerD );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    avePerK /= cntK;
    avePerD /= cntD;

    REPORT::rpt.Message( 2, "      %1c     %5d       %12.5le   %12.5le %s\n",
                            'K', noAbsK+1, aveAbsK, maxAbsK, "     (abs)" );

    REPORT::rpt.Message( 2, "      %1c     %5d       %12.5lf   %12.5lf %s\n\n",
                            ' ', noPerK+1, avePerK, maxPerK, "     ( % )" );


    REPORT::rpt.Message( 2, "      %1c     %5d       %12.5le   %12.5le %s\n",
                            'D', noAbsD+1, aveAbsD, maxAbsD, "     (abs)" );

    REPORT::rpt.Message( 2, "      %1c     %5d       %12.5lf   %12.5lf %s\n\n",
                            ' ', noPerD+1, avePerD, maxPerD, "     ( % )" );


    if( fabs(maxAbsK) > project->convKD )  conv = false;
    if( fabs(maxAbsD) > project->convKD )  conv = false;


    // determine relaxation parameter for NEWTON-RAPHSON ---------------------------------

    double relax;

    double maxKD = fabs(maxAbsK);
    if( fabs(maxAbsD) > maxKD )  maxKD = fabs(maxAbsD);

    switch( relaxMethod )
    {
      default:
        REPORT::rpt.Warning( kParameterFault,
                             "relaxation method %d not supported", relaxMethod );

      case 0:
        relax = 1.0;

        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n%-25s%s %12.4le\n\n",
                                " ", "relaxation (0): norm  =", norm,
                                " ", "                relax =", relax );
        break;

      case 2:
        relax = project->maxDeltaKD / maxKD;

        if( relax > 1.0 )  relax = 1.0;

        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n%-25s%s %12.4le\n\n",
                                " ", "relaxation (2): norm  =", norm,
                                " ", "                relax =", relax );
        break;

      case 3:
        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n",
                                " ", "relaxed time:   dt_KD =", dt_KD );

        if( dt_KD < dt )  conv = false;

        relax = project->maxDeltaKD / maxKD;

        dt_KD *= relax;

        if( dt_KD > dt )          dt_KD = dt;
        if( dt_KD < relaxDt_KD )  dt_KD = relaxDt_KD;

        if( steadyFlow ) relaxThdt_KD = 1.0 / dt_KD;
        else             relaxThdt_KD = 1.0 / dt_KD / th;

        if( relax > 1.0 )  relax = 1.0;

        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n%-25s%s %12.4le\n\n",
                                " ", "relaxation (3): norm  =", norm,
                                " ", "                relax =", relax );
        break;

      case 4:
        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n",
                                " ", "relaxed time:   dt_KD =", dt_KD );

        if( dt_KD < dt )  conv = false;

        relax = project->maxDeltaKD / maxKD;

        if( relax < project->relaxMax )
        {
          if( maxKDo < 0.0  ||  maxKD < maxKDo )  // initialisation: maxKDo = -1
          {
            if( relax < project->relaxMin ) relax = project->relaxMin;

            dt_KDo = dt_KD;
            maxKDo = maxKD;
          }
          else
          {
            if( relaxo > 1.01 * project->relaxMin )
            {
              // relaxed Newton-Raphson: restore K and D from previous iteration
              for( int i=0; i<np; i++ )
              {
                int n;
                n = GetEqno( node[i], 0 );
                if( n >= 0 )
                {
                  node[i]->v.K -= relaxo * xKDo[n];
                  xKD[n] = xKDo[n];
                }
                n = GetEqno( node[i], 1 );
                if( n >= 0 )
                {
                  node[i]->v.D -= relaxo * xKDo[n];
                  xKD[n] = xKDo[n];
                }
              }

              if( quarterShape )
              {
                for( int e=0; e<ne; e++ )
                {
                  ELEM* el = elem[e];

                  if( isFS(el->flag, ELEM::kRegion) )
                  {
                    if( el->Getncn() == 4 )
                    {
                      int no = el->Getno();
                      cent[no]->v.K -= relaxo * cxKo[no];
                      cent[no]->v.D -= relaxo * cxDo[no];
                      cxK[no] = cxKo[no];
                      cxD[no] = cxDo[no];
                    }
                  }
                }
              }
            }

            if( relax < project->relaxMin ) relax = project->relaxMin;

            dt_KDo = dt_KD;
            dt_KD *= relax;           // decrease dt_KD
            if( dt_KD < relaxDt_KD ) dt_KD = relaxDt_KD;

            maxKDo = maxKD;
          }
        }
        else
        {
          dt_KDo = dt_KD;

          if( relax > 1.0 ) dt_KD *= relax;       // increase dt_KD
          if( dt_KD > dt ) dt_KD = dt;

          relax  = 1.0;

          maxKDo = maxKD;

//          if( relax > relaxo )
//          {
//            relax = 0.1 * relaxo;

//            if( relax >= project->relaxMin )
//            {
//              // relaxed Newton-Raphson: restore K and D from previous iteration

//              for( int i=0; i<np; i++ )
//              {
//                int n;
//                n = GetEqno( node[i], 0 );
//                if( n >= 0 )
//                {
//                  node[i]->v.K -= relaxo * xKDo[n];
//                  xKD[n] = xKDo[n];
//                }
//                n = GetEqno( node[i], 1 );
//                if( n >= 0 )
//                {
//                  node[i]->v.D -= relaxo * xKDo[n];
//                  xKD[n] = xKDo[n];
//                }
//              }

//              if( quarterShape )
//              {
//                for( int e=0; e<ne; e++ )
//                {
//                  ELEM* el = elem[e];

//                  if( isFS(el->flag, ELEM::kRegion) )
//                  {
//                    if( el->Getncn() == 4 )
//                    {
//                      int no = el->Getno();
//                      cent[no]->v.K -= relaxo * cxKo[no];
//                      cent[no]->v.D -= relaxo * cxDo[no];
//                      cxK[no] = cxKo[no];
//                      cxD[no] = cxDo[no];
//                    }
//                  }
//                }
//              }

//              maxKDo = maxKD;
//            }
//            else
//            {
//              maxKDo = -1.0;
//            }
//          }
//          else
//          {
//            relax  = project->relaxMax;
//            maxKDo = maxKD;
//          }

//          if( relax > project->relaxMax )  relax = project->relaxMax;
//          if( relax < project->relaxMin )  relax = project->relaxMin;

//          if( relax < 0.999 * relaxo )  conv = false;
        }

        if( steadyFlow ) relaxThdt_KD = 1.0 / dt_KD;
        else             relaxThdt_KD = 1.0 / dt_KD / th;

        REPORT::rpt.Message( 2, "\n%-25s%s %12.4le\n%-25s%s %12.4le\n\n",
                                " ", "relaxation (4): norm  =", norm,
                                " ", "                relax =", relax );

        // relaxed Newton-Raphson: store xKD, cxK and cxD
        for( int i=0; i<np; i++ )
        {
          int n;
          n = GetEqno( node[i], 0 );
          if( n >= 0 )  xKDo[n] = xKD[n];

          n = GetEqno( node[i], 1 );
          if( n >= 0 )  xKDo[n] = xKD[n];
        }

        if( quarterShape )
        {
          for( int e=0; e<ne; e++ )
          {
            ELEM* el = elem[e];

            if( isFS(el->flag, ELEM::kRegion) )
            {
              if( el->Getncn() == 4 )
              {
                int no = el->Getno();
                cxKo[no] = cxK[no];
                cxDo[no] = cxD[no];
              }
            }
          }
        }
        break;
    }


    // update ----------------------------------------------------------------------------

    relaxo = relax;

    for( int i=0; i<np; i++ )
    {
      int n;
      n = GetEqno( node[i], 0 );
      if( n >= 0 )  node[i]->v.K += relax * xKD[n];

      n = GetEqno( node[i], 1 );
      if( n >= 0 )  node[i]->v.D += relax * xKD[n];
    }


    if( quarterShape )
    {
      for( int e=0; e<ne; e++ )
      {
        ELEM* el = elem[e];

        if( isFS(el->flag, ELEM::kRegion) )
        {
          if( el->Getncn() == 4 )
          {
            int no = el->Getno();
            cent[no]->v.K += relax * cxK[no];
            cent[no]->v.D += relax * cxD[no];
          }
        }
      }
    }


    // check for range of values ---------------------------------------------------------

    double minK = 0.0;
    double minD = 0.0;
    double maxK = 0.0;
    double maxD = 0.0;

    int first = true;

    int jk = 0;
    int jd = 0;

    for( int i=0; i<np; i++ )
    {
      if( first )
      {
        first = false;
        minK =
        maxK = node[i]->v.K;

        minD =
        maxD = node[i]->v.D;
      }
      else
      {
        if( node[i]->v.K < minK )  minK = node[i]->v.K;
        if( node[i]->v.K > maxK )  maxK = node[i]->v.K;

        if( node[i]->v.D < minD )  minD = node[i]->v.D;
        if( node[i]->v.D > maxD )  maxD = node[i]->v.D;
      }

      if( node[i]->v.K <= 0.0 )
      {
        if( GetEqno(node[i],0) >= 0 )  jk++;
      }

      if( node[i]->v.D <= 0.0 )
      {
        if( GetEqno(node[i],1) >= 0 )  jd++;
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: broadcast statistic
#   ifdef _MPI_
    jk   = project->subdom.Mpi_sum( jk );
    minK = project->subdom.Mpi_min( minK );
    maxK = project->subdom.Mpi_max( maxK );

    jd   = project->subdom.Mpi_sum( jd );
    minD = project->subdom.Mpi_min( minD );
    maxD = project->subdom.Mpi_max( maxD );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    // compute useful values for K and D where they are negative -------------------------
    if( jk || jd )
    {
      REPORT::rpt.Message( 3, "%-25s%d %s\n\n", " ", jk, "nodes with K out of range" );
      REPORT::rpt.Message( 3, "%-25s%d %s\n\n", " ", jd, "nodes with D out of range" );
    }

    REPORT::rpt.Message( 3, "%-25sminimum of K: %le\n", " ", minK );
    REPORT::rpt.Message( 3, "%-25smaximum of K: %le\n", " ", maxK );
    REPORT::rpt.Message( 3, "%-25sminimum of D: %le\n", " ", minD );
    REPORT::rpt.Message( 3, "%-25smaximum of D: %le\n", " ", maxD );


    Validate( project, np, node, ne, cent, elem );


    // compute midside values (linear interpolation) -------------------------------------

    if( linearShape )
    {
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

          left = el->nd[il]->v.D;
          rght = el->nd[ir]->v.D;
          el->nd[i]->v.D = 0.5 * (left + rght);
        }
      }
    }


    // compute time derivatives ----------------------------------------------------------

    rg->ReportCuPe( dt, project->vk );


    if( !steadyFlow )
    {
      double iTheta = 1.0  -  1.0 / project->timeint.thetaTurb;

      for( int i=0; i<np; i++ )
      {
        if(     isFS(node[i]->flag, NODE::kDry)
            ||  isFS(node[i]->flag, NODE::kMarsh) )
        {
          node[i]->v.dKdt =
          node[i]->v.dDdt = 0.0;
        }
        else
        {
          double K, pK, D, pD, pdKdt, pdDdt;
          VARS *v, *vo;

          v  = &(node[i]->v);
          vo = &(node[i]->vo);

          K     = v->K;
          pK    = vo->K;
          pdKdt = vo->dKdt;

          D     = v->D;
          pD    = vo->D;
          pdDdt = vo->dDdt;

          // compute derivatives at actual time step
          node[i]->v.dKdt = iTheta*pdKdt + thdt*(K - pK);
          node[i]->v.dDdt = iTheta*pdDdt + thdt*(D - pD);
        }
      }
    }


    // -----------------------------------------------------------------------------------

    if( conv  ||  it == project->actualCycit-1 )
    {
      if( REPORT::rpt.level == 1  &&  it == project->actualCycit-1 )
      {
        REPORT::rpt.Message( 1, "\n\n%-25s%s: %d\n\n",
                                " (EQS_KD2D::Execute)", "finished in iteration step", it+1 );

        REPORT::rpt.Message( 1, "\n\n%-25s%s\n\n %s\n\n",
                                " (EQS_KD2D::Execute)", "convergence parameters ...",
                                " variable   node          average        maximum" );

        REPORT::rpt.Message( 1, "      %1c     %5d       %12.5le   %12.5le %s\n",
                                'K', noAbsK+1, aveAbsK, maxAbsK, "     (abs)" );

        REPORT::rpt.Message( 1, "      %1c     %5d       %12.5lf   %12.5lf %s\n\n",
                                ' ', noPerK+1, avePerK, maxPerK, "     ( % )" );

        REPORT::rpt.Message( 1, "      %1c     %5d       %12.5le   %12.5le %s\n",
                                'D', noAbsD+1, aveAbsD, maxAbsD, "     (abs)" );

        REPORT::rpt.Message( 1, "      %1c     %5d       %12.5lf   %12.5lf %s\n\n",
                                ' ', noPerD+1, avePerD, maxPerD, "     ( % )" );
      }

      break;
    }
  }


  // -------------------------------------------------------------------------------------
  // finally:
  // compute eddy viscosity from revised turbulence parameters

  rg->Turbulence( project );


  // -------------------------------------------------------------------------------------

  MEMORY::memo.Detach( B );
  MEMORY::memo.Detach( xKD );

  if( cxK )  MEMORY::memo.Detach( cxK );
  if( cxD )  MEMORY::memo.Detach( cxD );
  if( cent ) MEMORY::memo.Detach( cent );
  if( xKDo ) MEMORY::memo.Detach( xKDo );
  if( cxKo ) MEMORY::memo.Detach( cxKo );
  if( cxDo ) MEMORY::memo.Detach( cxDo );


  // -------------------------------------------------------------------------------------

  if( !conv )        project->errLevel |= kErr_some_errors | kErr_no_conv_nr;
  if( diverged_cg )  project->errLevel |= diverged_cg | kErr_no_conv_cg;
}


//////////////////////////////////////////////////////////////////////////////////////////
// Check KD-values for validity (>= 0).
//////////////////////////////////////////////////////////////////////////////////////////

void EQS_KD2D::Validate( PROJECT* project, int np, NODE** node,
                         int ne, NODE** cent, ELEM** elem )
{
  // -------------------------------------------------------------------------------------
  // Check K and D for minmum flow velocity (Reynolds number) in surrounding nodes

  // define the limit of Re, where K and D will be set to their minimum
  double minRe = 10.0;

  // for each node: determine the maximum Reynolds number from surrounding nodes
  double* maxRe = (double*) MEMORY::memo.Array_nd( np, "EQS_KD2D::Validate(1)" );

  for( int n=0; n<np; n++ )
  {
    NODE* nd = node[n];
    int   no = nd->Getno();

    double U = nd->v.U;
    double V = nd->v.V;
    double H = nd->v.S - nd->z;

    maxRe[no] = 4.0 * H * sqrt( U*U + V*V ) / project->vk;

    // check also surrounding nodes, if Reynolds number is too small at node ndi
    if( maxRe[no] < minRe )
    {
      for( int j=0; j<nd->noel; j++ )
      {
        ELEM* el = nd->el[j];

        for( int k=0; k<el->Getncn(); k++ )
        {
          NODE* ndk = el->Getnode(k);

          double U  = ndk->v.U;
          double V  = ndk->v.V;
          double H  = ndk->v.S - ndk->z;
          double Re = H * 4.0 * sqrt( U*U + V*V ) / project->vk;

          if( Re > maxRe[no] )
          {
            maxRe[no] = Re;
            if( maxRe[no] > minRe )  break;
          }
        }

        if( maxRe[no] > minRe )  break;
      }
    }
  }

# ifdef _MPI_
  project->subdom.Mpi_max( maxRe );
# endif

  for( int n=0; n<np; n++ )
  {
    NODE* nd = node[n];
    int   no = nd->Getno();

    // set the turbulence at nodes with very small Reynolds numbers to the minimum
    if( maxRe[no] < minRe )
    {
      nd->v.K = project->minK;
      nd->v.D = project->minD;
    }
  }

  MEMORY::memo.Detach( maxRe );


  if( quarterShape )
  {
    for( int e=0; e<ne; e++ )
    {
      ELEM* el = elem[e];
      if( isFS( el->flag, ELEM::kRegion)  &&  el->Getncn() == 4 )
      {
        int   no = el->Getno();

        NODE* nd = cent[no];

        double U = nd->v.U;
        double V = nd->v.V;
        double H = nd->v.S - nd->z;

        double c_maxRe = 4.0 * H * sqrt( U*U + V*V ) / project->vk;

        // check also surrounding nodes, if Reynolds number is too small at node ndi
        if( c_maxRe < minRe )
        {
          for( int k=0; k<el->Getncn(); k++ )
          {
            NODE* ndk = el->Getnode(k);

            double U  = ndk->v.U;
            double V  = ndk->v.V;
            double H  = ndk->v.S - ndk->z;
            double Re = H * 4.0 * sqrt( U*U + V*V ) / project->vk;

            if( Re > c_maxRe )
            {
              c_maxRe = Re;
              if( c_maxRe > minRe )  break;
            }
          }
        }

        // set the turbulence at nodes with very small Reynolds numbers to the minimum
        if( c_maxRe < minRe )
        {
          nd->v.K = project->minK;
          nd->v.D = project->minD;
        }
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // try to interpolate nodes with negative values from surrounding nodes
  // Note: This interpolation is not so helpful in all cases.

  // ***PENDING*** MPI broadcasting of averaged nodal values Kave and Dave
  //               from surrounding elements.
/*
  GRID*   rg  = project->M2D->region;

  int*    cnK = (int*)    MEMORY::memo.Array_el( rg->Getne() );
  int*    cnD = (int*)    MEMORY::memo.Array_el( rg->Getne() );
  double* elK = (double*) MEMORY::memo.Array_el( rg->Getne() );
  double* elD = (double*) MEMORY::memo.Array_el( rg->Getne() );

  int iter;
  for( iter=0; iter<0; iter++ )
  {
    int further = false;

    for( int e=0; e<rg->Getne(); e++ )
    {
      ELEM* el = rg->Getelem(e);

      int no  = el->Getno();
      int ncn = el->getncn();

      cnK[no] = cnD[no] = 0;
      elK[no] = elD[no] = 0.0;

      for( int i=0; i<ncn; i++ )
      {
        NODE* nd = el->Getnode(i);

        if( isFS(nd->flag, NODE::kDry) )  continue;

        double K = nd->v.K;

        if( K > 0.0 )
        {
          elK[no] += K;
          cnK[no]++;
        }

        double D = nd->v.D;

        if( D > 0.0 )
        {
          elD[no] += D;
          cnD[no]++;
        }
      }
    }

    for( int n=0; n<np; n++ )
    {
      NODE* nd = node[n];

      if( nd->v.K <= 0.0 )
      {
        further = true;

        int    cntK = 0;
        double Kave = 0.0;

        for( int i=0; i<nd->noel; i++ )
        {
          ELEM* el = nd->el[i];
          int   no = el->Getno();

          if( cnK[no] )
          {
            Kave += elK[no];
            cntK += cnK[no];
          }
        }

        if( cntK )  nd->v.K = Kave / cntK;
      }

      if( nd->v.D <= 0.0 )
      {
        further = true;

        int    cntD = 0;
        double Dave = 0.0;

        for( int i=0; i<nd->noel; i++ )
        {
          ELEM* el = nd->el[i];
          int   no = el->Getno();

          if( cnD[no] )
          {
            Dave += elD[no];
            cntD += cnD[no];
          }
        }

        if( cntD )  nd->v.D = Dave / cntD;
      }
    }

    if( !further )  break;
  }

  MEMORY::memo.Detach( cnK );
  MEMORY::memo.Detach( cnD );
  MEMORY::memo.Detach( elK );
  MEMORY::memo.Detach( elD );

  if( iter )
  {
    REPORT::rpt.Message( "\n%-25s%s %d %s\n", " (EQS_KD2D::Validate)",
                         "interpolation of negative values in", iter, "iterations" );
  }
*/

  // -------------------------------------------------------------------------------------
  // check K and D for minimum values

  for( int n=0; n<np; n++ )
  {
    NODE* nd = node[n];

    if( nd->v.K < project->minK  ||  nd->v.D < project->minD )
    {
      nd->v.K = project->minK;
      nd->v.D = project->minD;
    }
  }

  if( quarterShape )
  {
    for( int e=0; e<ne; e++ )
    {
      ELEM* el = elem[e];
      if( isFS( el->flag, ELEM::kRegion)  &&  el->Getncn() == 4 )
      {
        int   no = el->Getno();

        NODE* nd = cent[no];

        if( nd->v.K < project->minK  ||  nd->v.D < project->minD )
        {
          nd->v.K = project->minK;
          nd->v.D = project->minD;
        }
      }
    }
  }
}
