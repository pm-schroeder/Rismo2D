// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_DZ
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
#include "Type.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsDz.h"


EQS_DZ::EQS_DZ() : EQS( 1, 1, 0 )
{
  neq      = 0;
  initDzQb = true;
  etaQb   = NULL;
}


EQS_DZ::~EQS_DZ()
{
}

// ---------------------------------------------------------------------------------------

void EQS_DZ::Execute( PROJECT* project, int shape )
{
  MODEL*   model  = project->M2D;

  GRID*    rg     = model->region;
  int      rgnp   = rg->Getnp();
  int      rgne   = rg->Getne();

  SED*     sed    = &project->sed;

  SUBDOM*  subdom = &project->subdom;

  int      diverged_cg = 0;
  char     text[500];

  double   dt;

  model->Incinit();

  // print information on actual iteration -----------------------------------------------
  project->PrintTheCycle( 1 );
  REPORT::rpt.PrintTime( 1 );

  // bed evolution -----------------------------------------------------------------------
  dzdt  = (double*) MEMORY::memo.Array_nd( rgnp, "EQS_DZ::Execute(1)" );
  for( int n=0; n<rgnp; n++ )  dzdt[n] = 0.0;

  // reduction factor to account for transport on rigid bed ------------------------------
  etaQb = (double*) MEMORY::memo.Array_nd( rgnp, "EQS_DZ::Execute(2)" );

  for( int n=0; n<rgnp; n++ )  etaQb[n] = 1.0;

  ////////////////////////////////////////////////////////////////////////////////////////
  // solve bottom evolution equation iteratively ...

  switch( sed->exnerEq )
  {
    //////////////////////////////////////////////////////////////////////////////////////
    // no change of bed
    //////////////////////////////////////////////////////////////////////////////////////
    case 0:
      break;

    //////////////////////////////////////////////////////////////////////////////////////
    // Finite-Volumes: element wise dzdt
    //////////////////////////////////////////////////////////////////////////////////////
    case 1:
    {
      double* A    = (double*) MEMORY::memo.Array_el( rgne, "EQS_DZ::Execute(3)" );  // element area
      double* Vtot = (double*) MEMORY::memo.Array_el( rgne, "EQS_DZ::Execute(4)" );  // change of volume (total)
      double* Vmax = (double*) MEMORY::memo.Array_el( rgne, "EQS_DZ::Execute(5)" );  // max change of volume
      double* dzdt = (double*) MEMORY::memo.Array_el( rgne, "EQS_DZ::Execute(6)" );  // change of elevation
      double* dz   = (double*) MEMORY::memo.Array_nd( rgnp, "EQS_DZ::Execute(7)" );  // change of elevation at node
      double* wght = (double*) MEMORY::memo.Array_nd( rgnp, "EQS_DZ::Execute(8)" );  // weighting factor per node

      double** V   = (double**) MEMORY::memo.Dmatrix( rgne, kMaxNodes2D );  // node wise change

      for( int it=0; it<project->actualCycit; it++ )
      {
        // compute change of sediment volume per element ---------------------------------
        if( it == 0 )
        {
          for( int e=0; e<rgne; e++ )
          {
            ELEM* elem = rg->Getelem(e);
            Vtot[e] = sed->Balance( elem, shape, etaQb, V[e], &A[e], &Vmax[e] );
          }
        }
        else
        {
          for( int e=0; e<rgne; e++ )
          {
            ELEM* elem = rg->Getelem(e);
            Vtot[e] = sed->Balance( elem, shape, etaQb, V[e] );
          }
        }

        /* erosion/deposition due to mud exchange with water body ------------------------
        if( sed->M > 0.0 )
        {
          for( int e=0; e<rgne; e++ )
          {
            ELEM* elem = rg->Getelem(e);

            int ncn = elem->getncn();
            int nnd = elem->getnnd();

            SHAPE* lShape = elem->getLShape();
            SHAPE* qShape = elem->getQShape();

            // local element coordinates -------------------------------------------------
            double x[kMaxNodes2D], y[kMaxNodes2D];

            x[0] = elem->nd[0]->x;
            y[0] = elem->nd[0]->y;

            for( int i=1; i<nnd; i++ )
            {
              x[i] = elem->nd[i]->x - x[0];
              y[i] = elem->nd[i]->y - y[0];
            }
            x[0] = y[0] = 0.0;


            // integrate over element area: A[e] -----------------------------------------
            for( int g=0; g<lShape->ngp; g++ )
            {
              double* dfdxPtr = lShape->dfdx[g];
              double* dfdyPtr = lShape->dfdy[g];

              double  trafo[2][2];

              double  detj   = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
              double  weight = detj * lShape->weight[g];

              double* m = lShape->f[g];
              double* n = qShape->f[g];

              double Hg = 0.0;
              double Mg = 0.0;

              for( int i=0; i<ncn; i++ )
              {
                NODE* nd = elem->nd[i];

                Hg += m[i] * (nd->v.S - nd->z);

                if( nd->zor > nd->zero )  Mg += m[i] * sed->M;
              }

              if( Hg < 0.0 )  Hg = 0.0;

              double Ug = 0.0;
              double Vg = 0.0;
              double Cg = 0.0;

              for( int i=0; i<ncn; i++ )
              {
                NODE* nd = elem->nd[i];

                Ug += n[i] * nd->v.U;
                Vg += n[i] * nd->v.V;
                Cg += n[i] * nd->v.C;
              }

              double Us   = sqrt( Ug*Ug + Vg*Vg );

              double Utau = sed->GetUtau( Us, Hg, 0.0, project );
              double tau  = project->rho * Utau * Utau;

              double Es = 0.0;     // erosion
              double Ds = 0.0;     // deposition

              if( sed->tauc > 0.0  &&  tau > sed->tauc  &&  Mg > 0.0 )
                Es = Mg * (tau/sed->tauc - 1.0);

              if( sed->taus > 0.0  &&  tau < sed->taus  &&  Cg > 0.0 )
                Ds = sed->us * (1.0 - tau/sed->taus) * Cg;

              Vtot[e] += weight * (Es - Ds) / sed->rhob;
            }
          }
        }*/

        // compute morphological time step -----------------------------------------------
        for( int e=0; e<rgne; e++ )  dzdt[e] = Vtot[e] / A[e] / (1.0 - sed->por);

        dt = MorphTime( project, rgne, dzdt );

        // set reduction factor to account for transport on rigid bed --------------------
        int    cnt    = 0;
        double maxAbs = 0.0;
        double aveAbs = 0.0;

        if( shape == kLinear )
        {
          for( int e=0; e<rgne; e++ )
          {
            ELEM*  elem = rg->Getelem(e);
            NODE** nd   = elem->nd;
            int    ncn  = elem->Getncn();

            if( dt*Vtot[e] < -Vmax[e] )     // not enough sediment volume available ?
            {                               // note: Vtot[e] < 0 in case of erosion
              double a = Vmax[e];
              double b = 0.0;

              for( int i=0; i<ncn; i++ )
              {
                if( V[e][i] > 0.0 )  a += dt*V[e][i];
                else                 b -= dt*V[e][i];
              }

              for( int i=0; i<ncn; i++ )
              {
                int no = nd[i]->Getno();

                if( a/b < etaQb[no] )
                {
                  double delta = etaQb[no] - a/b;
                  cnt++;
                  aveAbs += delta;
                  if( delta > maxAbs )  maxAbs = delta;

                  etaQb[no] = a/b;
                }
              }
            }
          }
        }
        else
        {
          for( int e=0; e<rgne; e++ )
          {
            ELEM*  elem = rg->Getelem(e);
            NODE** nd   = elem->nd;
            int    nnd  = elem->Getnnd();

            if( dt*Vtot[e] < -Vmax[e] )     // not enough sediment volume available ?
            {                               // note: Vtot[e] < 0 in case of erosion
              double a = Vmax[e];
              double b = 0.0;

              for( int i=0; i<nnd; i++ )
              {
                if( V[e][i] > 0.0 )  a += dt*V[e][i];
                else                 b -= dt*V[e][i];
              }

              for( int i=0; i<nnd; i++ )
              {
                int no = nd[i]->Getno();

                if( a/b < etaQb[no] )
                {
                  double delta = etaQb[no] - a/b;
                  cnt++;
                  aveAbs += delta;
                  if( delta > maxAbs )  maxAbs = delta;

                  etaQb[no] = a/b;
                }
              }
            }
          }
        }
        //////////////////////////////////////////////////////////////////////////////////
#       ifdef _MPI_
        maxAbs = subdom->Mpi_max( maxAbs );
        aveAbs = subdom->Mpi_sum( aveAbs );
        cnt    = subdom->Mpi_sum( cnt );
#       endif
        //////////////////////////////////////////////////////////////////////////////////

        if( cnt )  aveAbs /= cnt;

        sprintf( text, "\n\n%-25s%s\n\n %s\n\n",
                        " (EQS_DZ::Execute)", "convergence parameters ...",
                        "  variable         average        maximum" );
        REPORT::rpt.Message( 1, text );

        sprintf( text, "    %5s        %12.5le   %12.5le %s\n",
                      "etaQb", aveAbs, maxAbs, "     (abs)" );
        REPORT::rpt.Message( 1, text );

        // check for convergence ---------------------------------------------------------
        if( maxAbs < 1.0e-6 )  break;
      }

      // compute the changed bed elevation in element average ----------------------------
      for( int e=0; e<rgne; e++ )
      {
        ELEM* elem = rg->Getelem(e);
        elem->dz = dt * dzdt[e];
      }

      // compute changed bed elevation at nodes ------------------------------------------
      for( int n=0; n<rgnp; n++ )
      {
        dz[n] = wght[n] = 0.0;
      }

      for( int e=0; e<rgne; e++ )
      {
        ELEM* elem = rg->Getelem(e);

        int ncn = elem->Getncn();

        for( int i=0; i<ncn; i++ )
        {
          NODE* nd = elem->Getnode(i);
          int   no = nd->Getno();

          // don't change Z in case of erosion on rigid bed ------------------------------
          if( elem->dz < 0.0  &&  sed->Erode(nd) <= 0.0 )  continue;

          dz[no]   += A[e] * elem->dz;
          wght[no] += A[e];
        }
      }

      for( int n=0; n<rgnp; n++ )
      {
        if( wght[n] > 0.0 )  dz[n] /= wght[n];
      }

      // detach temporarily used memory --------------------------------------------------
      MEMORY::memo.Detach( A );
      MEMORY::memo.Detach( Vtot );
      MEMORY::memo.Detach( Vmax );
      MEMORY::memo.Detach( dzdt );
      MEMORY::memo.Detach( dz );
      MEMORY::memo.Detach( wght );
      MEMORY::memo.Detach( V );

    } break;

    //////////////////////////////////////////////////////////////////////////////////////
    // Finite-Volumes: node wise dzdt: flow through surrounding edges
    //////////////////////////////////////////////////////////////////////////////////////
//    case 2:
//    {
//      for( int e=0; e<rgne; e++ )
//      {
//        ELEM* elem = rg->Getelem(e);
//      }
//    } break;


    //////////////////////////////////////////////////////////////////////////////////////
    // Finite-Elements
    //////////////////////////////////////////////////////////////////////////////////////
    case 3:
    {
      // set up equation numbers and allocate memory -------------------------------------
//      project->fix[0]   = 0;
//      project->elemKind = ELEM::kBound | ELEM::kRegion;
      project->fix[0]   = BCON::kInlet;
      project->elemKind = ELEM::kRegion;

      SetEqno( model, 1, 0, 0, project->fix, project->elemKind );

      if( shape == kLinear )  coefs = kBottomEvol_L;
      else                    coefs = kBottomEvol;

      double* B = (double*) MEMORY::memo.Array_eq( neq, "EQS_DZ::Execute(9)" ); // right hand side
      double* X = (double*) MEMORY::memo.Array_eq( neq, "EQS_DZ::Execute(10)" );// change of Z

      initStructure = true;

      for( int it=0; it<project->actualCycit; it++ )
      {
        for( int e=0; e<neq; e++ )  X[e] = 0.0;

        // solve equation system ---------------------------------------------------------
        diverged_cg = Solve( model, neq, B, X, project );

        // copy solution to dzdt ---------------------------------------------------------
        for( int n=0; n<rgnp; n++ )
        {
          NODE* nd = rg->Getnode(n);
          int   no = nd->Getno();

          int eqno = GetEqno( nd, 0 );

          if( eqno >= 0 )  dzdt[no] += X[eqno];
          else             dzdt[no]  = 0.0;
        }

        // compute morphological time step -----------------------------------------------
        dt = MorphTime( project, rgnp, dzdt );

        // account for rigid bed condition -----------------------------------------------
//        for( int n=0; n<rgnp; n++ )
//        {
//          NODE* nd = rg->Getnode(n);
//          int   no = nd->Getno();

//          if( isFS(nd->flag, NODE::kCornNode) )
//          {
//            if( dzdt[no] < 0.0 )     // in case of erosion ...
//            {
//              etaQb[no] = -sed->Erode(nd) / dzdt[no] / dt;
//              if( etaQb[no] > 1.0 )  etaQb[no] = 1.0;
//            }
//          }
//        }

        // statistics of correction vector -----------------------------------------------
        int    nc     = 0;
        int    noAbs  = -1;
        double maxAbs = 0.0;
        double aveAbs = 0.0;

        for( int n=0; n<rgnp; n++ )
        {
          NODE* nd = rg->Getnode(n);
          int eqno = GetEqno( nd, 0 );

          if( eqno >= 0 )
          {
            nc++;
            aveAbs += X[eqno];

            if( noAbs < 0 )
            {
              noAbs  = nd->Getname();
              maxAbs = fabs( X[eqno] );
            }
            else if( fabs(X[eqno]) > maxAbs )
            {
              noAbs  = nd->Getname();
              maxAbs = fabs( X[eqno] );
            }
          }
        }

        //////////////////////////////////////////////////////////////////////////////////
#       ifdef _MPI_
        maxAbs = subdom->Mpi_max( maxAbs );
        aveAbs = subdom->Mpi_sum( aveAbs );
        nc     = subdom->Mpi_sum( nc );
#       endif
        //////////////////////////////////////////////////////////////////////////////////

        if( nc )  aveAbs /= nc;

        sprintf( text, "\n\n%-25s%s\n\n %s\n\n",
                        " (EQS_DZ::Execute)", "convergence parameters ...",
                        " variable   node          average        maximum" );
        REPORT::rpt.Message( 1, text );

        sprintf( text, "     %4s   %7d     %12.5le   %12.5le %s\n",
                      "dzdt", noAbs, aveAbs, maxAbs, "     (abs)" );
        REPORT::rpt.Message( 1, text );

        // check for convergence ---------------------------------------------------------
        if( maxAbs < 1.0e-6 )  break;
      }

      MEMORY::memo.Detach( B );
      MEMORY::memo.Detach( X );

    } break;


    //////////////////////////////////////////////////////////////////////////////////////
    // choice not implemented
    //////////////////////////////////////////////////////////////////////////////////////

    default:
      sprintf( text, "\n\n%-25s%s\n\n",
                      " (EQS_DZ::Execute)", "no suitable choice for bed evolution" );
      REPORT::rpt.Warning( kParameterFault, text );
      break;

  } // end of switch


  ////////////////////////////////////////////////////////////////////////////////////////

  // change bottom elevation -------------------------------------------------------------
  for( int n=0; n<rgnp; n++ )
  {
    NODE* nd = rg->Getnode(n);
    int   no = nd->Getno();

    if( isFS(nd->flag, NODE::kCornNode) )
    {
      nd->zor += dzdt[no] * dt;

      if( nd->zor < nd->zero )
      {
        dzdt[no] = 0.0;
        nd->zor = nd->zero;
      }

      nd->dz += dzdt[no] * dt;
      nd->z   = nd->zor;
    }
  }

  // set NODE::qb[] ----------------------------------------------------------------------
  for( int n=0; n<rgnp; n++ )
  {
    NODE* nd = rg->Getnode(n);
    int   no = nd->Getno();

    nd->v.Qb *= etaQb[no];
  }

  // interpolate bed elevation at midside nodes ------------------------------------------
  for( int e=0; e<rg->Getne(); e++ )
  {
    ELEM* el = rg->Getelem(e);

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    for( int i=ncn; i<nnd; i++ )
    {
      int il, ir;
      el->GetQShape()->getCornerNodes( i, &il, &ir );

      NODE* ndm = el->Getnode(i);
      NODE* ndl = el->Getnode(il);
      NODE* ndr = el->Getnode(ir);

      ndm->zor = 0.5 * (ndl->zor + ndr->zor);
      ndm->z   = 0.5 * (ndl->z + ndr->z);
      ndm->dz  = 0.5 * (ndl->dz + ndr->dz);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  MEMORY::memo.Detach( dzdt );
  MEMORY::memo.Detach( etaQb );

  REPORT::rpt.Message( 1, "\n%-25s%s %12.3lf [s]\n",
                          " (EQS_DZ::Execute)", "morphological time step =",
                          project->timeint.incTime.Getsec() );
}

//////////////////////////////////////////////////////////////////////////////////////////
/*
void EQS_DZ::FV( PROJECT* project )
{
  // loop on all elements ----------------------------------------------------------------

  // -------------------------------------------------------------------------------------
  // initializations

  SED*   sed = &project->sed;
  double por = 1.0 - project->sed.por;
  double area = elem->area();

  SHAPE* lShape = elem->getLShape();

  int ncn = lShape->nnd;

  if( force )  for( int i=0; i<maxEleq; i++ )  force[i] = 0.0;

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      double* estifmPtr = estifm[i];
      for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


  // compute local coordinates relative to first node ------------------------------------

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<ncn; i++ )
  {
    x[i]   = elem->nd[i]->x - *x;
    y[i]   = elem->nd[i]->y - *y;
  }

  x[0] = y[0] = 0.0;


  double adz[kMaxNodes2D], dzdt[kMaxNodes2D];

  for( int i=0; i<ncn; i++ )
  {
    int no  = elem->nd[i]->Getno();
    adz[i]  = this->etaQb[no];
    dzdt[i] = this->dzdt[no];
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // compute transport qb through element edges (Finite-Volumes)

  SHAPE* lineShape = SHAPE::get( kLine, 2 );

  switch( ncn )
  {
    case 3:
      for( int i=0; i<ncn; i++ )
      {
        int j = (i+1) % ncn;
        int k = (i+2) % ncn;

        NODE*  ndi = elem->nd[i];
        NODE*  ndj = elem->nd[j];

        int    noi = ndi->Getno();
        int    noj = ndj->Getno();

        double qbi = ndi->v.qb;
        double sxi = sed->sx[noi];
        double syi = sed->sy[noi];

        double qbj = ndj->v.qb;
        double sxj = sed->sx[noj];
        double syj = sed->sy[noj];

        for( int g=0; g<lineShape->ngp; g++ )
        {
          double* m      = lineShape->f[g];
          double* dm     = lineShape->dfdx[g];
          double  weight = lineShape->weight[g];

          double nx =  dm[0]*y[i] + dm[1]*y[j];
          double ny = -dm[0]*x[i] - dm[1]*x[j];

          double qb = m[0] * qbi  +  m[1] * qbj;
          double sx = m[0] * sxi  +  m[1] * sxj;
          double sy = m[0] * syi  +  m[1] * syj;

          force[k] -= weight * adz[k] * qb * (sx*nx + sy*ny);         // RHS
        }

        force[k] -= area * por * dzdt[k];                             // RHS

        estifm[k][k] += area * por;
      }
      break;

    case 4:
      break;
}
*/
//////////////////////////////////////////////////////////////////////////////////////////

double EQS_DZ::MorphTime( PROJECT* project, int n, double* dz )
{
  TIMEINT* tmint  = &project->timeint;
  SUBDOM*  subdom = &project->subdom;
  MODEL*   model  = project->M2D;
  GRID*    rg     = model->region;

  // determine maximum change in bottom elevation max( dz[] ) ----------------------------
  double maxdzdt = 0.0;

  for( int i=0; i<n; i++ )
  {
    NODE* nd = rg->Getnode(i);

    if( isFS(nd->flag, NODE::kCornNode) )
    {
      if( fabs(dz[i]) > maxdzdt )  maxdzdt = fabs(dz[i]);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////
# ifdef _MPI_
  maxdzdt = subdom->Mpi_max( maxdzdt );
# endif
  ////////////////////////////////////////////////////////////////////////////////////////

  // determine morphological time step ---------------------------------------------------
  double dt = project->sed.maxDz / maxdzdt;

  if( dt < tmint->relaxTimeFlow.Getsec() )  dt = tmint->relaxTimeFlow.Getsec();
  if( dt > 1.5 * tmint->incTime.Getsec() )  dt = 1.5 * tmint->incTime.Getsec();

  tmint->incTime.Setsec( dt );

  if( tmint->actualTime + (tmint->incTime * 1.1) > tmint->nextTime )
  {
    tmint->incTime = tmint->nextTime - tmint->actualTime;
    dt = tmint->incTime.Getsec();
  }
  else if( tmint->incTime > (tmint->nextTime - tmint->actualTime) * 0.45 )
  {
    tmint->incTime = (tmint->nextTime - tmint->actualTime) * 0.5;
    dt = tmint->incTime.Getsec();
  }
  else if( tmint->incTime.Getsec() <= 1.0e-6 )
  {
    tmint->incTime.Set( "00:00:00.000001" );
    dt = tmint->incTime.Getsec();
  }

  return dt;
}


// ---------------------------------------------------------------------------------------

void EQS_DZ::QbDiff( PROJECT* project, int shape )
{
  MODEL*   model   = project->M2D;

  GRID*    rg      = model->region;
  int      rgnp    = rg->Getnp();
  int      rgne    = rg->Getne();

  SED*     sed     = &project->sed;

  SUBDOM*  subdom  = &project->subdom;


  // print information on actual iteration -----------------------------------------------

  project->PrintTheCycle( 1 );
  REPORT::rpt.PrintTime( 1 );


  ////////////////////////////////////////////////////////////////////////////////////////
  // compute bottom evolution from difference:  qb - qbz

  // statistis: maximum and averaged change in qb ----------------------------------------

  double minQb = 0.0;
  double maxQb = 0.0;
  double aveQb = 0.0;

  double minDz = 0.0;
  double maxDz = 0.0;
  double aveDz = 0.0;

  int nptot = 0;

  for( int n=0; n<rgnp; n++ )
  {
    NODE* nd = rg->Getnode(n);
    int   no = nd->Getno();

    if( isFS(nd->flag, NODE::kCornNode) && !isFS(nd->flag, NODE::kDry) )
    {
      double delQb = nd->v.Qb - nd->qbo;
      double delDz = -fabs(nd->v.Qb) * delQb;

      aveQb += delQb;
      if( delQb < minQb )  minQb = delQb;
      if( delQb > maxQb )  maxQb = delQb;

      if( delDz > 0.0  ||  nd->zor > nd->zero + 0.5 * sed->maxDz )
      {
        aveDz += delDz;
        if( delDz < minDz )  minDz = delDz;
        if( delDz > maxDz )  maxDz = delDz;
      }

      if( !isFS(nd->flag, NODE::kInface_DN) )  nptot++;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
# ifdef _MPI_
  minQb = subdom->Mpi_min( minQb );
  maxQb = subdom->Mpi_max( maxQb );
  aveQb = subdom->Mpi_sum( aveQb );

  minDz = subdom->Mpi_min( minDz );
  maxDz = subdom->Mpi_max( maxDz );
  aveDz = subdom->Mpi_sum( aveDz );

  nptot = subdom->Mpi_sum( nptot );
# endif
  ////////////////////////////////////////////////////////////////////////////////////////

  if( nptot > 0 )
  {
    aveQb /= nptot;
    aveDz /= nptot;
  }

  // statistics of correction vector -----------------------------------------------------

  char text[500];

  sprintf( text, "\n\n%-25s%s\n\n%-25s  %12s   %12s   %12s\n",
                  " (EQS_DZ::QbDiff)", "convergence parameters ...",
                  " ",                 "average   ", "minimum   ", "maximum   " );
  REPORT::rpt.Message( 1, text );

  sprintf( text, "%25s  %12.4le   %12.4le   %12.4le\n",
                 "delta_qb  ", aveQb, minQb, maxQb );
  REPORT::rpt.Message( 1, text );

  sprintf( text, "%25s  %12.4le   %12.4le   %12.4le\n",
                 "qb * delta_qb  ", aveDz, minDz, maxDz );
  REPORT::rpt.Message( 1, text );


  // change bottom elevation -------------------------------------------------------------

  if( fabs(minDz) > fabs(maxDz) )  maxDz = minDz;

  if( initDzQb  ||  fabs(factDzQb*maxDz) > sed->maxDz )
  {
    initDzQb = false;
    factDzQb = sed->maxDz / fabs(maxDz);
  }

  sprintf( text, "\n\n%-25s%s %12.5le\n\n",
                  " (EQS_DZ::QbDiff)", "maximum change in bottom elevation =",
                  factDzQb * maxDz );
  REPORT::rpt.Message( 1, text );


  for( int n=0; n<rgnp; n++ )
  {
    NODE* nd = rg->Getnode(n);
    int   no = nd->Getno();

    if( isFS(nd->flag,NODE::kCornNode) && !isFS(nd->flag,NODE::kDry) )
    {
      double dz = -factDzQb * fabs(nd->v.Qb) * (nd->v.Qb - nd->qbo);

      int evol = true;

      // no erosion/deposition on inlet boundaries with equilibrium condition ------------
      if( isFS(nd->bc.kind, BCON::kRateC)
          &&  nd->bc.val->Qb[0] < 0.0 )                     evol = false;

      // no erosion of points with steep sloping bed  ------------------------------------
      if( isFS(sed->slope,SED::kSLOPE_Max)  &&  dz < 0.0
          &&  sed->dzmx[no] > sed->maxSlope )               evol = false;

      // no erosion of unerodible horizont -----------------------------------------------
      if( nd->zor < nd->zero+maxDz/100.0  &&  dz < 0.0 )    evol = false;

      // no bottom evolution on dry nodes (flag NODE::kMarsh) ----------------------------
      if( isFS(nd->flag, NODE::kMarsh) )                    evol = false;


      if( evol )
      {
        nd->zor += dz;

        if( nd->zor < nd->zero )
        {
          dz = 0.0;
          nd->zor = nd->zero;
        }

        nd->z = nd->zor;
        nd->dz += dz;
      }
    }
  }


  // interpolate bottom elevation at midside nodes ---------------------------------------

  for( int e=0; e<rg->Getne(); e++ )
  {
    ELEM* el = rg->Getelem(e);

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    for( int i=ncn; i<nnd; i++ )
    {
      int il, ir;
      el->GetQShape()->getCornerNodes( i, &il, &ir );

      NODE* ndm = el->Getnode(i);
      NODE* ndl = el->Getnode(il);
      NODE* ndr = el->Getnode(ir);

      ndm->zor = 0.5 * (ndl->zor + ndr->zor);
      ndm->z   = 0.5 * (ndl->z + ndr->z);
    }
  }
}
