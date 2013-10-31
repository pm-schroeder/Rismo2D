// ======================================================================================
//
// Copyright (C) 1992-2007  by  P.M. SCHROEDER
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
// ---------------------------------------------------------------------------------------
// The module init.c contains several initialization procedures:
//
//    void initKD()            initialize values of K and D
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Type.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Project.h"

#include "Grid.h"
#include "Model.h"


//#define DEBUG


// ---------------------------------------------------------------------------------------
// Function to initialize values for the turbulent kinetic energy (K) and
// the dissipation of K (D) for all nodes of the flow domain.
//
// In the 2D depth-averaged case K and D are approximated with an
// algebraic relation from the flow depth and the bottom shear stresses.
//
//
// Michael Schroeder, November 1992
//                    January  1994
// ---------------------------------------------------------------------------------------

void GRID::InitKD( PROJECT* project )
{
  char   text[200];

  int    turb = project->actualTurb;

  KDCONST* KD = &project->KD;

  double cm   = KD->cm;             // constants of the k-epsilon model
  double cd   = KD->cd;
  double c1D  = KD->c1D;
  double c2D  = KD->c2D;


  ////////////////////////////////////////////////////////////////////////////////////////
  // initialize K,D with algebraic model

  if( isFS(turb, BCONSET::kVtAlgebraic) )
  {
    sprintf( text, "\n (GRID::InitKD)          %s\n",
                   "initializing K,D with algebraic model" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt   = (int*)    MEMORY::memo.Array_nd( np );
    double* stEst = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n]   = 0;
      stEst[n] = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          int no = nd[i]->Getno();
          stEst[no] += type->KDest * type->st;
          cnt[no]++;
        }
      }
    }

    for( int n=0; n<np; n++ )
    {
      if( isFS(node[n].bc.kind, BCON::kSetKD) )  continue;

      if( cnt[n] )
      {
        stEst[n] /= cnt[n];

        double U = node[n].v.U;
        double V = node[n].v.V;
        double H = node[n].v.S - node[n].z;
        if( H < project->hmin ) H = project->hmin;

        double Utau = sqrt( node[n].cf * (U*U + V*V) );

        double cK = 1.0 / sqrt( node[n].cf );

        node[n].v.K = sqrt( cK * stEst[n] / cm / cd ) * Utau*Utau;
        node[n].v.D = cK * Utau*Utau*Utau / H;
      }
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( stEst );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // experimental: initialize K,D with algebraic shear stress model

  else if( isFS(turb, BCONSET::kVtAlgShear) )
  {
    sprintf( text, "\n (GRID::InitKD)          %s\n",
                   "initializing K,D with algebraic shear stress model" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt   = (int*) MEMORY::memo.Array_nd( np );
    double* stEst = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n]   = 0;
      stEst[n] = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          int no = nd[i]->Getno();
          stEst[no] += type->KDest * type->st;
          cnt[no]++;
        }
      }
    }

    double* phi = project->M2D->Phi2D();

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        stEst[n] /= cnt[n];

        double U  = node[n].v.U;
        double V  = node[n].v.V;
        double H  = node[n].v.S - node[n].z;
        double cf = node[n].cf;

        double Utau = sqrt( cf * (U*U + V*V) );

        double cK   = 1.0 / sqrt( cf );
        double PKv  = cK * Utau*Utau*Utau / H;

        double cD   = c2D * sqrt( cm/stEst[n] / sqrt(cf*cf*cf) );
        double PDv  = cD * (Utau*Utau/H) * (Utau*Utau/H);

        double alfa = PDv*PDv / phi[n] / cm / cd;
        double beta = c1D * PKv;
        double delc = c2D - c1D;

        double ac = delc * delc;
        double bc = 2.0 * delc * beta;
        double cc = beta * beta - alfa;
        double dc = alfa * PKv;

        double x[3];
        int real = project->Cube( ac, bc, cc, dc, x );

        if( real == 1 )
        {
          node[n].v.D = x[0];
        }
        else
        {
          node[n].v.D = x[2];

          if( node[n].v.D < PKv )
          {
            node[n].v.D = x[0];
            if( node[n].v.D < PKv )  node[n].v.D = x[1];
          }
        }

        if( node[n].v.D < PKv )
        {
          node[n].v.K = sqrt( cK * stEst[n] / cm / cd ) * Utau*Utau;
          node[n].v.D = cK * Utau*Utau*Utau / H;
        }
        else
        {
          //node[n].v.K = sqrt( node[n].v.D * (node[n].v.D-PKv) / phi[n] / cm / cd );
          node[n].v.K = (c2D*node[n].v.D - c1D*(node[n].v.D-PKv)) * node[n].v.D / PDv;
        }
      }
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( stEst );
    MEMORY::memo.Detach( phi );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // initialize K,D with algebraic shear stress model

  //else if( isFS(turb, BCONSET::kVtAlgShear) )
  //{
  //  sprintf( text, "\n (GRID::InitKD)          %s\n",
  //                 "initializing K,D with algebraic shear stress model" );
  //  REPORT::rpt.Output( text, 4 );

  //  int*    cnt   = (int*)    MEMORY::memo.Array_nd( np );
  //  double* lm    = (double*) MEMORY::memo.Array_nd( np );
  //  double* stEst = (double*) MEMORY::memo.Array_nd( np );

  //  for( int n=0; n<np; n++ )
  //  {
  //    cnt[n]   = 0;
  //    lm[n]    = 0;
  //    stEst[n] = 0.0;
  //  }

  //  for( int e=0; e<ne; e++ )
  //  {
  //    ELEM* el = &elem[e];

  //    if( !isFS(el->flag, ELEM::kDry) )
  //    {
  //      TYPE* type = TYPE::Getid( el->type );

  //      NODE** nd  = el->nd;
  //      int    nnd = el->getnnd();

  //      for( int i=0; i<nnd; i++ )
  //      {
  //        int no = nd[i]->Getno();
  //        lm[no]    += type->lm;
  //        stEst[no] += type->KDest * type->st;
  //        cnt[no]++;
  //      }
  //    }
  //  }

  //  double* phi = project->M2D->Phi2D();

  //  for( int n=0; n<np; n++ )
  //  {
  //    if( cnt[n] )
  //    {
  //      stEst[n] /= cnt[n];
  //      lm[n]    /= cnt[n];

  //      double U  = node[n].v.U;
  //      double V  = node[n].v.V;
  //      double H  = node[n].v.S - node[n].z;
  //      double cf = node[n].cf;

  //      double Utau = sqrt( cf * (U*U + V*V) );

  //      double cK   = 1.0 / sqrt( cf );
  //      double PKv  = cK * Utau*Utau*Utau / H;

  //      double lmH  = lm[n];
  //      double lmH4 = lmH*lmH*lmH*lmH;

  //      double ac =  1.0;
  //      double bc =  0.0;
  //      double cc = -lmH4 * phi[n];
  //      double dc = -lmH4 * PKv;

  //      double vt, x[3];
  //      int real = project->Cube( ac, bc, cc, dc, x );

  //      if( real == 1 )
  //      {
  //        vt = x[0];
  //        if( vt < project->minVt )  vt = project->minVt;
  //      }
  //      else
  //      {
  //        vt = x[0];

  //        if( vt < 0.0 )
  //        {
  //          if( x[1] > 0.0 )  vt = x[1];
  //          if( x[2] > 0.0 )  vt = x[2];
  //        }

  //        if( vt < project->minVt )  vt = project->minVt;
  //      }

  //      node[n].v.K = vt * vt / lmH / lmH / sqrt( cm * cd );
  //      node[n].v.D = vt * vt * vt / lmH4;
  //      node[n].vt  = vt;
  //    }
  //  }

  //  MEMORY::memo.Detach( cnt );
  //  MEMORY::memo.Detach( lm );
  //  MEMORY::memo.Detach( stEst );
  //  MEMORY::memo.Detach( phi );
  //}


  ////////////////////////////////////////////////////////////////////////////////////////
  // Prandtl's mixing length model

  else if( isFS(turb, BCONSET::kVtMixingLength) )
  {
    sprintf( text, "\n (GRID::InitKD)          %s\n",
                   "initializing K,D with Prandtl's mixing length model" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
    double* lm  = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n] = 0;
      lm[n]  = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          int no = nd[i]->Getno();
          lm[no] += type->lm;
          cnt[no]++;
        }
      }
    }

    double* phi = project->M2D->Phi2D();

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        double H  = node[n].v.S - node[n].z;

        double lmH = H * lm[n]/cnt[n];

        double vt  = lmH * lmH * sqrt( phi[n] );
        double tmp = vt / lmH;

        node[n].v.K = tmp * tmp / sqrt( cm*cd );
        node[n].v.D = tmp * tmp * tmp / lmH;
      }
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( lm );
    MEMORY::memo.Detach( phi );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // large eddy model

  else if( isFS(turb, BCONSET::kVtLES) )
  {
    sprintf( text, "\n (GRID::InitKD)          %s\n",
                   "initializing K,D with large eddy model" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
    double* lm  = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n] = 0;
      lm[n]  = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        double xmin = nd[0]->x;
        double xmax = nd[0]->x;
        double ymin = nd[0]->y;
        double ymax = nd[0]->y;

        for( int i=1; i<nnd; i++ )
        {
          if( nd[i]->x < xmin )  xmin = nd[i]->x;
          if( nd[i]->x > xmax )  xmax = nd[i]->x;

          if( nd[i]->y < ymin )  ymin = nd[i]->y;
          if( nd[i]->y > ymax )  ymax = nd[i]->y;
        }

        double dl = sqrt( (xmax-xmin)*(ymax-ymin) );

        for( int i=0; i<nnd; i++ )
        {
          int no = nd[i]->Getno();
          lm[no] += type->lm * dl;
          cnt[no]++;
        }
      }
    }

    double* phi = project->M2D->Phi2D();

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        lm[n] /= cnt[n];

        double vt  = lm[n] * lm[n] * sqrt( phi[n] );
        double tmp = vt / lm[n];

        node[n].v.K = tmp * tmp / sqrt( cm*cd );
        node[n].v.D = tmp * tmp * tmp / lm[n];
      }
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( lm );
    MEMORY::memo.Detach( phi );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // initialize K,D with minimum values

  else
  {
    sprintf( text, "\n (GRID::InitKD)          %s\n",
                   "initializing K,D with minimum values" );
    REPORT::rpt.Output( text, 4 );

    for( int n=0; n<np; n++ )
    {
      node[n].v.K = project->minK;
      node[n].v.D = project->minD;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
//  minVt will be hold only in EQS::RegionUVS2D(), ... 02-03-2010, SC
//  if( isFS(turb, BCONSET::kVtMin) )
//  {
//    for( int n=0; n<np; n++ )
//    {
//      if( node[n].v.K < project->minK )  node[n].v.K = project->minK;
//      if( node[n].v.D < project->minD )  node[n].v.D = project->minD;
//    }
//  }
}
