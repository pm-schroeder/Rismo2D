// =================================================================================================
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
// =================================================================================================
//
#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Memory.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsDisp.h"


EQS_DISP::EQS_DISP() : EQS( 1, 1, 0 )
{
  neq  = 0;
}


EQS_DISP::~EQS_DISP()
{
}


// -------------------------------------------------------------------------------------------------

void EQS_DISP::Execute( PROJECT* project )
{
  MODEL* model = project->M2D;
  GRID*  rg    = model->region;
  int    np    = rg->Getnp();
  int    ne    = rg->Getne();

  int    diverged_cg = 0;

  model->Incinit();

  initStructure = true;

  project->PrintTheCycle( 1 );
  REPORT::rpt.PrintTime( 1 );

  // set up equation numbers -----------------------------------------------------------------------
  project->fix[0]   = BCON::kFixV;
  project->elemKind = ELEM::kRegion;

  //SetEqno( model, 1, 1, 0, project->fix, project->elemKind );
  SetEqno( model, 1, 0, 0, project->fix, project->elemKind );


  // -----------------------------------------------------------------------------------------------
  // compute secondary flow from curvature of stream lines

  double* rho = project->M2D->Curv2D();                // curvatur of streamlines
  double* wgt = (double*) MEMORY::memo.Array_nd(np);

  double minUSf   = project->minUSf;
  double maxTanSf = project->maxTanSf;
  double kappa    = project->kappa;

  for( int n=0; n<np; n++ )
  {
    NODE* nd = rg->Getnode(n);
    nd->Vsec = 0.0;
    wgt[n]   = 0.0;
  }

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = rg->Getelem(e);

    TYPE*  type = TYPE::Getid( el->type );
    double beta = type->betaSf;
    double area = el->area();

    if( beta > 0.0 )
    {
      for( int i=0; i<el->Getnnd(); i++ )
      {
        NODE* nd = el->Getnode(i);
        int   no = nd->Getno();

        if( GetEqno(nd,0) >= 0 )
        {
          double mm   = kappa / sqrt( nd->cf );
          double alfa = beta * (mm + 0.5) / kappa / kappa / mm;

          double H = nd->v.S - nd->z;
          if( H < 0.0 )  H = 0.0;

          double U  = nd->v.U;
          double V  = nd->v.V;
          double Us = sqrt( U*U + V*V );

          double tanSf = 0.0;
          if( Us > minUSf ) tanSf = alfa * H * rho[no];

          // limit maximum secondary velocity --------------------------------------------------------
          if( fabs(tanSf) > maxTanSf ) tanSf *= maxTanSf / fabs(tanSf);

          nd->Vsec += tanSf * Us;
          wgt[no]  += area;
        }
        else
        {
          nd->Vsec = 0.0;
        }
      }
    }
  }

  for( int n=0; n<np; n++ )
  {
    NODE* nd = rg->Getnode(n);
    if( wgt[n] > 0.0 ) nd->Vsec /= wgt[n];
  }

  MEMORY::memo.Detach( wgt );
  MEMORY::memo.Detach( rho );


  // -----------------------------------------------------------------------------------------------
  // solve differential equation for secondary flow (depending on parameter mueSf)

  if( project->mueSf > 0.0 )
  {
    double* X = (double*) MEMORY::memo.Array_eq(neq);
    double* B = (double*) MEMORY::memo.Array_eq(neq);

    for( int n=0; n<neq; n++ )  X[n] = B[n] = 0.0;

    for( int n=0; n<np; n++ )
    {
      NODE* nd = rg->Getnode(n);
      int   en = GetEqno( nd, 0 );

      if( en >= 0 ) X[en] = nd->Vsec;
    }

    // solve equations -----------------------------------------------------------------------------
    diverged_cg = Solve( model, neq, B, X, project );

    // compute secondary flow velocity -------------------------------------------------------------
    for( int n=0; n<np; n++ )
    {
      NODE* nd = rg->Getnode(n);
      int   en = GetEqno( nd, 0 );

      if( en >= 0 )
      {
        double U  = nd->v.U;
        double V  = nd->v.V;
        double Us = sqrt( U*U + V*V );

        double tanSf = 0.0;
        if( Us > minUSf ) tanSf = X[en] / Us;

        if( fabs(tanSf) > maxTanSf ) nd->Vsec = Us * maxTanSf * tanSf / fabs(tanSf);
        else                         nd->Vsec = Us * tanSf;
      }
      else
      {
        nd->Vsec = 0.0;
      }
    }

    MEMORY::memo.Detach( B );
    MEMORY::memo.Detach( X );
  }

  // smooth secondary flow and  force betaSf -------------------------------------------------------
  double *V = (double*) MEMORY::memo.Array_eq(np);
  double *W = (double*) MEMORY::memo.Array_eq(np);

  int smoothpasses = 5;

  for( int pass=1; pass<=smoothpasses; pass++ )
  {
    for( int n=0; n<np; n++ ) V[n] = W[n] = 0.0;

    for( int e=0; e<rg->Getne(); e++ )
    {
      ELEM  *el   = rg->Getelem(e);
      TYPE  *type = TYPE::Getid( el->type );

      if( type->betaSf > 0.0 )
      {
        double area = el->area();
        int    ncn  = el->Getncn();
        double Vsec = 0.0;
        int    cnt  = 0;

        for( int i=0; i<ncn; i++ )
        {
          NODE *nd = el->Getnode(i);

          if( GetEqno(nd,0) >= 0 )
          {
            Vsec += nd->Vsec;
            cnt++;
          }
        }
        if( cnt > 0 ) Vsec /= cnt;

        for( int i=0; i<ncn; i++ )
        {
          NODE *nd = el->Getnode(i);

          int   no = nd->Getno();
          V[no] += area * Vsec;
          W[no] += area;
        }
      }
    }

    for( int n=0; n<np; n++ )
    {
      NODE *nd = rg->Getnode(n);

      if( GetEqno(nd,0) >= 0  &&  W[n] > 0.0 ) nd->Vsec = V[n] / W[n];
      else                                     nd->Vsec = 0.0;
    }
  }

  MEMORY::memo.Detach( V );
  MEMORY::memo.Detach( W );

  // compute secondary flow at midside nodes -------------------------------------------------------
  for( int e=0; e<rg->Getne(); e++ )
  {
    ELEM *el  = rg->Getelem(e);
    int   ncn = el->Getncn();
    int   nnd = el->Getnnd();

    for( int i=ncn; i<nnd; i++ )
    {
      // get left and right corner node to midside node i ------------------------------------------
      int    il, ir;
      double left, rght;

      el->GetQShape()->getCornerNodes( i, &il, &ir );

      left = el->nd[il]->Vsec;
      rght = el->nd[ir]->Vsec;
      el->nd[i]->Vsec = 0.5 * (left + rght);
    }
  }

  // compute dispersion coefficients ---------------------------------------------------------------
  for( int n=0; n<np; n++ )
  {
    NODE* nd = rg->Getnode(n);

    double U  = nd->v.U;
    double V  = nd->v.V;
    double H  = nd->v.S - nd->z;
    double Us = sqrt( U*U + V*V );

    double cf = nd->cf;

    if( Us > minUSf )
    {
      // ===========================================================================================
      // changed on 07.04.2007
      // dispersion computed in relation to H*Us*Us

      //double mm   = project->kappa / sqrt( cf );
      //double aHUR = nd->Vsec;

      //double Dxx = Us * Us / mm / mm;
      //double Dyy = aHUR * aHUR / 3.0;
      //double Dxy = aHUR * Us / (2.0*mm + 1.0);

      //// rotation to mean flow direction
      //double sina =  V / Us;
      //double cosa =  U / Us;

      //nd->Duu = cosa*cosa*Dxx - 2.0*sina*cosa*Dxy + sina*sina*Dyy;
      //nd->Dvv = sina*sina*Dxx + 2.0*sina*cosa*Dxy + cosa*cosa*Dyy;
      //nd->Duv = sina*cosa*(Dxx-Dyy) + (cosa*cosa - sina*sina)*Dxy;

      double kcf = project->kappa / sqrt( cf );
      double rvs = nd->Vsec / Us;

      nd->Dxx = 1.0 / kcf / kcf;
      nd->Dyy = rvs * rvs / 3.0;
      nd->Dxy = rvs / (2.0*kcf + 1.0);

      // ===========================================================================================
    }
    else
    {
      nd->Dxx = 0.0;
      nd->Dxy = 0.0;
      nd->Dyy = 0.0;
    }
  }

  char text[80];

  sprintf( text, "\n (EQS_DISP::Execute)     %s\n",
                 "dispersion coefficients determined" );
  REPORT::rpt.Output( text, 4 );

  if( diverged_cg )  project->errLevel |= diverged_cg | kErr_no_conv_cg;
}
