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

// #define DEBUG

#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsKD2D.h"


int EQS_KD2D::Coefs( ELEM*    elem,
                     PROJECT* project,
                     double** estifm,
                     double*  force )
{
  // -------------------------------------------------------------------------------------
  // call the base class to perform some initializations

  EQS::Coefs( elem, project, estifm, force );


  // -------------------------------------------------------------------------------------
  // nothing to do with dry elements and boundary elements

  if(    isFS(elem->flag, ELEM::kDry)
      || isFS(elem->flag, ELEM::kBound) ) return 0;


  // -------------------------------------------------------------------------------------

  if( linearShape )
  {
    if( isFS(project->actualTurb, BCONSET::kVtAnisotrop) )
      RegionLAI( elem, project, estifm, force );

    else
      RegionL( elem, project, estifm, force );
  }

  else if( quarterShape )
  {
    // number of corner nodes
    int ncn = elem->Getncn();

    // quart is a virtual quarter of the element: copy parameters from elem to quart
    ELEM quart = *elem;

    // set up relation between the element nodes and its quarter nodes
    int relate[4][4];
    if( ncn == 3 )
    {
      relate[0][0] = 0; relate[0][1] = 3; relate[0][2] = 5;
      relate[1][0] = 1; relate[1][1] = 4; relate[1][2] = 3;
      relate[2][0] = 2; relate[2][1] = 5; relate[2][2] = 4;
      relate[3][0] = 3; relate[3][1] = 4; relate[3][2] = 5;
    }
    else
    {
      relate[0][0] = 7; relate[0][1] = 0; relate[0][2] = 4; relate[0][3] = 8;
      relate[1][0] = 4; relate[1][1] = 1; relate[1][2] = 5; relate[1][3] = 8;
      relate[2][0] = 5; relate[2][1] = 2; relate[2][2] = 6; relate[2][3] = 8;
      relate[3][0] = 6; relate[3][1] = 3; relate[3][2] = 7; relate[3][3] = 8;
    }

    // loop on 4 quarter elements: setup estifm and force
    if( ncn == 3 )
    {
      int Keq[3];       // equation index for K
      int Deq[3];       // equation index for D

      for( int i=0; i<4; i++ )
      {
        quart.nd[0] = elem->nd[relate[i][0]];
        quart.nd[1] = elem->nd[relate[i][1]];
        quart.nd[2] = elem->nd[relate[i][2]];

        Keq[0] = relate[i][0];   Deq[0] = relate[i][0] + 6;
        Keq[1] = relate[i][1];   Deq[1] = relate[i][1] + 6;
        Keq[2] = relate[i][2];   Deq[2] = relate[i][2] + 6;

        RegionQ( &quart, project, estifm, force, Keq, Deq );
      }
    }
    else
    {
      int Keq[4];       // equation index for K
      int Deq[4];       // equation index for D

      for( int i=0; i<4; i++ )
      {
        quart.nd[0] = elem->nd[relate[i][0]];
        quart.nd[1] = elem->nd[relate[i][1]];
        quart.nd[2] = elem->nd[relate[i][2]];
        quart.nd[3] = cent[elem->Getno()];

        Keq[0] = relate[i][0];   Deq[0] = relate[i][0] + 8;
        Keq[1] = relate[i][1];   Deq[1] = relate[i][1] + 8;
        Keq[2] = relate[i][2];   Deq[2] = relate[i][2] + 8;
        Keq[3] = 16;             Deq[3] = 17;              // center node equation

        RegionQ( &quart, project, estifm, force, Keq, Deq );
      }

      // eliminate the center node equations 17 and 16 from the equation system
      double pivot;

      pivot = estifm[17][17];
      if( fabs(pivot) < 1.0e-30 )
        REPORT::rpt.Error( kParameterFault,
                           "pivot too small (EQS_KD2D::Coefs - 1)" );

      for( int i=0; i<17; i++ )
      {
        double factor = estifm[i][17] / pivot;

        force[i] -= factor * force[17];

        for( int j=0; j<17; j++ )
        {
          estifm[i][j] -= factor * estifm[17][j];
        }
      }

      pivot = estifm[16][16];

      if( fabs(pivot) < 1.0e-30 )
        REPORT::rpt.Error( kParameterFault,
                           "pivot too small (EQS_KD2D::Coefs - 2)" );

      for( int i=0; i<16; i++ )
      {
        double factor = estifm[i][16] / pivot;

        force[i] -= factor * force[16];

        for( int j=0; j<16; j++ )
        {
          estifm[i][j] -= factor * estifm[16][j];
        }
      }
    }
  }

  else
  {
    if( isFS(project->actualTurb, BCONSET::kVtAnisotrop) )
      RegionAI( elem, project, estifm, force );

    else
      Region( elem, project, estifm, force );
  }

  return 1;
}


void EQS_KD2D::Region( ELEM*    elem,
                       PROJECT* project,
                       double** estifm,
                       double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  TYPE* type = TYPE::Getid( elem->type );

  double cmd = project->KD.cm * project->KD.cd;   // constants of k-epsilon model
  double sK  = project->KD.sK;
  double sD  = project->KD.sD;
  double c1D = project->KD.c1D;
  double c2D = project->KD.c2D;
  double est = type->KDest;
  double st  = type->st;

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int ngp = qShape->ngp;                          // number of GAUSS points
  int nnd = qShape->nnd;                          // total number of nodes
  int ncn = lShape->nnd;                          // number of corner nodes

  int startD = nnd;


  // kinematic viscosity -----------------------------------------------------------------

  double vk = project->vk;


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve k-epsilon equations

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double* dfdxPtr = qShape->dfdx[g];
    double* dfdyPtr = qShape->dfdy[g];

    double trafo[2][2];

    double detj   = qShape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );
    double weight = detj * qShape->weight[g];


    // -------------------------------------------------------------------------------------
    // compute values of quadratic shape functions at GP g

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = qShape->f[g];

    for( int j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // ------------------------------------------------------------------------------------
    // values of linear shape functions at GP g

    double  dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    dfdxPtr = lShape->dfdx[g];
    dfdyPtr = lShape->dfdy[g];

    double* m = lShape->f[g];

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP i
    //             horizontal velocities   : U and V
    //             flow depth              : H
    //             turbulent kinetic energy: K
    //             dissipation of K        : D


    // integrate H, cf and vt with linear shape ------------------------------------------

    double H    = 0.0;
    double cf   = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE *node = elem->nd[j];

      double ndH  = node->v.S - node->z;
      if( ndH <= 0.0 )  ndH = 0.0;

      H  += m[j] * ndH;
      cf += m[j] * node->cf;
    }


    // integrate U, V, K, D and vt with quadratic shape ----------------------------------

    double U, dUdx, dUdy;
    double V, dVdx, dVdy;
    double K, dKdx, dKdy, dKdt;
    double D, dDdx, dDdy, dDdt;

    U = dUdx = dUdy = 0.0;
    V = dVdx = dVdy = 0.0;
    K = dKdx = dKdy = dKdt = 0.0;
    D = dDdx = dDdy = dDdt = 0.0;

    for( int j=0; j<nnd; j++ )
    {
      NODE *node = elem->nd[j];

      double ndU  = node->v.U;
      double ndV  = node->v.V;

      U    +=    n[j] * ndU;
      dUdx += dndx[j] * ndU;
      dUdy += dndy[j] * ndU;

      V    +=    n[j] * ndV;
      dVdx += dndx[j] * ndV;
      dVdy += dndy[j] * ndV;

      double ndK = node->v.K - project->minK;
      double ndD = node->v.D - project->minD;

      K    +=    n[j] * ndK;
      dKdx += dndx[j] * ndK;
      dKdy += dndy[j] * ndK;

      D    +=    n[j] * ndD;
      dDdx += dndx[j] * ndD;
      dDdy += dndy[j] * ndD;

      dKdt +=    n[j] * node->v.dKdt;
      dDdt +=    n[j] * node->v.dDdt;
    }

    if( K < 0.0  ||  D < 0.0 )  K = D = 0.0;

    K += project->minK;
    D += project->minD;


    // -----------------------------------------------------------------------------------

    double vt = 0.0;

    for( int j=0; j<nnd; j++ )  vt += n[j] * elem->nd[j]->vt;

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < project->minVtKD )  vt = project->minVtKD;
    }


    // -----------------------------------------------------------------------------------
    // compute production and dissipation terms

    double Fi    = 2.0*dUdx*dUdx + 2.0*dVdy*dVdy + (dUdy+dVdx)*(dUdy+dVdx);

    double Us    = sqrt( U*U + V*V );
    double Utau  = sqrt( cf ) * Us;

    double cK    = 1.0 / sqrt( cf );
    double PKv   = 0.0;
    if( H > 0.0 )  PKv = cK * Utau*Utau*Utau / H;

//  double cD    = c2D * sqrt(cmd) / sqrt(st*est) / pow(cf, 0.75);
    double cD    = c2D * sqrt( cmd/est/st / sqrt(cf*cf*cf) );
    double PDv   = 0.0;
    if( H > 0.0 )  PDv = cD * (Utau*Utau/H) * (Utau*Utau/H);

    double PK    =  cmd * K * K * Fi;
    double dPKdK =  cmd * 2.0 * K * Fi;

//  double PD    =  c1D * cmd * K * K * Fi;
    double PD    =  c1D * PK;
//  double dPDdK =  c1D * cmd * 2.0 * K * Fi;
    double dPDdK =  c1D * dPKdK;

    double DK    =  D * D;
    double dDKdD =  2.0 * D;

//  double DD    =  c2D * D * D;
    double DD    =  c2D * DK;
//  double dDDdD =  c2D * 2.0 * D;
    double dDDdD =  c2D * dDKdD;


    // -----------------------------------------------------------------------------------
    // compute diffusion terms

    double vtsK = vt / sK  +  vk;
    double vtsD = vt / sD  +  vk;


    // -----------------------------------------------------------------------------------
    // compute K and D equation and coefficients of NEWTON-RAPHSON matrix

    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // ---------------------------------------------------------------------------------
      // compute K-equation

      f   = H * D * dKdt;                        // time
      f  += H * D * (U * dKdx  +  V * dKdy);     // advection

      fx  = H * D * vtsK * dKdx;                 // diffusion
      fy  = H * D * vtsK * dKdy;

      f  += H * dDdx * vtsK * dKdx;
      f  += H * dDdy * vtsK * dKdy;

      f  += H * (DK - PK - D*PKv);               // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute D-equation

      f   = H * K * dDdt;                        // time
      f  += H * K * (U * dDdx  +  V * dDdy);     // advection

      fx  = H * K * vtsD * dDdx;                 // diffusion
      fy  = H * K * vtsD * dDdy;

      f  += H * dKdx * vtsD * dDdx;
      f  += H * dKdy * vtsD * dDdy;

      f  += H * (DD - PD - K*PDv);               // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + startD;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f + dndx[j] * fx + dndy[j] * fy;
      }
    }


    if( estifm )
    {
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfx_, dfy_, dfxx, dfyy;

      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      // K-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * D * relaxThdt_KD;
      df_x  =  weight * H * D * U;
      df_y  =  weight * H * D * V;

      dfxx  =
      dfyy  =  weight * H * D * vtsK;

      df_x +=  weight * H * dDdx * vtsK;
      df_y +=  weight * H * dDdy * vtsK;

      df__ -=  weight * H * dPKdK;

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] =                 dfxx * dndx[j];
        ty[j] =                                    dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        double* estifmPtr = estifm[j];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // D-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * dKdt;
      df__ +=  weight * H * (U*dKdx + V*dKdy);

      dfx_  =  weight * H * vtsK * dKdx;
      dfy_  =  weight * H * vtsK * dKdy;

      df_x  =  weight * H * vtsK * dKdx;
      df_y  =  weight * H * vtsK * dKdy;

      df__ +=  weight * H * (dDKdD - PKv);

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        double* estifmPtr = estifm[j] + startD;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // K-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * dDdt;
      df__ +=  weight * H * (U * dDdx  +  V * dDdy);

      dfx_  =  weight * H * vtsD * dDdx;
      dfy_  =  weight * H * vtsD * dDdy;

      df_x  =  weight * H * vtsD * dDdx;
      df_y  =  weight * H * vtsD * dDdy;

      df__ -=  weight * H * (dPDdK + PDv);

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        double* estifmPtr = estifm[j + startD];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // D-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * K * relaxThdt_KD;
      df_x  =  weight * H * K * U;
      df_y  =  weight * H * K * V;

      dfxx  =
      dfyy  =  weight * H * K * vtsD;

      df_x +=  weight * H * dKdx * vtsD;
      df_y +=  weight * H * dKdy * vtsD;

      df__ +=  weight * H * dDDdD;

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] =                 dfxx * dndx[j];
        ty[j] =                                    dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        double* estifmPtr = estifm[j + startD] + startD;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void EQS_KD2D::RegionAI( ELEM*    elem,
                         PROJECT* project,
                         double** estifm,
                         double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  TYPE* type = TYPE::Getid( elem->type );

  double cmd = project->KD.cm * project->KD.cd;  // constants of k-epsilon model
  double sK  = project->KD.sK;
  double sD  = project->KD.sD;
  double c1D = project->KD.c1D;
  double c2D = project->KD.c2D;
  double est = type->KDest;
  double st  = type->st;

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int ngp = qShape->ngp;                   // number of GAUSS points
  int nnd = qShape->nnd;                   // total number of nodes
  int ncn = lShape->nnd;                   // number of corner nodes

  int startD = nnd;


  // kinematic viscosity -----------------------------------------------------------------

  double vk = project->vk;


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve k-epsilon equations

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double* dfdxPtr = qShape->dfdx[g];
    double* dfdyPtr = qShape->dfdy[g];

    double trafo[2][2];

    double detj   = qShape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );
    double weight = detj * qShape->weight[g];


    // -------------------------------------------------------------------------------------
    // compute values of quadratic shape functions at GP g

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = qShape->f[g];

    for( int j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // ------------------------------------------------------------------------------------
    // compute values of linear shape functions at GP g

    dfdxPtr = lShape->dfdx[g];
    dfdyPtr = lShape->dfdy[g];

    double dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    double* m = lShape->f[g];

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP i
    //             horizontal velocities   : U and V
    //             flow depth              : H
    //             turbulent kinetic energy: K
    //             dissipation of K        : D


    // integrate H, cf and vt with linear shape ------------------------------------------

    double H    = 0.0;
    double dHdx = 0.0;
    double dHdy = 0.0;
    double cf   = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE *node = elem->nd[j];

      double ndH  = node->v.S - node->z;
      if( ndH <= 0.0 )  ndH = 0.0;

      H    +=    m[j] * ndH;
      dHdx += dmdx[j] * ndH;
      dHdy += dmdy[j] * ndH;

      cf   +=    m[j] * node->cf;
    }

//    if( H <= project->hmin )  H = project->hmin;


    // integrate U, V, K, D and vt with quadratic shape ----------------------------------

    double U, dUdx, dUdy;
    double V, dVdx, dVdy;
    double K, dKdx, dKdy, dKdt, Ko;
    double D, dDdx, dDdy, dDdt, Do;

    U = dUdx = dUdy = 0.0;
    V = dVdx = dVdy = 0.0;
    K = dKdx = dKdy = dKdt = 0.0;
    D = dDdx = dDdy = dDdt = 0.0;

    for( int j=0; j<nnd; j++ )
    {
      NODE *node = elem->nd[j];

      double ndU = node->v.U;
      double ndV = node->v.V;

      U    +=    n[j] * ndU;
      dUdx += dndx[j] * ndU;
      dUdy += dndy[j] * ndU;

      V    +=    n[j] * ndV;
      dVdx += dndx[j] * ndV;
      dVdy += dndy[j] * ndV;

      double ndK = node->v.K - project->minK;
      double ndD = node->v.D - project->minD;

      K    +=    n[j] * ndK;
      dKdx += dndx[j] * ndK;
      dKdy += dndy[j] * ndK;

      D    +=    n[j] * ndD;
      dDdx += dndx[j] * ndD;
      dDdy += dndy[j] * ndD;

      dKdt +=    n[j] * node->v.dKdt;
      dDdt +=    n[j] * node->v.dDdt;
    }

    if( K < 0.0  ||  D < 0.0 )  K = D = 0.0;

    K += project->minK;
    D += project->minD;


    // -----------------------------------------------------------------------------------

    double vtxx = 0.0;
    double vtxy = 0.0;
    double vtyy = 0.0;

    for( int j=0; j<nnd; j++ )
    {
      NODE* node = elem->nd[j];

      vtxx += n[j] * node->exx * node->vt;
      vtxy += n[j] * node->exy * node->vt;
      vtyy += n[j] * node->eyy * node->vt;
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vtxx < project->minVtKD )  vtxx = project->minVtKD;
      if( vtyy < project->minVtKD )  vtyy = project->minVtKD;
    }


    // -----------------------------------------------------------------------------------
    // compute production and dissipation terms

    double Fi    = 2.0*dUdx*dUdx + 2.0*dVdy*dVdy + (dUdy+dVdx)*(dUdy+dVdx);

    double Us    = sqrt( U*U + V*V );
    double Utau  = sqrt( cf ) * Us;

    double cK    = 1.0 / sqrt( cf );
    double PKv   = 0.0;
    if( H > 0.0 )  PKv = cK * Utau*Utau*Utau / H;

    double cD    = c2D * sqrt( cmd/est/st / sqrt(cf*cf*cf) );
    double PDv   = 0.0;
    if( H > 0.0 )  PDv = cD * (Utau*Utau/H) * (Utau*Utau/H);

    double PK    =  cmd * K * K * Fi;
    double dPKdK =  cmd * 2.0 * K * Fi;
//  double dPKdD =  0.0;

//  double PD    =  c1D * cmd * K * K * Fi;
    double PD    =  c1D * PK;
//  double dPDdK =  c1D * cmd * 2.0 * K * Fi;
    double dPDdK =  c1D * dPKdK;
//  double dPDdD =  0.0;

    double DK    =  D * D;
//  double dDKdK =  0.0;
    double dDKdD =  2.0 * D;

//  double DD    =  c2D * D * D;
    double DD    =  c2D * DK;
//  double dDDdK =  0.0;
//  double dDDdD =  c2D * 2.0 * D;
    double dDDdD =  c2D * dDKdD;


    // -----------------------------------------------------------------------------------
    // compute diffusion terms

    double vtsKxx = vtxx / sK  +  vk;
    double vtsKxy = vtxy / sK  +  vk;
    double vtsKyy = vtyy / sK  +  vk;

    double vtsDxx = vtxx / sD  +  vk;
    double vtsDxy = vtxy / sD  +  vk;
    double vtsDyy = vtyy / sD  +  vk;


    // -----------------------------------------------------------------------------------
    // compute K and D equation and coefficients of NEWTON-RAPHSON matrix

    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // ---------------------------------------------------------------------------------
      // compute K-equation

      f   = H * D * dKdt;                             // time
      f  += H * D * (U * dKdx  +  V * dKdy);          // advection

      fx  = H * D * (vtsKxx * dKdx + vtsKxy * dKdy);  // diffusion
      fy  = H * D * (vtsKxy * dKdx + vtsKyy * dKdy);

      f  += H * dDdx * (vtsKxx * dKdx + vtsKxy * dKdy);
      f  += H * dDdy * (vtsKxy * dKdx + vtsKyy * dKdy);

//      f  += dHdx * D * (vtsKxx * dKdx + vtsKxy * dKdy);
//      f  += dHdy * D * (vtsKxy * dKdx + vtsKyy * dKdy);

      f  += H * (DK - PK - D*PKv);                   // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f + dndx[j] * fx + dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute D-equation

      f   = H * K * dDdt;                             // time
      f  += H * K * (U * dDdx  +  V * dDdy);          // advection

      fx  = H * K * (vtsDxx * dDdx + vtsDxy * dDdy);  // diffusion
      fy  = H * K * (vtsDxy * dDdx + vtsDyy * dDdy);

      f  += H * dKdx * (vtsDxx * dDdx + vtsDxy * dDdy);
      f  += H * dKdy * (vtsDxy * dDdx + vtsDyy * dDdy);

//      f  += dHdx * K * (vtsDxx * dDdx + vtsDxy * dDdy);
//      f  += dHdy * K * (vtsDxy * dDdx + vtsDyy * dDdy);

        f  += H * (DD - PD - K*PDv);                    // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + startD;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f + dndx[j] * fx + dndy[j] * fy;
      }
    }


    if( estifm )
    {
      double  df__, dfx_, dfy_, df_x, df_y, dfxx, dfxy, dfyx, dfyy;
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double* estifmPtr;

      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      // K-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * D * relaxThdt_KD;
      df_x  =  weight * H * D * U;
      df_y  =  weight * H * D * V;

      dfxx  =  weight * H * D * vtsKxx;
      dfxy  =
      dfyx  =  weight * H * D * vtsKxy;
      dfyy  =  weight * H * D * vtsKyy;

      df_x +=  weight * H * dDdx * vtsKxx;
      df_y +=  weight * H * dDdx * vtsKxy;

      df_x +=  weight * H * dDdy * vtsKxy;
      df_y +=  weight * H * dDdy * vtsKyy;

//      df_x +=  weight * dHdx * D * vtsKxx;
//      df_y +=  weight * dHdx * D * vtsKxy;
//
//      df_x +=  weight * dHdy * D * vtsKxy;
//      df_y +=  weight * dHdy * D * vtsKyy;

      df__ +=  weight * H * (/*dDKdK*/ - dPKdK);

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] =                 dfxx * dndx[j]  +  dfxy * dndy[j];
        ty[j] =                 dfyx * dndx[j]  +  dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // D-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * dKdt;
      df__ +=  weight * H * (U*dKdx + V*dKdy);

      dfx_  =  weight * H * (vtsKxx*dKdx + vtsKxy*dKdy);
      dfy_  =  weight * H * (vtsKxy*dKdx + vtsKyy*dKdy);

      df_x  =  weight * H * (vtsKxx*dKdx + vtsKxy*dKdy);
      df_y  =  weight * H * (vtsKxy*dKdx + vtsKyy*dKdy);

//      df__ +=  weight * dHdx * (vtsKxx*dKdx + vtsKxy*dKdy);
//      df__ +=  weight * dHdy * (vtsKxy*dKdx + vtsKyy*dKdy);

      df__ +=  weight * H * (dDKdD /*- dPKdD*/ - PKv);

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j] + startD;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // K-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * dDdt;
      df__ +=  weight * H * (U*dDdx + V*dDdy);

      dfx_  =  weight * H * (vtsDxx*dDdx + vtsDxy*dDdy);
      dfy_  =  weight * H * (vtsDxy*dDdx + vtsDyy*dDdy);

      df_x  =  weight * H * (vtsDxx*dDdx + vtsDxy*dDdy);
      df_y  =  weight * H * (vtsDxy*dDdx + vtsDyy*dDdy);

//      df__ +=  weight * dHdx * (vtsDxx*dDdx + vtsDxy*dDdy);
//      df__ +=  weight * dHdy * (vtsDxy*dDdx + vtsDyy*dDdy);

      df__ +=  weight * H * (/*dDDdK*/ - dPDdK - PDv);

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j+ startD];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // D-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * K * relaxThdt_KD;
      df_x  =  weight * H * K * U;
      df_y  =  weight * H * K * V;

      dfxx  =  weight * H * K * vtsDxx;
      dfxy  =
      dfyx  =  weight * H * K * vtsDxy;
      dfyy  =  weight * H * K * vtsDyy;

      df_x +=  weight * H * dKdx * vtsDxx;
      df_y +=  weight * H * dKdx * vtsDxy;

      df_x +=  weight * H * dKdy * vtsDxy;
      df_y +=  weight * H * dKdy * vtsDyy;

//      df_x +=  weight * dHdx * K * vtsDxx;
//      df_y +=  weight * dHdx * K * vtsDxy;
//
//      df_x +=  weight * dHdy * K * vtsDxy;
//      df_y +=  weight * dHdy * K * vtsDyy;

      df__ +=  weight * H * (dDDdD /*- dPDdD*/);

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] =                 dfxx * dndx[j]  +  dfxy * dndy[j];
        ty[j] =                 dfyx * dndx[j]  +  dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + startD] + startD;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void EQS_KD2D::RegionL( ELEM*    elem,
                        PROJECT* project,
                        double** estifm,
                        double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  TYPE* type = TYPE::Getid( elem->type );

  double cmd = project->KD.cm * project->KD.cd;  // constants of k-epsilon model
  double sK  = project->KD.sK;
  double sD  = project->KD.sD;
  double c1D = project->KD.c1D;
  double c2D = project->KD.c2D;
  double est = type->KDest;
  double st  = type->st;

  SHAPE* lShape = elem->GetLShape();
//  SHAPE* qShape = elem->GetQShape();

  int ngp = lShape->ngp;                         // number of GAUSS points
  int ncn = lShape->nnd;                         // number of corner nodes
//  int nnd = qShape->nnd;                         // total number of nodes

  int startD = elem->GetQShape()->nnd;


  // kinematic viscosity -----------------------------------------------------------------

  double vk = project->vk;


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<ncn; i++ )
  {
    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve k-epsilon equations

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double trafo[2][2];

    double* dfdxPtr = lShape->dfdx[g];
    double* dfdyPtr = lShape->dfdy[g];

    double detj   = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
    double weight = detj * lShape->weight[g];


    // -----------------------------------------------------------------------------------
    // compute values of linear shape functions at GP g

    double  dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];
    double* m = lShape->f[g];

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j] + trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j] + trafo[1][1] * dfdyPtr[j];
    }


    // -------------------------------------------------------------------------------------
    // compute values of quadratic shape functions at GP g

//    double  dndx[kMaxNodes2D], dndy[kMaxNodes2D];
//    double* n = qShape->f[g];
//
//    dfdxPtr = qShape->dfdx[g];
//    dfdyPtr = qShape->dfdy[g];
//
//    for( int j=0; j<nnd; j++ )
//    {
//      dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
//      dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
//    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP g
    //             horizontal velocities   : U and V
    //             flow depth              : H
    //             turbulent kinetic energy: K
    //             dissipation of K        : D


    // integrate H, cf and vt with linear shape ------------------------------------------

    double cf, H;
    H  = 0.0;
    cf = 0.0;

    double K, dKdx, dKdy, dKdt;
    double D, dDdx, dDdy, dDdt;
    K = dKdx = dKdy = dKdt = 0.0;
    D = dDdx = dDdy = dDdt = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE* node = elem->nd[j];

      double ndH = node->v.S - node->z;
      if( ndH <= 0.0 )  ndH = 0.0;

      H  += m[j] * ndH;
      cf += m[j] * node->cf;

      double ndK = node->v.K;
      double ndD = node->v.D;

      K    +=    m[j] * ndK;
      dKdx += dmdx[j] * ndK;
      dKdy += dmdy[j] * ndK;

      D    +=    m[j] * ndD;
      dDdx += dmdx[j] * ndD;
      dDdy += dmdy[j] * ndD;

      dKdt +=    m[j] * node->v.dKdt;
      dDdt +=    m[j] * node->v.dDdt;
    }


    double U, dUdx, dUdy;
    double V, dVdx, dVdy;
    U = dUdx = dUdy = 0.0;
    V = dVdx = dVdy = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE* node = elem->nd[j];

      double ndU = node->v.U;
      double ndV = node->v.V;

      U    +=    m[j] * ndU;
      dUdx += dmdx[j] * ndU;
      dUdy += dmdy[j] * ndU;

      V    +=    m[j] * ndV;
      dVdx += dmdx[j] * ndV;
      dVdy += dmdy[j] * ndV;
    }

//    for( int j=0; j<nnd; j++ )
//    {
//      NODE* node = elem->nd[j];
//
//      double ndU = node->v.U;
//      double ndV = node->v.V;
//
//      U    +=    n[j] * ndU;
//      dUdx += dndx[j] * ndU;
//      dUdy += dndy[j] * ndU;
//
//      V    +=    n[j] * ndV;
//      dVdx += dndx[j] * ndV;
//      dVdy += dndy[j] * ndV;
//    }


    // -----------------------------------------------------------------------------------

    double vt = 0.0;

    for( int j=0; j<ncn; j++ )  vt += m[j] * elem->nd[j]->vt;

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < project->minVtKD )  vt = project->minVtKD;
    }


    // -----------------------------------------------------------------------------------
    // compute K and D equation and coefficients of NEWTON-RAPHSON matrix

    double Fi    = 2.0*dUdx*dUdx + 2.0*dVdy*dVdy + (dUdy+dVdx)*(dUdy+dVdx);

    double Us    = sqrt( U*U + V*V );
    double Utau  = sqrt( cf ) * Us;

    double cK    = 1.0 / sqrt( cf );
    double PKv   = 0.0;
    if( H > 0.0 )  PKv = cK * Utau*Utau*Utau / H;

    double cD    = c2D * sqrt( cmd/est/st / sqrt(cf*cf*cf) );
    double PDv   = 0.0;
    if( H > 0.0 )  PDv = cD * (Utau*Utau/H) * (Utau*Utau/H);

    double PK    =  cmd * K * K * Fi;
    double dPKdK =  cmd * 2.0 * K * Fi;

    double PD    =  c1D * PK;
    double dPDdK =  c1D * dPKdK;

    double DK    =  D * D;
    double dDKdD =  2.0 * D;

    double DD    =  c2D * DK;
    double dDDdD =  c2D * dDKdD;


    // compute diffusion terms -----------------------------------------------------------

    double vtsK = vt / sK  +  vk;
    double vtsD = vt / sD  +  vk;


    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // ---------------------------------------------------------------------------------
      // compute K-equation

      f   = H * D * dKdt;                        // time
      f  += H * D * (U * dKdx  +  V * dKdy);     // advection

      fx  = H * D * vtsK * dKdx;                 // diffusion
      fy  = H * D * vtsK * dKdy;

      f  += H * dDdx * vtsK * dKdx;
      f  += H * dDdy * vtsK * dKdy;

      f  += H * (DK - PK - D*PKv);               // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f  +  dmdx[j] * fx  +  dmdy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute D-equation

      f   = H * K * dDdt;                        // time
      f  += H * K * (U * dDdx  +  V * dDdy);     // advection

      fx  = H * K * vtsD * dDdx;                 // diffusion
      fy  = H * K * vtsD * dDdy;

      f  += H * dKdx * vtsD * dDdx;
      f  += H * dKdy * vtsD * dDdy;

      f  += H * (DD - PD - K*PDv);               // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + startD;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f  +  dmdx[j] * fx  +  dmdy[j] * fy;
      }
    }


    if( estifm )
    {
      double t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfx_, dfy_, dfxx, dfyy;

      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      // K-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * D * relaxThdt_KD;
      df_x  =  weight * H * D * U;
      df_y  =  weight * H * D * V;

      dfxx  =
      dfyy  =  weight * H * D * vtsK;

      df_x +=  weight * H * dDdx * vtsK;
      df_y +=  weight * H * dDdy * vtsK;

      df__ -=  weight * H * dPKdK;

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j];
        ty[j] =                                    dfyy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j];

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // D-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * dKdt;
      df__ +=  weight * H * (U*dKdx + V*dKdy);

      dfx_  =  weight * H * vtsK * dKdx;
      dfy_  =  weight * H * vtsK * dKdy;

      df_x  =  weight * H * vtsK * dKdx;
      df_y  =  weight * H * vtsK * dKdy;

      df__ +=  weight * H * (dDKdD - PKv);

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j] + startD;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // K-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * dDdt;
      df__ +=  weight * H * (U * dDdx  +  V * dDdy);

      dfx_  =  weight * H * vtsD * dDdx;
      dfy_  =  weight * H * vtsD * dDdy;

      df_x  =  weight * H * vtsD * dDdx;
      df_y  =  weight * H * vtsD * dDdy;

      df__ -=  weight * H * (dPDdK + PDv);

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j + startD];

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // D-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * K * relaxThdt_KD;
      df_x  =  weight * H * K * U;
      df_y  =  weight * H * K * V;

      dfxx  =
      dfyy  =  weight * H * K * vtsD;

      df_x +=  weight * H * dKdx * vtsD;
      df_y +=  weight * H * dKdy * vtsD;

      df__ +=  weight * H * dDDdD;

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j];
        ty[j] =                                    dfyy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j + startD] + startD;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void EQS_KD2D::RegionLAI( ELEM*    elem,
                          PROJECT* project,
                          double** estifm,
                          double*  force )
{
  int i, j, k;


  // -------------------------------------------------------------------------------------
  // initializations

  TYPE* type = TYPE::Getid( elem->type );

  double cmd = project->KD.cm * project->KD.cd;  // constants of k-epsilon model
  double sK  = project->KD.sK;
  double sD  = project->KD.sD;
  double c1D = project->KD.c1D;
  double c2D = project->KD.c2D;
  double est = type->KDest;
  double st  = type->st;

  SHAPE* lShape = elem->GetLShape();

  int ngp = lShape->ngp;                   // number of GAUSS points
  int ncn = lShape->nnd;                   // number of corner nodes

  int nnd    = elem->Getnnd();             // total number of nodes
  int startD = nnd;


  // kinematic viscosity -----------------------------------------------------------------

  double vk = project->vk;


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( i=1; i<nnd; i++ )
  {
    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve k-epsilon equations

  for( i=0; i<ngp; i++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double trafo[2][2];

    double* dfdxPtr = lShape->dfdx[i];
    double* dfdyPtr = lShape->dfdy[i];

    double detj   = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
    double weight = detj * lShape->weight[i];


    // -----------------------------------------------------------------------------------
    // compute values of linear shape functions at GP i

    double* m = lShape->f[i];
    double  dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    for( j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j] + trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j] + trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP i
    //             horizontal velocities   : U and V
    //             flow depth              : H
    //             turbulent kinetic energy: K
    //             dissipation of K        : D


    // integrate H, cf and vt with linear shape ------------------------------------------

    double cf, H;

    double U, dUdx, dUdy;
    double V, dVdx, dVdy;
    double K, dKdx, dKdy, dKdt;
    double D, dDdx, dDdy, dDdt;

    H  = 0.0;
    cf = 0.0;

    U = dUdx = dUdy = 0.0;
    V = dVdx = dVdy = 0.0;
    K = dKdx = dKdy = dKdt = 0.0;
    D = dDdx = dDdy = dDdt = 0.0;

    for( j=0; j<ncn; j++ )
    {
      NODE* node;
      double ndU, ndV, ndK, ndD;

      node = elem->nd[j];

      H  += m[j] * (elem->nd[j]->v.S - elem->nd[j]->z);
      cf += m[j] * elem->nd[j]->cf;

      ndU  = node->v.U;
      ndV  = node->v.V;
      ndK  = node->v.K;
      ndD  = node->v.D;

      U    +=    m[j] * ndU;
      dUdx += dmdx[j] * ndU;
      dUdy += dmdy[j] * ndU;

      V    +=    m[j] * ndV;
      dVdx += dmdx[j] * ndV;
      dVdy += dmdy[j] * ndV;

      K    +=    m[j] * ndK;
      dKdx += dmdx[j] * ndK;
      dKdy += dmdy[j] * ndK;

      D    +=    m[j] * ndD;
      dDdx += dmdx[j] * ndD;
      dDdy += dmdy[j] * ndD;

      dKdt +=    m[j] * node->v.dKdt;
      dDdt +=    m[j] * node->v.dDdt;
    }

    if( H <= 0.0 )  H = 0.0; //H = project->hmin;

    double vtpk = 0.0;

    double vtxx = 0.0;
    double vtxy = 0.0;
    double vtyy = 0.0;

    if(    isFS(project->actualTurb, BCONSET::kVtIterat)
        && isFS(project->actualTurb, BCONSET::kVtPrandtlKol) )
    {
      for( j=0; j<ncn; j++ )
      {
        NODE* node = elem->nd[j];

        double ndK = node->v.K;
        double ndD = node->v.D;

        if( fabs(ndD) > 1.0e-9 )
        {
          vtpk += m[j] * ndK*ndK/ndD;
          vtxx += m[j] * node->exx * ndK*ndK/ndD;
          vtxy += m[j] * node->exy * ndK*ndK/ndD;
          vtyy += m[j] * node->eyy * ndK*ndK/ndD;
        }
      }

      vtpk *= cmd;
      vtxx *= cmd;
      vtxy *= cmd;
      vtyy *= cmd;
    }
    else
    {
      for( j=0; j<ncn; j++ )
      {
        NODE* node = elem->nd[j];

        vtpk += m[j] * node->vt;
        vtxx += m[j] * node->exx * node->vt;
        vtxy += m[j] * node->exy * node->vt;
        vtyy += m[j] * node->eyy * node->vt;
      }
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vtxx < project->minVtKD )  vtxx = project->minVtKD;
      if( vtyy < project->minVtKD )  vtyy = project->minVtKD;
    }

    // -----------------------------------------------------------------------------------
    // compute production and dissipation terms

    double Fi    = 2.0*dUdx*dUdx + 2.0*dVdy*dVdy + (dUdy+dVdx)*(dUdy+dVdx);

    double Us    = sqrt( U*U + V*V );
    double Utau  = sqrt( cf ) * Us;

    double cK    = 1.0 / sqrt( cf );
    double PKv   = 0.0;
    if( H > 0.0 )  PKv = cK * Utau*Utau*Utau / H;

//  double cD    = c2D * sqrt(cmd) / sqrt(st*est) / pow(cf, 0.75);
    double cD    = c2D * sqrt( cmd/est/st / sqrt(cf*cf*cf) );
    double PDv   = 0.0;
    if( H > 0.0 )  PDv = cD * (Utau*Utau/H) * (Utau*Utau/H);

    double PK    =  vtpk * Fi;
    double dPKdK =  0.0;
    double dPKdD =  0.0;

    double PD    =  c1D * D * vtpk * Fi;
    double dPDdK =  0.0;
    double dPDdD =  c1D * vtpk * Fi;

    double DD    =  c2D * D * D;
    double dDDdK =  0.0;
    double dDDdD =  c2D * 2 * D;


    // compute diffusion terms -----------------------------------------------------------

    double vtsKxx = vtxx / sK  +  vk;
    double vtsKxy = vtxy / sK  +  vk;
    double vtsKyy = vtyy / sK  +  vk;

    double vtsDxx = vtxx / sD  +  vk;
    double vtsDxy = vtxy / sD  +  vk;
    double vtsDyy = vtyy / sD  +  vk;


    if( force )
    {
      double f, fx, fy;

      // ---------------------------------------------------------------------------------
      // compute K-equation

      f   = H * dKdt;                                 // time
      f  += H * (U * dKdx  +  V * dKdy);              // advection

      fx  = H * (vtsKxx * dKdx + vtsKxy * dKdy);      // diffusion
      fy  = H * (vtsKxy * dKdx + vtsKyy * dKdy);

      f  += H * (D - PK - PKv);                       // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      double* forcePtr = force;

      for( j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f + dmdx[j] * fx + dmdy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute D-equation

      f   = H * K * dDdt;                             // time
      f  += H * K * (U * dDdx  +  V * dDdy);          // advection

      fx  = H * K * (vtsDxx * dDdx + vtsDxy * dDdy);  // diffusion
      fy  = H * K * (vtsDxy * dDdx + vtsDyy * dDdy);

        f  += H * (DD - PD - K * PDv);                  // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + startD;

      for( j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f + dmdx[j] * fx + dmdy[j] * fy;
      }
    }


    if( estifm )
    {
      double t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double df__, df_x, df_y, dfx_, dfy_, dfxx, dfxy, dfyx, dfyy;

      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      // K-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * relaxThdt_KD;
      df_x  =  weight * H * U;
      df_y  =  weight * H * V;

      dfxx  =  weight * H * vtsKxx;
      dfxy  =  weight * H * vtsKxy;

      dfyx  =  weight * H * vtsKxy;
      dfyy  =  weight * H * vtsKyy;

      df__ -=  weight * H * dPKdK;

      for( j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j]  +  dfxy * dmdy[j];
        ty[j] =                 dfyx * dmdx[j]  +  dfyy * dmdy[j];
      }

      for( j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j];

        for( k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // D-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * (1.0 - dPKdD);

      for( j=0; j<ncn; j++ )
      {
        t[j] = df__ * m[j];
      }

      for( j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j] + startD;

        for( k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }


      // K-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * dDdt;
      df__ +=  weight * H * (U * dDdx  +  V * dDdy);

      dfx_  =  weight * H * (vtsDxx * dDdx + vtsDxy * dDdy);
      dfy_  =  weight * H * (vtsDxy * dDdx + vtsDyy * dDdy);

      df__ +=  weight * H * (dDDdK - dPDdK - PDv);

      for( j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j + startD];

        for( k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // D-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * K * relaxThdt_KD;
      df_x  =  weight * H * K * U;
      df_y  =  weight * H * K * V;

      dfxx  =  weight * H * K * vtsDxx;
      dfxy  =  weight * H * K * vtsDxy;

      dfyx  =  weight * H * K * vtsDxy;
      dfyy  =  weight * H * K * vtsDyy;

      df__ +=  weight * H * (dDDdD - dPDdD);

      for( j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j]  +  dfxy * dmdy[j];
        ty[j] =                 dfyx * dmdx[j]  +  dfyy * dmdy[j];
      }

      for( j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j + startD] + startD;

        for( k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void EQS_KD2D::RegionQ( ELEM*    elem,
                        PROJECT* project,
                        double** estifm,
                        double*  force,
                        int*     Keq,
                        int*     Deq )
{
  // -------------------------------------------------------------------------------------
  // initializations

  TYPE* type = TYPE::Getid( elem->type );

  double cmd = project->KD.cm * project->KD.cd;   // constants of k-epsilon model
  double sK  = project->KD.sK;
  double sD  = project->KD.sD;
  double c1D = project->KD.c1D;
  double c2D = project->KD.c2D;
  double est = type->KDest;
  double st  = type->st;

  SHAPE* lShape = elem->GetLShape();

  int ngp = lShape->ngp;                          // number of GAUSS points
  int ncn = lShape->nnd;                          // number of corner nodes


  // kinematic viscosity -----------------------------------------------------------------

  double vk = project->vk;


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<ncn; i++ )
  {
    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve k-epsilon equations

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double* dfdxPtr = lShape->dfdx[g];
    double* dfdyPtr = lShape->dfdy[g];

    double trafo[2][2];

    double detj   = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
    double weight = detj * lShape->weight[g];


    // ------------------------------------------------------------------------------------
    // values of linear shape functions at GP g

    double  dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];
    double* m = lShape->f[g];

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP i
    //             horizontal velocities   : U and V
    //             flow depth              : H
    //             turbulent kinetic energy: K
    //             dissipation of K        : D


    // integrate H and cf with linear shape ----------------------------------------------

    double H    = 0.0;
    double cf   = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE *node = elem->nd[j];

      double ndH  = node->v.S - node->z;
      if( ndH <= 0.0 )  ndH = 0.0;

      H  += m[j] * ndH;
      cf += m[j] * node->cf;
    }


    // integrate U, V, K, D and vt with linear shape -------------------------------------

    double U, dUdx, dUdy;
    double V, dVdx, dVdy;
    double K, dKdx, dKdy, dKdt;
    double D, dDdx, dDdy, dDdt;

    U = dUdx = dUdy = 0.0;
    V = dVdx = dVdy = 0.0;
    K = dKdx = dKdy = dKdt = 0.0;
    D = dDdx = dDdy = dDdt = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE *node = elem->nd[j];

      double ndU  = node->v.U;
      double ndV  = node->v.V;

      U    +=    m[j] * ndU;
      dUdx += dmdx[j] * ndU;
      dUdy += dmdy[j] * ndU;

      V    +=    m[j] * ndV;
      dVdx += dmdx[j] * ndV;
      dVdy += dmdy[j] * ndV;

      double ndK = node->v.K;
      double ndD = node->v.D;

      K    +=    m[j] * ndK;
      dKdx += dmdx[j] * ndK;
      dKdy += dmdy[j] * ndK;

      D    +=    m[j] * ndD;
      dDdx += dmdx[j] * ndD;
      dDdy += dmdy[j] * ndD;

      dKdt +=    m[j] * node->v.dKdt;
      dDdt +=    m[j] * node->v.dDdt;
    }


    // -----------------------------------------------------------------------------------

    double vt = 0.0;

    for( int j=0; j<ncn; j++ )  vt += m[j] * elem->nd[j]->vt;

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < project->minVtKD )  vt = project->minVtKD;
    }


    // -----------------------------------------------------------------------------------
    // compute production and dissipation terms

    double Fi    = 2.0*dUdx*dUdx + 2.0*dVdy*dVdy + (dUdy+dVdx)*(dUdy+dVdx);

    double Us    = sqrt( U*U + V*V );
    double Utau  = sqrt( cf ) * Us;

    double cK    = 1.0 / sqrt( cf );
    double PKv   = 0.0;
    if( H > 0.0 )  PKv = cK * Utau*Utau*Utau / H;

//  double cD    = c2D * sqrt(cmd) / sqrt(st*est) / pow(cf, 0.75);
    double cD    = c2D * sqrt( cmd/est/st / sqrt(cf*cf*cf) );
    double PDv   = 0.0;
    if( H > 0.0 )  PDv = cD * (Utau*Utau/H) * (Utau*Utau/H);

    double PK    =  cmd * K * K * Fi;
    double dPKdK =  cmd * 2.0 * K * Fi;

//  double PD    =  c1D * cmd * K * K * Fi;
    double PD    =  c1D * PK;
//  double dPDdK =  c1D * cmd * 2.0 * K * Fi;
    double dPDdK =  c1D * dPKdK;

    double DK    =  D * D;
    double dDKdD =  2.0 * D;

//  double DD    =  c2D * D * D;
    double DD    =  c2D * DK;
//  double dDDdD =  c2D * 2.0 * D;
    double dDDdD =  c2D * dDKdD;


    // -----------------------------------------------------------------------------------
    // compute diffusion terms

    double vtsK = vt / sK  +  vk;
    double vtsD = vt / sD  +  vk;


    // -----------------------------------------------------------------------------------
    // compute K and D equation and coefficients of NEWTON-RAPHSON matrix

    if( force )
    {
      double  f, fx, fy;

      // ---------------------------------------------------------------------------------
      // compute K-equation

      f   = H * D * dKdt;                        // time
      f  += H * D * (U * dKdx  +  V * dKdy);     // advection

      fx  = H * D * vtsK * dKdx;                 // diffusion
      fy  = H * D * vtsK * dKdy;

      f  += H * dDdx * vtsK * dKdx;
      f  += H * dDdy * vtsK * dKdy;

      f  += H * (DK - PK - D*PKv);               // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      for( int j=0; j<ncn; j++ )
      {
        int r = Keq[j];
        force[r] -= m[j] * f  +  dmdx[j] * fx  +  dmdy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute D-equation

      f   = H * K * dDdt;                        // time
      f  += H * K * (U * dDdx  +  V * dDdy);     // advection

      fx  = H * K * vtsD * dDdx;                 // diffusion
      fy  = H * K * vtsD * dDdy;

      f  += H * dKdx * vtsD * dDdx;
      f  += H * dKdy * vtsD * dDdy;

      f  += H * (DD - PD - K*PDv);               // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      for( int j=0; j<ncn; j++ )
      {
        int r = Deq[j];
        force[r] -= m[j] * f  +  dmdx[j] * fx  +  dmdy[j] * fy;
      }
    }


    if( estifm )
    {
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfx_, dfy_, dfxx, dfyy;

      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      // K-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * D * relaxThdt_KD;
      df_x  =  weight * H * D * U;
      df_y  =  weight * H * D * V;

      dfxx  =
      dfyy  =  weight * H * D * vtsK;

      df_x +=  weight * H * dDdx * vtsK;
      df_y +=  weight * H * dDdy * vtsK;

      df__ -=  weight * H * dPKdK;

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j];
        ty[j] =                                    dfyy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        int r = Keq[j];

        for( int k=0; k<ncn; k++ )
        {
          int c = Keq[k];
          estifm[r][c] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // D-derivative of K-equation ------------------------------------------------------

      df__  =  weight * H * dKdt;
      df__ +=  weight * H * (U*dKdx + V*dKdy);

      dfx_  =  weight * H * vtsK * dKdx;
      dfy_  =  weight * H * vtsK * dKdy;

      df_x  =  weight * H * vtsK * dKdx;
      df_y  =  weight * H * vtsK * dKdy;

      df__ +=  weight * H * (dDKdD - PKv);

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        int r = Keq[j];

        for( int k=0; k<ncn; k++ )
        {
          int c = Deq[k];
          estifm[r][c] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // K-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * dDdt;
      df__ +=  weight * H * (U * dDdx  +  V * dDdy);

      dfx_  =  weight * H * vtsD * dDdx;
      dfy_  =  weight * H * vtsD * dDdy;

      df_x  =  weight * H * vtsD * dDdx;
      df_y  =  weight * H * vtsD * dDdy;

      df__ -=  weight * H * (dPDdK + PDv);

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        int r = Deq[j];

        for( int k=0; k<ncn; k++ )
        {
          int c = Keq[k];
          estifm[r][c] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // D-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * K * relaxThdt_KD;
      df_x  =  weight * H * K * U;
      df_y  =  weight * H * K * V;

      dfxx  =
      dfyy  =  weight * H * K * vtsD;

      df_x +=  weight * H * dKdx * vtsD;
      df_y +=  weight * H * dKdy * vtsD;

      df__ +=  weight * H * dDDdD;

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j];
        ty[j] =                                    dfyy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        int r = Deq[j];

        for( int k=0; k<ncn; k++ )
        {
          int c = Deq[k];
          estifm[r][c] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void EQS_KD2D::RegionQAI( ELEM*    elem,
                          PROJECT* project,
                          double** estifm,
                          double*  force )
{
}
