// =======================================================================================
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
// =======================================================================================

#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsD2D.h"


int EQS_D2D::Coefs( ELEM*    elem,
                    PROJECT* project,
                    double** estifm,
                    double*  force )
{
  if( isFS(elem->flag, ELEM::kDry)  ||  isFS(elem->flag, ELEM::kBound) )  return 0;

  Region( elem, project, estifm, force );

  return 1;
}


void EQS_D2D::Region( ELEM*    elem,
                      PROJECT* project,
                      double** estifm,
                      double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  TYPE* type = TYPE::Getid( elem->type );

  double cmd = project->KD.cm * project->KD.cd;   // constants of k-epsilon model
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

  if( force )
  {
    for( int i=0; i<maxEleq; i++ )  force[i] = 0.0;
  }

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      double* estifmPtr = estifm[i];

      for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


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

    double* m = lShape->f[g];


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

//  if( H <= 0.0 )  H = 0.0; //H = project->hmin;


    // integrate U, V, K, D and vt with quadratic shape ----------------------------------

    double U, dUdx, dUdy;
    double V, dVdx, dVdy;
    double K, dKdx, dKdy;
    double D, dDdx, dDdy, dDdt;

    U = dUdx = dUdy = 0.0;
    V = dVdx = dVdy = 0.0;
    K = dKdx = dKdy = 0.0;
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

      double ndK = node->v.K;
      double ndD = node->v.D;

      K    +=    n[j] * ndK;
      dKdx += dndx[j] * ndK;
      dKdy += dndy[j] * ndK;

      D    +=    n[j] * ndD;
      dDdx += dndx[j] * ndD;
      dDdy += dndy[j] * ndD;

      dDdt +=    n[j] * node->v.dDdt;
    }

    if( K < 0.0 )  K = project->minK;
    if( D < 0.0 )  D = project->minD;


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

//  double cD    = c2D * sqrt(cmd) / sqrt(st*est) / pow(cf, 0.75);
    double cD    = c2D * sqrt( cmd/est/st / sqrt(cf*cf*cf) );
    double PDv   = 0.0;
    if( H > 0.0 )  PDv = cD * (Utau*Utau/H) * (Utau*Utau/H);

    double PD    =  c1D * cmd * K * K * Fi;


    // -----------------------------------------------------------------------------------
    // compute diffusion terms

    double vtsD = vt / sD  +  vk;


    // -----------------------------------------------------------------------------------
    // compute D equation and coefficients of NEWTON-RAPHSON matrix

    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // ---------------------------------------------------------------------------------
      // compute D-equation

      f   = H * K * dDdt;                        // time
      f  += H * K * (U * dDdx  +  V * dDdy);     // advection

      fx  = H * K * vtsD * dDdx;                 // diffusion
      fy  = H * K * vtsD * dDdy;

      f  += H * dKdx * vtsD * dDdx;
      f  += H * dKdy * vtsD * dDdy;

      f  += H * (c2D*D*D - PD - K*PDv);          // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f + dndx[j] * fx + dndy[j] * fy;
      }
    }


    if( estifm )
    {
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfxx, dfyy;
      double* estifmPtr;

      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      // D-derivative of D-equation ------------------------------------------------------

      df__  =  weight * H * K * relaxThdt_KD;
      df_x  =  weight * H * K * U;
      df_y  =  weight * H * K * V;

      dfxx  =
      dfyy  =  weight * H * K * vtsD;

      df_x +=  weight * H * dKdx * vtsD;
      df_y +=  weight * H * dKdy * vtsD;

      df__ +=  weight * H * c2D * 2.0 * D;

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] =                 dfxx * dndx[j];
        ty[j] =                                    dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }
    }
  }
}
