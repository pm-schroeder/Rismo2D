// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_K2D
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
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsK2D.h"


int EQS_K2D::Coefs( ELEM*    elem,
                    PROJECT* project,
                    double** estifm,
                    double*  force )
{
  if( isFS(elem->flag, ELEM::kDry)  ||  isFS(elem->flag, ELEM::kBound) )  return 0;

  Region( elem, project, estifm, force );

  return 1;
}


void EQS_K2D::Region( ELEM*    elem,
                      PROJECT* project,
                      double** estifm,
                      double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  TYPE* type = TYPE::Getid( elem->type );

  double cmd = project->KD.cm * project->KD.cd;   // constants of k-epsilon model
  double sK  = project->KD.sK;
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
  // use GAUSS point integration to solve k-equation

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
    double K, dKdx, dKdy, dKdt;
    double D, dDdx, dDdy;

    U = dUdx = dUdy = 0.0;
    V = dVdx = dVdy = 0.0;
    K = dKdx = dKdy = dKdt = 0.0;
    D = dDdx = dDdy = 0.0;

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

      dKdt +=    n[j] * node->v.dKdt;
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

    double cK    = 1.0 / sqrt( cf );
    double PKv   = 0.0;
    if( H > 0.0 )  PKv = cK * Utau*Utau*Utau / H;

    double PK    =  cmd * K * K * Fi;
    double dPKdK =  cmd * 2.0 * K * Fi;


    // -----------------------------------------------------------------------------------
    // compute diffusion terms

    double vtsK = vt / sK  +  vk;


    // -----------------------------------------------------------------------------------
    // compute K equation and coefficients of NEWTON-RAPHSON matrix

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

      f  += H * (D*D - PK - D*PKv);              // production and dissipation

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
        estifmPtr = estifm[j];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }
    }
  }
}
