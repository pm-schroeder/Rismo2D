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
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsKL2D.h"


int EQS_KL2D::Coefs( ELEM*    elem,
                     PROJECT* project,
                     double** estifm,
                     double*  force )
{
  if(    isFS(elem->flag, ELEM::kDry)
      || isFS(elem->flag, ELEM::kBound) ) return 0;

  Region( elem, project, estifm, force );

  return 1;
}


void EQS_KL2D::Region( ELEM*    elem,
                       PROJECT* project,
                       double** estifm,
                       double*  force )
{
  int i, j, k;


  // -------------------------------------------------------------------------------------
  // initializations
  // -------------------------------------------------------------------------------------

  TYPE* type = TYPE::Getid( elem->type );

  double cm  = project->KD.cm;
  double cd  = project->KD.cd;
  double sK  = project->KD.sK;

  SHAPE* lShape = elem->GetLShape();

  int ngp = lShape->ngp;                   // number of GAUSS points
  int ncn = lShape->nnd;                   // number of corner nodes

  int nnd    = elem->Getnnd();             // total number of nodes


  if( force )
  {
    for( i=0; i<maxEleq; i++ )  force[i] = 0.0;
  }

  if( estifm )
  {
    for( i=0; i<maxEleq; i++ )
    {
      double* estifmPtr = estifm[i];

      for( j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


  // kinematic viscosity -----------------------------------------------------------------

  double vk = project->vk;


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node
  // -------------------------------------------------------------------------------------

  double xmin, xmax, ymin, ymax;
  double x[kMaxNodes2D], y[kMaxNodes2D];

  xmin = xmax = x[0] = elem->nd[0]->x;
  ymin = ymax = y[0] = elem->nd[0]->y;

  for( i=1; i<nnd; i++ )
  {
    if( elem->nd[i]->x < xmin )  xmin = elem->nd[i]->x;
    if( elem->nd[i]->x > xmax )  xmax = elem->nd[i]->x;

    if( elem->nd[i]->y < ymin )  ymin = elem->nd[i]->y;
    if( elem->nd[i]->y > ymax )  ymax = elem->nd[i]->y;

    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;

  double dl = sqrt( (xmax-xmin)*(ymax-ymin) );
  double lm = type->lm;
  double ls = lm * dl;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve k-epsilon equations
  // -------------------------------------------------------------------------------------

  for( i=0; i<ngp; i++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions
    // -----------------------------------------------------------------------------------

    double trafo[2][2];

    double* dfdxPtr = lShape->dfdx[i];
    double* dfdyPtr = lShape->dfdy[i];

    double detj   = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
    double weight = detj * lShape->weight[i];


    // -----------------------------------------------------------------------------------
    // compute values of linear shape functions at GP
    // -----------------------------------------------------------------------------------

    double* m = lShape->f[i];
    double  dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    for( j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j] + trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j] + trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP
    //             horizontal velocities   : U and V
    //             flow depth              : H
    //             turbulent kinetic energy: K
    // -----------------------------------------------------------------------------------

    double cf, H;

    double U, dUdx, dUdy;
    double V, dVdx, dVdy;
    double K, dKdx, dKdy, dKdt;

    H  = 0.0;
    cf = 0.0;

    U = dUdx = dUdy = 0.0;
    V = dVdx = dVdy = 0.0;
    K = dKdx = dKdy = dKdt = 0.0;

    for( j=0; j<ncn; j++ )
    {
      NODE* node;
      double ndU, ndV, ndK;

      node = elem->nd[j];

      H  += m[j] * (elem->nd[j]->v.S - elem->nd[j]->z);
      cf += m[j] * elem->nd[j]->cf;

      ndU  = node->v.U;
      ndV  = node->v.V;
      ndK  = node->v.K;

      U    +=    m[j] * ndU;
      dUdx += dmdx[j] * ndU;
      dUdy += dmdy[j] * ndU;

      V    +=    m[j] * ndV;
      dVdx += dmdx[j] * ndV;
      dVdy += dmdy[j] * ndV;

      K    +=    m[j] * ndK;
      dKdx += dmdx[j] * ndK;
      dKdy += dmdy[j] * ndK;

      dKdt +=    m[j] * node->v.dKdt;
    }

    if( H <= 0.0 ) H = project->hmin;


    // -----------------------------------------------------------------------------------
    // compute K equation
    // -----------------------------------------------------------------------------------

    double cK   = 1.0 / sqrt( cf );

    double Ures = sqrt( U*U + V*V );
    double Ust  = sqrt( cf ) * Ures;

    double PKv  = cK * Ust*Ust*Ust / H;

    double Fi   = 2*dUdx*dUdx + 2*dVdy*dVdy + (dUdy+dVdx)*(dUdy+dVdx);


    // eddy viscosity vt and turbulent length ls -----------------------------------------

    double vt = 0.0;
    double ls = lm / sqrt( sqrt(cm*cm*cm/cd) );

    if( isFS(project->actualTurb, BCONSET::kVtIterat) )
    {
      for( j=0; j<ncn; j++ )
      {
        double ndK = elem->nd[j]->v.K;
        double ndD = elem->nd[j]->v.D;

        if( ndD > 0.0 ) vt += m[j] * ndK * ndK / ndD;
      }

      vt *= cm * cd;
    }
    else
    {
      for( j=0; j<ncn; j++ )  vt += m[j] * elem->nd[j]->vt;
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < project->minVtKD )  vt = project->minVtKD;
    }


    // production ------------------------------------------------------------------------
/*
    double PK    =  cm * K * K / D * Fi;
    double dPKdK =  0.0;// cm * 2 * K / D * Fi;
*/
    double PK    =  vt * Fi;
    double dPKdK =  0.0;


    // dissipation -----------------------------------------------------------------------

    double D = cd * sqrt( K * K * K ) / ls;


    // diffusion -------------------------------------------------------------------------

    double vtsK = vt / sK  +  vk;


    if( force )
    {
      double f, fx, fy;

      f   = dKdt;                          // time
      f  += U * dKdx  +  V * dKdy;         // advection

      fx  = dKdx * vtsK;                   // diffusion
      fy  = dKdy * vtsK;

      f  += D - PK - PKv;                  // production and dissipation

      f  *= weight;
      fx *= weight;
      fy *= weight;

      double* forcePtr = force;

      for( j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f + dmdx[j] * fx + dmdy[j] * fy;
      }
    }


    if( estifm )
    {
      double t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double df__, df_x, df_y, dfxx, dfyy;

      df__  =  weight * relaxThdt_KD;
      df_x  =  weight * U;
      df_y  =  weight * V;

      dfxx  =
      dfyy  =  weight * vtsK;

      df__ -=  weight * dPKdK;

      for( j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j];
        ty[j] =                                    dfyy * dmdy[j];
      }

      for( j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j];

        for( k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }
    }
  }
}
