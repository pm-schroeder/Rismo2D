// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_Disp
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

#include "EqsDisp.h"


int EQS_DISP::Coefs( ELEM*    elem,
                     PROJECT* project,
                     double** estifm,
                     double*  force )
{
  if(    isFS(elem->flag, ELEM::kDry)
      || isFS(elem->flag, ELEM::kBound) ) return 0;

  Region( elem, project, estifm, force );

  return 1;
}


void EQS_DISP::Region( ELEM*    elem,
                       PROJECT* project,
                       double** estifm,
                       double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations
  // -------------------------------------------------------------------------------------

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int ncn = lShape->nnd;            // number of corner nodes
  int nnd = qShape->nnd;            // number of nodes
  int ngp = qShape->ngp;            // number of GAUSS points


  TYPE* type = TYPE::Getid( elem->type );

  double betaSf   = type->betaSf;
  double mueSf    = project->mueSf;
  double maxTanSf = project->maxTanSf;
  double minUSf   = project->minUSf;
  double kappa    = project->kappa;


  if( force )
  {
    for( int i=0; i<maxEleq; i++ ) force[i] = 0.0;
  }

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      double* estifmPtr = estifm[i];

      for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node
  // -------------------------------------------------------------------------------------

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    x[i] = elem->nd[i]->x - *x;
    y[i] = elem->nd[i]->y - *y;
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve momentum and continuity
  // -------------------------------------------------------------------------------------

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with quadratic shape functions
    // -----------------------------------------------------------------------------------

    double  trafo[2][2];
    double* dfdxPtr = lShape->dfdx[g];
    double* dfdyPtr = lShape->dfdy[g];

    double detj = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );

    double weight = detj * lShape->weight[g];


    // -----------------------------------------------------------------------------------
    // compute values of shape functions at GP i
    // -----------------------------------------------------------------------------------

    double* m = lShape->f[g];

    double dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j] + trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j] + trafo[1][1] * dfdyPtr[j];
    }

    // -----------------------------------------------------------------------------------

    dfdxPtr = qShape->dfdx[g];
    dfdyPtr = qShape->dfdy[g];

    double* n = qShape->f[g];

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    for( int j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j] + trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j] + trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives
    // -----------------------------------------------------------------------------------

    double H  = 0.0;
    double cf = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      NODE* node = elem->nd[i];

      H  += m[i] * ( node->v.S - node->z );
      cf += m[i] * ( node->cf );
    }

    if( H < 0.0 )  H = 0.0;


    // -----------------------------------------------------------------------------------

    double U    = 0.0;
    double dUdx = 0.0;
    double dUdy = 0.0;

    double V    = 0.0;
    double dVdx = 0.0;
    double dVdy = 0.0;


    for( int i=0; i<nnd; i++ )
    {
      NODE* node = elem->nd[i];

      double ndU = node->v.U;
      U    +=    n[i] * ndU;
      dUdx += dndx[i] * ndU;
      dUdy += dndy[i] * ndU;

      double ndV = node->v.V;
      V    +=    n[i] * ndV;
      dVdx += dndx[i] * ndV;
      dVdy += dndy[i] * ndV;
    }


    double detx = U * dVdx  -  V * dUdx;
    double dety = U * dVdy  -  V * dUdy;

    double Us   = sqrt( U*U + V*V );

    double uvdet = U * detx  +  V * dety;
    double rho   = 0.0;

    // limit minimum velocity ------------------------------------------------------------
    if( Us > minUSf )  rho = uvdet / Us / Us / Us;

    // limit maximum curvature -----------------------------------------------------------
    //if( rho > 0.01 )  rho = 0.01;

    double sx = 0.0;
    double sy = 0.0;

    if( Us > 1.0e-6 )
    {
      sx = U / Us;
      sy = V / Us;
    }

    double mm   = kappa / sqrt( cf );
    double alfa = betaSf * (mm + 0.5) / kappa / kappa / mm;

    // limit maximum secondary velocity --------------------------------------------------
    double tanSf = alfa * H * rho;
    if( fabs(tanSf) > maxTanSf  ) tanSf *= maxTanSf / fabs(tanSf);


    // -----------------------------------------------------------------------------------
    // right hand side (RHS)
    // -----------------------------------------------------------------------------------

    if( force )
    {
      double f = weight * tanSf * Us;

      //for( int j=0; j<nnd; j++ )  force[j] += n[j] * f;
      for( int j=0; j<ncn; j++ )  force[j] += m[j] * f;
    }


    // -----------------------------------------------------------------------------------
    // coefficients (LHS)
    // -----------------------------------------------------------------------------------

    if( estifm )
    {
      double lambda = 0.0;
//      if( Us > 1.0e-6 )  lambda = mueSf * mm * H / 2.0 / kappa / kappa / Us; ???
      if( Us > 1.0e-6 )  lambda = mueSf * mm * H / kappa;

      double df__ = weight;
      double df_x = weight * lambda * sx;
      double df_y = weight * lambda * sy;

      //for( int j=0; j<nnd; j++ )
      //{
      //  double* estifmPtr = estifm[j];

      //  for( int k=0; k<nnd; k++ )
      //  {
      //    estifmPtr[k] += n[j] * ( df__ * n[k]  +  df_x * dndx[k]  +  df_y * dndy[k] );
      //  }
      //}
      for( int j=0; j<ncn; j++ )
      {
        double* estifmPtr = estifm[j];

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j] * ( df__ * m[k]  +  df_x * dmdx[k]  +  df_y * dmdy[k] );
        }
      }
    }
  }
}
