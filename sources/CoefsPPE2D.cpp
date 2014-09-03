// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_PPE2D
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

// #define DEBUG 1

#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsPPE2D.h"


int EQS_PPE2D::Coefs( ELEM*    elem,
                      PROJECT* project,
                      double** estifm,
                      double*  force )
{
  if( isFS(elem->flag, ELEM::kDry) ) return 0;


  if( isFS(elem->flag, ELEM::kBound) )  Bound( elem, project, estifm, force );
  else                                  Region( elem, project, estifm, force );

  return 1;
}


void EQS_PPE2D::Bound( ELEM*    elem,
                       PROJECT* project,
                       double** estifm,
                       double*  force )
{
  int     i, j;
  double* estifmPtr;

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int nnd   = qShape->nnd;
  int eqidY = nnd;
  int eqidP = 2 * nnd;

  if( force ) for( i=0; i<maxEleq; i++ ) force[i] = 0.0;

  if( estifm )
  {
    for( i=0; i<maxEleq; i++ )
    {
      estifmPtr = estifm[i];

      for( j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


  if( !estifm ) return;

  if( isFS(elem->flag, ELEM::kOutlet) )  return;


  NODE   *node[3];

  node[0] = elem->nd[0];      /* corner nodes */
  node[1] = elem->nd[1];
  node[2] = elem->nd[2];      /* midside node */


  // row index -----------------------------------------------------------------------------------

  int r0X = 0;           int r1X = 1;               int r2X = 2;
  int r0Y = eqidY;       int r1Y = eqidY + 1;       int r2Y = eqidY + 2;


  // column index --------------------------------------------------------------------------------

  int c0P = eqidP;       int c1P = eqidP + 1;


  for( j=0; j<qShape->ngp; j++ )       // loop on GAUSS points
  {
    double weight, nx, ny;

    double* m  = lShape->f[j];         // linear shape
    double* n  = qShape->f[j];         // quadratic shape
    double* dn = qShape->dfdx[j];


    // -------------------------------------------------------------------------------------------
    // compute normal vector
    // notice: the normal is not reduced to unit length

    nx =  dn[0]*node[0]->y + dn[1]*node[1]->y + dn[2]*node[2]->y;
    ny = -dn[0]*node[0]->x - dn[1]*node[1]->x - dn[2]*node[2]->x;


    // weight of Gauss point j -------------------------------------------------------------------

    weight = qShape->weight[j];


    // compute estifm ----------------------------------------------------------------------------

    estifm[r0X][c0P] -= n[0] * weight * nx * m[0];
    estifm[r0X][c1P] -= n[0] * weight * nx * m[1];

    estifm[r1X][c0P] -= n[1] * weight * nx * m[0];
    estifm[r1X][c1P] -= n[1] * weight * nx * m[1];

    estifm[r2X][c0P] -= n[2] * weight * nx * m[0];
    estifm[r2X][c1P] -= n[2] * weight * nx * m[1];

    estifm[r0Y][c0P] -= n[0] * weight * ny * m[0];
    estifm[r0Y][c1P] -= n[0] * weight * ny * m[1];

    estifm[r1Y][c0P] -= n[1] * weight * ny * m[0];
    estifm[r1Y][c1P] -= n[1] * weight * ny * m[1];

    estifm[r2Y][c0P] -= n[2] * weight * ny * m[0];
    estifm[r2Y][c1P] -= n[2] * weight * ny * m[1];
  }


  // apply transformation ------------------------------------------------------------------------

  Rotate2D( nnd, eqidY, elem->nd, estifm, force );
}


void EQS_PPE2D::Region( ELEM*    elem,
                        PROJECT* project,
                        double** estifm,
                        double*  force )
{
  int     i, j, k;
  double* forcePtr;
  double* estifmPtr;
  double f;
  double weight, detj;
  double wh, wdhdx, wdhdy;
  double h, dhdt, dhdx, dhdy;
  double u, dudx;
  double v, dvdy;
  double trafo[2][2];
  double x[kMaxNodes2D], y[kMaxNodes2D];
  double dmdx[kMaxNodes2D], dmdy[kMaxNodes2D],
         dndx[kMaxNodes2D], dndy[kMaxNodes2D];


  // ---------------------------------------------------------------------------------------------
  // initializations

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int ngp = qShape->ngp;         // number of GAUSS points
  int ncn = lShape->nnd;         // number of corner nodes
  int nnd = qShape->nnd;         // number of all nodes


  if( force ) for( i=0; i<maxEleq; i++ ) force[i] = 0.0;

  if( estifm )
  {
    for( i=0; i<maxEleq; i++ )
    {
      estifmPtr = estifm[i];

      for( j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


  int eqidY = nnd;
  int eqidP = 2 * nnd;

  NODE** nd = elem->nd;


  // ---------------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  x[0] = nd[0]->x;
  y[0] = nd[0]->y;

  for( i=1; i<nnd; i++ )
  {
    x[i] = nd[i]->x - *x;
    y[i] = nd[i]->y - *y;
  }

  x[0] = y[0] = 0.0;


  // ---------------------------------------------------------------------------------------------
  // GAUSS point integration loop

  for( i=0; i<ngp; i++ )
  {
    // -------------------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with shape functions

    double* dfdxPtr = qShape->dfdx[i];
    double* dfdyPtr = qShape->dfdy[i];

    detj = qShape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );

    weight = detj * qShape->weight[i];


    // -------------------------------------------------------------------------------------------
    // compute values of quadaratic shape functions at GP i

    double* n = qShape->f[i];

    for( j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // -------------------------------------------------------------------------------------------
    // compute values of linear shape functions at GP i

    dfdxPtr = lShape->dfdx[i];
    dfdyPtr = lShape->dfdy[i];

    double* m = lShape->f[i];

    for( j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // -------------------------------------------------------------------------------------------
    // compute linear parameters at GP i

    h = dhdt = dhdx = dhdy = 0.0;

    for( j=0; j<ncn; j++ )
    {
      double ndH;

      ndH  = nd[j]->v.S     - nd[j]->z;

      h    +=    m[j] * ndH;
      dhdx += dmdx[j] * ndH;
      dhdy += dmdy[j] * ndH;

      dhdt +=    m[j] * nd[j]->v.dSdt;
    }


    // -------------------------------------------------------------------------------------------
    // compute right side

    if( force )
    {
      // -----------------------------------------------------------------------------------------
      // compute qadratic parameters at GP i

      u = dudx = v = dvdy = 0.0;

// #  double duhdx = 0.0;
// #  double dvhdy = 0.0;

      for( j=0; j<nnd; j++ )
      {
// # ------------------------------------------------------------------------
// #    the following version, the minimisation of U*H and V*H instead of
// #    U and V, does not lead to acceptable results !!
// #
// #    double H = nd[j]->v.S - nd[j]->z;

// #    duhdx += dndx[j] * (H * nd[j]->v.U);
// #    dvhdy += dndy[j] * (H * nd[j]->v.V);

        u    +=    n[j] * nd[j]->v.U;
        dudx += dndx[j] * nd[j]->v.U;

        v    +=    n[j] * nd[j]->v.V;
        dvdy += dndy[j] * nd[j]->v.V;
      }


      // 3. equation -----------------------------------------------------------------------------

      f = weight * ( dhdt + h*(dudx+dvdy) + u*dhdx + v*dhdy );
// #  f = weight * ( dhdt + duhdx + dvhdy );


      forcePtr = force + eqidP;

      for( j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f;
      }
    }


    // -------------------------------------------------------------------------------------------
    // compute coefficients of element stiffness matrix

    if( estifm )
    {
      // 1. equation, x-derivative of phi --------------------------------------------------------

      for( j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j];

        for( k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j] * weight * n[k];
        }

        estifmPtr = estifm[j] + eqidP;

        for( k=0; k<ncn; k++ )
        {
          estifmPtr[k] += dndx[j] * weight * m[k];
        }
      }


      // 2. equation, y-derivative of phi --------------------------------------------------------

      for( j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + eqidY] + eqidY;

        for( k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j] * weight * n[k];
        }

        estifmPtr = estifm[j + eqidY] + eqidP;

        for( k=0; k<ncn; k++ )
        {
          estifmPtr[k] += dndy[j] * weight * m[k];
        }
      }


      // 3. equation, phi ------------------------------------------------------------------------

      wdhdx = weight * dhdx;
      wdhdy = weight * dhdy;
      wh    = weight * h;

      for( j=0; j<ncn; j++ )
      {
        // ----------------------------
        estifmPtr = estifm[j + eqidP];

        for( k=0; k<nnd; k++ )
        {
          estifmPtr[k] += m[j] * (wh * dndx[k] + wdhdx * n[k]);
// #      estifmPtr[k] += weight * m[j] * dndx[k];
        }

        // -------------------------------------
        estifmPtr = estifm[j + eqidP] + eqidY;

        for( k=0; k<nnd; k++ )
        {
          estifmPtr[k] += m[j] * (wh * dndy[k] + wdhdy * n[k]);
// #      estifmPtr[k] += weight * m[j] * dndy[k];
        }

        // -------------------------------------
        //estifmPtr = estifm[j + startP] + startP;

        //for( k=0; k<ncn; k++ )
        //{
        //  estifmPtr[k] += m[j] * weight * diagWeightPPE;
        //}
      }
    }
  }


  // apply transformation ------------------------------------------------------------------------

  Rotate2D( nnd, eqidY, nd, estifm, force );
}
