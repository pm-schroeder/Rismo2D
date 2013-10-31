// ======================================================================================
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
// ======================================================================================
//
// -----------------------------------------------------------------------------------------------
// compute element stiffness matrix
// solve for momentum equations (U, V) and continuity equation
//
// methods:   EQS_PPE2D::coefs()
//            EQS_PPE2D::bound()    one-dimensional boundary elements
//            EQS_PPE2D::region()   two-dimensional region elements
//
// Michael Schroeder in September 1992
//                      November  1999
// -----------------------------------------------------------------------------------------------

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

  int nnd    = qShape->nnd;
  int startY = nnd;
  int startP = 2 * nnd;

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

  int r0X = 0;            int r1X = 1;                int r2X = 2;
  int r0Y = startY;       int r1Y = startY + 1;       int r2Y = startY + 2;


  // column index --------------------------------------------------------------------------------

  int c0P = startP;       int c1P = startP + 1;


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

  Rotate2D( nnd, elem->nd, 3, estifm, force );
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


  int startY = nnd;
  int startP = 2 * nnd;

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


      forcePtr = force + startP;

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

        estifmPtr = estifm[j] + startP;

        for( k=0; k<ncn; k++ )
        {
          estifmPtr[k] += dndx[j] * weight * m[k];
        }
      }


      // 2. equation, y-derivative of phi --------------------------------------------------------

      for( j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + startY] + startY;

        for( k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j] * weight * n[k];
        }

        estifmPtr = estifm[j + startY] + startP;

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
        estifmPtr = estifm[j + startP];

        for( k=0; k<nnd; k++ )
        {
          estifmPtr[k] += m[j] * (wh * dndx[k] + wdhdx * n[k]);
// #      estifmPtr[k] += weight * m[j] * dndx[k];
        }

        // -------------------------------------
        estifmPtr = estifm[j + startP] + startY;

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

  Rotate2D( nnd, nd, 3, estifm, force );
}
