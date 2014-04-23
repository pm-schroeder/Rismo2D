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
// ---------------------------------------------------------------------------------------
// Compute rotation of flow field:  rot = dVdx - dUdy
// The finite elements collocation method is used.
//
// Michael Schroeder in July 2006
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Shape.h"
#include "Memory.h"
#include "Type.h"
#include "Vars.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"


double* MODEL::Rot2D()
{
  int np = region->Getnp();
  int ne = region->Getne();


  // -------------------------------------------------------------------------------------
  // allocate memory for rotation of flow field
  // -------------------------------------------------------------------------------------

  double* rot = (double*) MEMORY::memo.Array_nd( np );
  double* wgt = (double*) MEMORY::memo.Array_nd( np );

  for( int n=0; n<np; n++ )  rot[n] = wgt[n] = 0.0;


  // -------------------------------------------------------------------------------------
  // loop on elements
  // -------------------------------------------------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* elem = region->Getelem(e);

    SHAPE* shape = elem->GetQShape();

    int ngp = shape->ngp;          // number of GAUSS points
    int nnd = shape->nnd;          // number of corner nodes


    // -----------------------------------------------------------------------------------
    // compute coordinates relative to first node

    double x[kMaxNodes2D], y[kMaxNodes2D];

    x[0] = elem->nd[0]->x;
    y[0] = elem->nd[0]->y;

    for( int i=1; i<nnd; i++ )
    {
      x[i] = elem->nd[i]->x - *x;
      y[i] = elem->nd[i]->y - *y;
    }
    x[0] = y[0] = 0.0;


    // -----------------------------------------------------------------------------------
    // use GAUSS point integration to solve momentum and continuity

    for( int g=0; g<ngp; g++ )
    {
      // form JACOBIAN transformation matrix ---------------------------------------------

      double  trafo[2][2];

      double* dfdxPtr = shape->dfdx[g];
      double* dfdyPtr = shape->dfdy[g];

      double detj   = shape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );
      double weight = detj * shape->weight[g];

      // compute values of shape functions at GP g ---------------------------------------

      double* n = shape->f[g];

      double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

      for( int i=0; i<nnd; i++ )
      {
        dndx[i] = trafo[0][0] * dfdxPtr[i] + trafo[0][1] * dfdyPtr[i];
        dndy[i] = trafo[1][0] * dfdxPtr[i] + trafo[1][1] * dfdyPtr[i];
      }

      // compute flow parameters and their derivatives -----------------------------------

      double dUdy = 0.0;
      double dVdx = 0.0;

      for( int i=0; i<nnd; i++ )
      {
        NODE* node = elem->nd[i];

        dUdy += dndy[i] * node->v.U;
        dVdx += dndx[i] * node->v.V;
      }

      // compute rotation ----------------------------------------------------------------

      double f = weight * ( dVdx - dUdy );

      for( int i=0; i<nnd; i++ )
      {
        NODE* nd = elem->nd[i];
        int   no = nd->Getno();

        rot[no] += f;
        wgt[no] += weight;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // solve for rotation

  for( int n=0; n<np; n++ )
  {
    rot[n] /= wgt[n];
  }

  // -------------------------------------------------------------------------------------

  MEMORY::memo.Detach( wgt );

  return rot;
}
