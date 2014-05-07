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
// compute diagonal entries of lumped mass matrix for surface and region
//
// Michael Schroeder in September 1994
//
// 24.09.1997   sc   2D linear elements
// 22.12.1999   sc
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Shape.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"

#include "Grid.h"


void GRID::LumpedMassMatrix( double** lmm )
{
  int i, j;


  // -------------------------------------------------------------------------------------
  // initializations

  *lmm = (double*) MEMORY::memo.Array_nd( np );

  for( i=0; i<np; i++ )  (*lmm)[i] = 0.0;


  // -------------------------------------------------------------------------------------
  // loop on all elements

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = &elem[e];

    SHAPE* lShape = el->GetLShape();

    int ngp = lShape->ngp;         // number of GAUSS points
    int ncn = lShape->nnd;         // number of corner nodes


    NODE** nd = el->nd;


    // -----------------------------------------------------------------------------------
    // compute coordinates relative to first node

    double x[kMaxNodes2D], y[kMaxNodes2D];

    x[0] = nd[0]->x;
    y[0] = nd[0]->y;

    for( i=1; i<ncn; i++ )
    {
      x[i] = nd[i]->x - x[0];
      y[i] = nd[i]->y - y[0];
    }

    x[0] = y[0] = 0.0;


    // -----------------------------------------------------------------------------------
    // use GAUSS point integration to solve for region elements

    for( i=0; i<ngp; i++ )
    {
      double* m = lShape->f[i];

      // ---------------------------------------------------------------------------------
      // transformation of integration area

      double* dfdx = lShape->dfdx[i];      // derivative of shape function
      double* dfdy = lShape->dfdy[i];      // in local co-ordinates

      double trafo[2][2];

      double detj   = lShape->jacobi2D( ncn, dfdx, dfdy, x, y, trafo );
      double weight = detj * lShape->weight[i];


      // ---------------------------------------------------------------------------------
      // compute integral

      for( j=0; j<ncn; j++ )
      {
        (*lmm)[nd[j]->Getno()] += m[j] * weight;
      }
    }
  }
}
