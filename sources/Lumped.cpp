// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class GRID
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

// -------------------------------------------------------------------------------------------------
// compute diagonal entries of lumped mass matrix for surface and region
// -------------------------------------------------------------------------------------------------

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
