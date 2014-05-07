// ======================================================================================
//                                      G R I D
// ======================================================================================
// This method GRID::VeloGrad() computes all velocity gradients.
// ======================================================================================
//
// Copyright (C) 1992-2007  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version)
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 05.04.2007     sc     first implementation / first concept
//
// ======================================================================================
//
#include "Defs.h"
#include "Report.h"
#include "Memory.h"
#include "Shape.h"
#include "Node.h"
#include "Elem.h"
#include "Project.h"
#include "Model.h"

#include "Grid.h"


void GRID::VeloGrad( PROJECT* project, double* dUdx, double* dUdy,
                                       double* dVdx, double* dVdy )
{
  GRID* rg   = project->M2D->region;    // pointer to the grid
  int   rgnp = rg->Getnp();             // number of nodes in the grid
  int   rgne = rg->Getne();             // number of elements in the grid


  // allocate memory for temporary use ---------------------------------------------------

  double* A = (double*) MEMORY::memo.Array_nd( rgnp );


  // loop over elements ------------------------------------------------------------------

  for( int e=0; e<rgne; e++ )
  {
    ELEM* elem = rg->Getelem(e);

    int ncn = elem->Getncn();
    int nnd = elem->Getnnd();

    SHAPE* shape = elem->GetQShape();

    // local element coordinates ---------------------------------------------------------
    double x[kMaxNodes2D], y[kMaxNodes2D];

    x[0] = elem->nd[0]->x;
    y[0] = elem->nd[0]->y;

    for( int i=1; i<nnd; i++ )
    {
      x[i] = elem->nd[i]->x - x[0];
      y[i] = elem->nd[i]->y - y[0];
    }
    x[0] = y[0] = 0.0;

    // integrate over element area: A[e] -------------------------------------------------
    for( int g=0; g<shape->ngp; g++ )
    {
      double* dfdxPtr = shape->dfdx[g];
      double* dfdyPtr = shape->dfdy[g];

      double  trafo[2][2];

      double  detj   = shape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
      double  weight = detj * shape->weight[g];

      // compute values of quadratic shape functions at GP g
      double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

      for( int j=0; j<nnd; j++ )
      {
        dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
        dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
      }

      double Ux = 0.0;
      double Uy = 0.0;

      double Vx = 0.0;
      double Vy = 0.0;

      for( int j=0; j<nnd; j++ )
      {
        NODE* node = elem->nd[j];

        double ndU = node->v.U;
        double ndV = node->v.V;

        Ux += dndx[j] * ndU;
        Uy += dndy[j] * ndU;

        Vx += dndx[j] * ndV;
        Vy += dndy[j] * ndV;
      }

      for( int j=0; j<nnd; j++ )
      {
        NODE* node = elem->nd[j];
        int   no   = node->Getno();

        dUdx[no] += weight * Ux;
        dUdy[no] += weight * Uy;
        dVdx[no] += weight * Vx;
        dVdy[no] += weight * Vy;

        A[no]    += weight;
      }
    }
  }


  // compute gradients from integrated values; divide by area A[] ------------------------

  for( int n=0; n<rgnp; n++ )
  {
    NODE* node = rg->Getnode(n);
    int   no   = node->Getno();

    dUdx[no] /= A[no];
    dUdy[no] /= A[no];
    dVdx[no] /= A[no];
    dVdy[no] /= A[no];
  }


  // detach temporarily used memory ------------------------------------------------------

  MEMORY::memo.Detach( A );
}
