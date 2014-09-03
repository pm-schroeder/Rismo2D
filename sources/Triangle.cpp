// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class SHAPE
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
// compute values of all shape functions at Gauss point g and its weight
// element type: TRIANGLE
// -------------------------------------------------------------------------------------------------

#include "Defs.h"
#include "Shape.h"


void SHAPE::shapeOfTriangle()
{
  double x, y;
  double xgp = 0.0;
  double ygp = 0.0;

  for( int g=0; g<ngp; g++ )
  {
    switch( degree )
    {
      case 3:
        xgp = loc_TRI3[g][0];
        ygp = loc_TRI3[g][1];
        weight[g] = weight_TRI3[g];
        break;

      case 5:
        xgp = loc_TRI5[g][0];
        ygp = loc_TRI5[g][1];
        weight[g] = weight_TRI5[g];
        break;

      case 7:
        xgp = loc_TRI7[g][0];
        ygp = loc_TRI7[g][1];
        weight[g] = weight_TRI7[g];
        break;
    }

    switch( nnd )
    {
      // linear shape function, corner nodes: 0,1,2
      case 3:
        for( int i=0; i<3; i++ )
        {
          f[g][i]    = lTriangle( i, xgp, ygp );
          dfdx[g][i] = lTriangleDx( i );
          dfdy[g][i] = lTriangleDy( i );
        }
        break;

      // hyperbolic shape function,  corner nodes: 0,1,2
      //                             center node:  3
      case 4:
        for( int i=0; i<4; i++ )
        {
          f[g][i]    = bTriangle( i, xgp, ygp );
          dfdx[g][i] = bTriangleDx( i, xgp, ygp );
          dfdy[g][i] = bTriangleDy( i, xgp, ygp );
        }
        break;

      // quadratic shape function,  corner nodes: 0,1,2
      //                           midside nodes: 3,4,5
      case 6:
        for( int i=0; i<6; i++ )
        {
          f[g][i]    = qTriangle( i, xgp, ygp );
          dfdx[g][i] = qTriangleDx( i, xgp, ygp );
          dfdy[g][i] = qTriangleDy( i, xgp, ygp );
        }
        break;
    }
  }


  for( int n=0; n<7; n++ )
  {
    localTriangle( n, &x, &y );

    switch( nnd )
    {
      // linear shape function, corner nodes: 0,1,2
      case 3:
        for( int i=0; i<3; i++ )
        {
          dndx[n][i] = lTriangleDx( i );
          dndy[n][i] = lTriangleDy( i );
        }
        break;

      // hyperbolic shape function,  corner nodes: 0,1,2
      //                             center node:  3
      case 4:
        for( int i=0; i<4; i++ )
        {
          dndx[n][i] = bTriangleDx( i, x, y );
          dndy[n][i] = bTriangleDy( i, x, y );
        }
        break;

      // quadratic shape function,  corner nodes: 0,1,2
      //                           midside nodes: 3,4,5
      case 6:
        for( int i=0; i<6; i++ )
        {
          dndx[n][i] = qTriangleDx( i, x, y );
          dndy[n][i] = qTriangleDy( i, x, y );
        }
        break;
    }
  }
}


// function to obtain local coordinates of nodes ---------------------------------------------------

int SHAPE::localTriangle( int node, double *xi, double *eta )
{
  switch( nnd )
  {
    case 3:
    case 6:
      switch( node )
      {
        case  0:                                      // corner nodes
          *xi = 0.0;      *eta = 0.0;       break;
        case  1:
          *xi = 1.0;      *eta = 0.0;       break;
        case  2:
          *xi = 0.0;      *eta = 1.0;       break;

        case  3:                                      // midside nodes
          *xi = 0.5;      *eta = 0.0;       break;
        case  4:
          *xi = 0.5;      *eta = 0.5;       break;
        case  5:
          *xi = 0.0;      *eta = 0.5;       break;

        default:
          return -1;
      }
      break;

    case 4:
      switch( node )
      {
        case  0:                                      // corner nodes
          *xi = 0.0;      *eta = 0.0;       break;
        case  1:
          *xi = 1.0;      *eta = 0.0;       break;
        case  2:
          *xi = 0.0;      *eta = 1.0;       break;
        case  3:                                      // center node
          *xi = 1.0/3.0;  *eta = 1.0/3.0;   break;

        default:
          return -1;
      }
      break;
  }

  return 0;
}


// function to compute values of a linear shape function (triangle) --------------------------------

double SHAPE::lTriangle( int node, double x, double y )
{
  switch( node )
  {
    case 0:
      return 1.0 - x - y;
    case 1:
      return x;
    case 2:
      return y;
  }

  return 0.0;
}


// function to compute dfdx of a linear shape function (triangle) --------------

double SHAPE::lTriangleDx( int node )
{
  switch( node )
  {
    case 0:
      return -1.0;
    case 1:
      return  1.0;
    case 2:
      return  0.0;
  }

  return 0.0;
}


// function to compute dfdy of a linear shape function (triangle) --------------

double SHAPE::lTriangleDy( int node )
{
  switch(node)
  {
    case 0:
      return -1.0;
    case 1:
      return  0.0;
    case 2:
      return  1.0;
  }

  return 0.0;
}


// function to compute values of a bubble shape function (triangle) -------------

double SHAPE::bTriangle( int node, double x, double y )
{
  switch( node )
  {
    case 0:  return (1.0 - x - y) * (1.0 - 9.0 * x * y);
    case 1:  return x * (1.0 - 9.0 * (1.0 - x - y) * y);
    case 2:  return y * (1.0 - 9.0 * (1.0 - x - y) * x);
    case 3:  return 27.0 * (1.0 - x - y) * x * y;
  }

  return 0.0;
}


// function to compute dfdx of a bubble shape function (triangle) --------------

double SHAPE::bTriangleDx( int node, double x, double y )
{
  switch( node )
  {
    case 0:  return -1.0 + 9.0 * y * (2.0 * x  +  y  -  1);
    case 1:  return  1.0 + 9.0 * y * (2.0 * x  +  y  -  1);
    case 2:  return        9.0 * y * (2.0 * x  +  y  -  1);
    case 3:  return      -27.0 * y * (2.0 * x  +  y  -  1);
  }

  return 0.0;
}


// function to compute dfdy of a bubble shape function (triangle) --------------

double SHAPE::bTriangleDy( int node, double x, double y )
{
  switch( node )
  {
    case 0:  return -1.0 + 9.0 * x * (2.0 * y  +  x  -  1);
    case 1:  return        9.0 * x * (2.0 * y  +  x  -  1);
    case 2:  return  1.0 + 9.0 * x * (2.0 * y  +  x  -  1);
    case 3:  return      -27.0 * x * (2.0 * y  +  x  -  1);
  }

  return 0.0;
}


// function to compute values of a quadratic shape function (triangle) ---------

double SHAPE::qTriangle( int node, double x, double y )
{
  int    i, j;
  double Li, Lj;

  switch( node )
  {
    case 0:
    case 1:
    case 2:
      Li = lTriangle (node, x, y);
      return (2.0*Li - 1.0) * Li;

    case 3:
    case 4:
    case 5:
      i = node % 3;
      j = (node + 1) % 3;
      Li = lTriangle (i, x, y);
      Lj = lTriangle (j, x, y);
      return 4.0 * Li * Lj;
  }

  return 0.0;
}


// function to compute dfdx of a quadratic shape function (triangle) -----------

double SHAPE::qTriangleDx( int node, double x, double y )
{
  int i, j;
  double Li, Lj, dLidx, dLjdx;

  switch( node )
  {
    case 0:
    case 1:
    case 2:
      Li = lTriangle (node, x, y);
      dLidx = lTriangleDx (node);
      return (4.0*Li - 1.0) * dLidx;

    case 3:
    case 4:
    case 5:
      i = node % 3;
      j = (node + 1) % 3;
      Li = lTriangle (i, x, y);
      Lj = lTriangle (j, x, y);
      dLidx = lTriangleDx (i);
      dLjdx = lTriangleDx (j);
      return 4.0 * (dLidx*Lj + Li*dLjdx);
  }

  return 0.0;
}


// function to compute dfdy of a quadratic shape function (triangle) -----------

double SHAPE::qTriangleDy( int node, double x, double y )
{
  int    i, j;
  double Li, Lj, dLidy, dLjdy;

  switch( node )
  {
    case 0:
    case 1:
    case 2:
      Li = lTriangle (node, x, y);
      dLidy = lTriangleDy (node);
      return (4.0*Li - 1.0) * dLidy;

    case 3:
    case 4:
    case 5:
      i = node % 3;
      j = (node + 1) % 3;
      Li = lTriangle (i, x, y);
      Lj = lTriangle (j, x, y);
      dLidy = lTriangleDy (i);
      dLjdy = lTriangleDy (j);
      return 4.0 * (dLidy*Lj + Li*dLjdy);
  }

  return 0.0;
}
