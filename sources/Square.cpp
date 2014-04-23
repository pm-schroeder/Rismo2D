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

/* ------------------------------------------------------------------------ *
 * compute values of all shape functions at Gauss point g and its weight    *
 * element type: SQUARE                                                     *
 * ------------------------------------------------------------------------ */

#include "Defs.h"
#include "Shape.h"


void SHAPE::shapeOfSquare()
{
  int dg = (this->degree + 1) / 2;

  for( int g=0; g<this->ngp; g++ )
  {
    int igp_x = g % dg;
    int igp_y = g / dg;

    double xgp = this->loc_1D[dg-1][igp_x];
    double ygp = this->loc_1D[dg-1][igp_y];

    this->weight[g] = this->weight_1D[dg-1][igp_x] *
                      this->weight_1D[dg-1][igp_y];

    switch( this->nnd )
    {
      // linear shape function, corner nodes: 0,1,2,3
       case 4:
        for( int i=0; i<4; i++ )
        {
          this->f[g][i]    = lSquare(   i, xgp, ygp );
          this->dfdx[g][i] = lSquareDx( i,      ygp );
          this->dfdy[g][i] = lSquareDy( i, xgp      );
        }
        break;

      // quadratic shape function,  corner nodes: 0,1,2,3
      //                           midside nodes: 4,5,6,7
      case 8:
        for( int i=0; i<8; i++ )
        {
          this->f[g][i]    = qSquare(   i, xgp, ygp );
          this->dfdx[g][i] = qSquareDx( i, xgp, ygp );
          this->dfdy[g][i] = qSquareDy( i, xgp, ygp );
        }
        break;
    }
  }


  for( int n=0; n<8; n++ )
  {
    double x, y;
    localSquare( n, &x, &y );

    switch( this->nnd )
    {
      // linear shape function, corner nodes: 0,1,2,3
       case 4:
        for( int i=0; i<4; i++ )
        {
          this->dndx[n][i] = lSquareDx( i,    y );
          this->dndy[n][i] = lSquareDy( i, x    );
        }
        break;

      // quadratic shape function,  corner nodes: 0,1,2,3
      //                           midside nodes: 4,5,6,7
      case 8:
        for( int i=0; i<8; i++ )
        {
          this->dndx[n][i] = qSquareDx( i, x, y );
          this->dndy[n][i] = qSquareDy( i, x, y );
        }
        break;
    }
  }
}

/* ------------------------------------------------------------------------ *
 * compute local coordinates of node                                        */

int SHAPE::localSquare( int node, double *xi, double *eta )
{
  switch(node)
  {
  case  0:
    *xi   = -1.0; *eta  = -1.0; break;
  case  1:
    *xi   =  1.0; *eta  = -1.0; break;
  case  2:
    *xi   =  1.0; *eta  =  1.0; break;
  case  3:
    *xi   = -1.0; *eta  =  1.0; break;

  case  4:
    *xi   =  0.0; *eta  = -1.0; break;
  case  5:
    *xi   =  1.0; *eta  =  0.0; break;
  case  6:
    *xi   =  0.0; *eta  =  1.0; break;
  case  7:
    *xi   = -1.0; *eta  =  0.0; break;

  default:
    return (-1);
  }
  return (0);
}

/* function to compute values of a linear shape function (square) */
double SHAPE::lSquare( int node, double x, double y )
{
  double xi, eta;       /* local coordinates of node */

  localSquare(node, &xi, &eta);

  return ( (1.0 + xi*x) * (1.0 + eta*y) / 4.0 );
}

/* function to compute dfdx of a quadratic shape function (square) */
double SHAPE::lSquareDx( int node, double y )
{
  double xi, eta;       /* local coordinates of node */

  localSquare(node, &xi, &eta);

  return ( xi * (1.0 + eta*y) / 4.0 );
}

/* function to compute dfdy of a quadratic shape function (square) */
double SHAPE::lSquareDy( int node, double x )
{
  double xi, eta;       /* local coordinates of node */

  localSquare(node, &xi, &eta);

  return( eta * (1.0 + xi*x) / 4.0 );
}

/* function to compute values of a quadratic shape function (square) */
double SHAPE::qSquare( int node, double x, double y )
{
  double xi, eta;       /* local coordinates of node */

  localSquare(node, &xi, &eta);

  switch(node)
  {
  case 0:
  case 1:
  case 2:
  case 3:
    return( (1.0 + xi*x) * (1.0 + eta*y) * (xi*x + eta*y -1.0) / 4.0 );

  case 4:
  case 6:
    return( (1.0 - x*x) * (1.0 + eta*y) / 2.0 );

  case 5:
  case 7:
    return( (1.0 + xi*x) * (1.0-y*y) / 2.0 );
  }

  return 0.0;
}

/* function to compute dfdx of a quadratic shape function (square) */
double SHAPE::qSquareDx( int node, double x, double y )
{
  double xi, eta;       /* local coordinates of node */

  localSquare(node, &xi, &eta);

  switch(node)
  {
  case 0:
  case 1:
  case 2:
  case 3:
    return( (1.0 + eta*y) * (2.0*x + xi*eta*y) / 4.0 );

  case 4:
  case 6:
    return( -x * (1.0 + eta*y) );

  case 5:
  case 7:
    return( xi * (1.0 - y*y) / 2.0 );
  }

  return 0.0;
}

/* function to compute dfdy of a quadratic shape function (square) */
double SHAPE::qSquareDy( int node, double x, double y )
{
  double xi, eta;       /* local coordinates of node */

  localSquare(node, &xi, &eta);

  switch(node)
  {
  case 0:
  case 1:
  case 2:
  case 3:
    return( (1.0 + xi*x) * (2.0*y + xi*eta*x) / 4.0 );

  case 4:
  case 6:
    return( eta * (1.0 - x*x) / 2.0 );

  case 5:
  case 7:
    return( -y * (1.0 + xi*x) );
  }

  return 0.0;
}
