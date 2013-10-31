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
// compute values of all shape functions at Gauss point g and its weight
// element type: LINE
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Shape.h"


void SHAPE::shapeOfLine()
{
  int    n, g;
  double x, xgp;

  n = this->ngp - 1;

  for( g=0; g<this->ngp; g++ )
  {
    xgp = this->loc_1D[n][g];

    this->weight[g] = this->weight_1D[n][g];

    switch( this->nnd )
    {
      case 2:      // linear shape function, corner nodes: 0,1
        this->f[g][0] = lLine( 0, xgp );
        this->f[g][1] = lLine( 1, xgp );

        this->dfdx[g][0] = lLineDx( 0 );
        this->dfdx[g][1] = lLineDx( 1 );

        break;

      case 3:      // quadratic shape function, corner nodes: 0,1
                   //                           midside node: 2
        this->f[g][0] = qLine( 0, xgp );
        this->f[g][1] = qLine( 1, xgp );
        this->f[g][2] = qLine( 2, xgp );

        this->dfdx[g][0] = qLineDx( 0, xgp );
        this->dfdx[g][1] = qLineDx( 1, xgp );
        this->dfdx[g][2] = qLineDx( 2, xgp );
        break;
    }
  }


  for( n=0; n<this->nnd; n++ )
  {
    localLine( n, &x );

    switch( this->nnd )
    {
      case 2:      // linear shape function, corner nodes: 0,1
        this->dndx[n][0] = lLineDx( 0 );
        this->dndx[n][1] = lLineDx( 1 );

        break;

      case 3:      // quadratic shape function, corner nodes: 0,1
                   //                           midside node: 2
        this->dndx[n][0] = qLineDx( 0, x );
        this->dndx[n][1] = qLineDx( 1, x );
        this->dndx[n][2] = qLineDx( 2, x );
        break;
    }
  }
}


void SHAPE::localLine( int node, double* xi )
{
  switch( node )
  {
    case  0:   *xi = -1.0;   break;
    case  1:   *xi =  1.0;   break;
    case  2:   *xi =  0.0;   break;
  }
}


double SHAPE::lLine( int node, double x )
{
       if( node == 0 )  return (1.0 - x) / 2.0;
  else if( node == 1 )  return (1.0 + x) / 2.0;
  else                  return (0.0);
}


double SHAPE::lLineDx( int node )
{
       if( node == 0 )  return -0.5;
  else if( node == 1 )  return  0.5;
  else                  return  0.0;
}


double SHAPE::qLine( int node, double x )
{
  switch( node )
  {
    case 0:   return x * (x - 1.0) / 2.0;
    case 1:   return x * (x + 1.0) / 2.0;
    case 2:   return 1.0 - x * x;
    default:  return 0.0;
  }
}


double SHAPE::qLineDx( int node, double x )
{
  switch( node )
  {
    case 0:   return x - 0.5;
    case 1:   return x + 0.5;
    case 2:   return -2.0 * x;
    default:  return 0.0;
  }
}
