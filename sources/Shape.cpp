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

// ---------------------------------------------------------------------------------------
// Determine values of shape functions at Gauss points.
//
// element types:
//     one-dimensional  elements     : LINE
//     two-dimensional  elements     : TRIANGLE, SQUARE
//
// shapes                            : linear and quadratic
//
// Gauss point integration           : 3rd, 5th or 7th order
//
// ---------------------------------------------------------------------------------------
// input parameters :  GPdeg         : degree of Gauss point integration
//                                     (possible values: 3, 5 or 7)
// ---------------------------------------------------------------------------------------
//
// Michael Schroeder in June    1992
//                      March   1993
//                      October 1999
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Memory.h"
#include "Report.h"

#include "Shape.h"


int    SHAPE::nShape = 0;
SHAPE* SHAPE::pShape = NULL;

double SHAPE::loc_1D[4][4];
double SHAPE::weight_1D[4][4];

double SHAPE::loc_TRI3[4][2];
double SHAPE::weight_TRI3[4];
double SHAPE::loc_TRI5[7][2];
double SHAPE::weight_TRI5[7];
double SHAPE::loc_TRI7[13][2];
double SHAPE::weight_TRI7[13];


SHAPE::SHAPE()
{
}


SHAPE::~SHAPE()
{
  if( nShape )  delete[]  pShape;
}


void SHAPE::init( int GPdeg )
{
  gp();               // initialize GAUSS point locations and weights
  types( GPdeg );     // initialize shape types


  // initialize values of shape functions at Gauss points --------------------------------

  for( int i=0; i<nShape; i++ )
  {
    SHAPE* s = &pShape[i];

    // compute values at Gauss points ----------------------------------------------------

    switch( s->ident )
    {
      // one-dimensional element
      case kLine:
        s->shapeOfLine();
        break;

      // two-dimensional elements
      case kTriangle:
        s->shapeOfTriangle();
        break;

      case kSquare:
        s->shapeOfSquare();
        break;
    }
  }
}


SHAPE* SHAPE::get( int ident, int order )
{
  for( int i=0; i<nShape; i++ )
  {
    SHAPE* s = &pShape[i];

    if( s->ident == ident  &&  s->order == order )  return s;
  }

  REPORT::rpt.Error( kParameterFault, "shape with ID=%d and ORDER=%d is not supported %s",
                                      ident, order, "(SHAPE::get)" );
  return NULL;
}


void SHAPE::types( int GPdeg )
{
  // allocate memory for shapes ----------------------------------------------------------

  nShape = 7;
  pShape = new SHAPE [nShape];

  if( !pShape )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (SHAPE::types - 1)" );


  // initialization ----------------------------------------------------------------------

  for( int i=0; i<nShape; i++ ) pShape[i].degree = GPdeg;


  SHAPE* s;

  switch( GPdeg )
  {
    int i;

    case 3:
      i=0; s = &pShape[i]; s->dim = 1; s->ident = kLine; s->ngp =  2; s->order = 1; s->nnd = 2;
      i++; s = &pShape[i]; s->dim = 1; s->ident = kLine; s->ngp =  2; s->order = 2; s->nnd = 3;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp =  4; s->order = 1; s->nnd = 3;
      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp =  4; s->order = 2; s->nnd = 6;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp =  4; s->order = 3; s->nnd = 4;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kQuad; s->ngp =  4; s->order = 1; s->nnd = 4;
      i++; s = &pShape[i]; s->dim = 2; s->ident = kQuad; s->ngp =  4; s->order = 2; s->nnd = 8;
      break;

    case 5:
      i=0; s = &pShape[i]; s->dim = 1; s->ident = kLine; s->ngp =  3; s->order = 1; s->nnd = 2;
      i++; s = &pShape[i]; s->dim = 1; s->ident = kLine; s->ngp =  3; s->order = 2; s->nnd = 3;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp =  7; s->order = 1; s->nnd = 3;
      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp =  7; s->order = 2; s->nnd = 6;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp =  7; s->order = 3; s->nnd = 4;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kQuad; s->ngp =  9; s->order = 1; s->nnd = 4;
      i++; s = &pShape[i]; s->dim = 2; s->ident = kQuad; s->ngp =  9; s->order = 2; s->nnd = 8;
      break;

    case 7:
      i=0; s = &pShape[i]; s->dim = 1; s->ident = kLine; s->ngp =  4; s->order = 1; s->nnd = 2;
      i++; s = &pShape[i]; s->dim = 1; s->ident = kLine; s->ngp =  4; s->order = 2; s->nnd = 3;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp = 13; s->order = 1; s->nnd = 3;
      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp = 13; s->order = 2; s->nnd = 6;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kTri;  s->ngp = 13; s->order = 3; s->nnd = 4;

      i++; s = &pShape[i]; s->dim = 2; s->ident = kQuad; s->ngp = 16; s->order = 1; s->nnd = 4;
      i++; s = &pShape[i]; s->dim = 2; s->ident = kQuad; s->ngp = 16; s->order = 2; s->nnd = 8;
      break;
  }
}


void SHAPE::gp()
{
  // location of 1D-Gauss points (1st-, 3rd-, 5th- and 7th-order) ------------------------

  loc_1D[0][0] =  0.0;
  loc_1D[0][1] =  0.0;
  loc_1D[0][2] =  0.0;
  loc_1D[0][3] =  0.0;
  loc_1D[1][0] = -1.0 / sqrt(3.0);
  loc_1D[1][1] =  1.0 / sqrt(3.0);
  loc_1D[1][2] =  0.0;
  loc_1D[1][3] =  0.0;
  loc_1D[2][0] = -sqrt(0.6);
  loc_1D[2][1] =  0.0;
  loc_1D[2][2] =  sqrt(0.6);
  loc_1D[2][3] =  0.0;
  loc_1D[3][0] = -0.861136311594053;
  loc_1D[3][1] = -0.339981043584856;
  loc_1D[3][2] =  0.339981043584856;
  loc_1D[3][3] =  0.861136311594053;


  // weight of Gauss points --------------------------------------------------------------

  weight_1D[0][0] =  2.0;
  weight_1D[0][1] =  0.0;
  weight_1D[0][2] =  0.0;
  weight_1D[0][3] =  0.0;
  weight_1D[1][0] =  1.0;
  weight_1D[1][1] =  1.0;
  weight_1D[1][2] =  0.0;
  weight_1D[1][3] =  0.0;
  weight_1D[2][0] =  5.0/9.0;
  weight_1D[2][1] =  8.0/9.0;
  weight_1D[2][2] =  5.0/9.0;
  weight_1D[2][3] =  0.0;
  weight_1D[3][0] =  0.347854845137454;
  weight_1D[3][1] =  0.652145154862546;
  weight_1D[3][2] =  0.652145154862546;
  weight_1D[3][3] =  0.347854845137454;


  // location of triangular Gauss points (3rd-order)--------------------------------------

  loc_TRI3[0][0] = 1.0/3.0;
  loc_TRI3[0][1] = 1.0/3.0;
  loc_TRI3[1][0] = 0.6;
  loc_TRI3[1][1] = 0.2;
  loc_TRI3[2][0] = 0.2;
  loc_TRI3[2][1] = 0.6;
  loc_TRI3[3][0] = 0.2;
  loc_TRI3[3][1] = 0.2;


  // weight of Gauss points --------------------------------------------------------------

  weight_TRI3[0] = -0.28125;
  weight_TRI3[1] =  0.260416666666667;
  weight_TRI3[2] =  0.260416666666667;
  weight_TRI3[3] =  0.260416666666667;


  // location of triangular Gauss points (5th-order) -------------------------------------

  loc_TRI5[0][0] = 1.0/3.0;
  loc_TRI5[0][1] = 1.0/3.0;
  loc_TRI5[1][0] = 0.0597158717897698;
  loc_TRI5[1][1] = 0.470142064105115;
  loc_TRI5[2][0] = 0.470142064105115;
  loc_TRI5[2][1] = 0.0597158717897698;
  loc_TRI5[3][0] = 0.470142064105115;
  loc_TRI5[3][1] = 0.470142064105115;
  loc_TRI5[4][0] = 0.797426985353087;
  loc_TRI5[4][1] = 0.101286507323456;
  loc_TRI5[5][0] = 0.101286507323456;
  loc_TRI5[5][1] = 0.797426985353087;
  loc_TRI5[6][0] = 0.101286507323456;
  loc_TRI5[6][1] = 0.101286507323456;


  // weight of Gauss points --------------------------------------------------------------

  weight_TRI5[0] = 0.1125;
  weight_TRI5[1] = 0.0661970763942531;
  weight_TRI5[2] = 0.0661970763942531;
  weight_TRI5[3] = 0.0661970763942531;
  weight_TRI5[4] = 0.0629695902724136;
  weight_TRI5[5] = 0.0629695902724136;
  weight_TRI5[6] = 0.0629695902724136;


  // location of triangular Gauss points (7th-order) -------------------------------------

  loc_TRI7[ 0][0] = 1.0/3.0;
  loc_TRI7[ 0][1] = 1.0/3.0;
  loc_TRI7[ 1][0] = 0.479308067841923;
  loc_TRI7[ 1][1] = 0.260345966079038;
  loc_TRI7[ 2][0] = 0.260345966079038;
  loc_TRI7[ 2][1] = 0.479308067841923;
  loc_TRI7[ 3][0] = 0.260345966079038;
  loc_TRI7[ 3][1] = 0.260345966079038;
  loc_TRI7[ 4][0] = 0.869739794195568;
  loc_TRI7[ 4][1] = 0.065130102902216;
  loc_TRI7[ 5][0] = 0.065130102902216;
  loc_TRI7[ 5][1] = 0.869739794195568;
  loc_TRI7[ 6][0] = 0.065130102902216;
  loc_TRI7[ 6][1] = 0.065130102902216;
  loc_TRI7[ 7][0] = 0.638444188569809;
  loc_TRI7[ 7][1] = 0.312865496004875;
  loc_TRI7[ 8][0] = 0.312865496004875;
  loc_TRI7[ 8][1] = 0.638444188569809;
  loc_TRI7[ 9][0] = 0.638444188569809;
  loc_TRI7[ 9][1] = 0.048690315425316;
  loc_TRI7[10][0] = 0.048690315425316;
  loc_TRI7[10][1] = 0.638444188569809;
  loc_TRI7[11][0] = 0.312865496004875;
  loc_TRI7[11][1] = 0.048690315425316;
  loc_TRI7[12][0] = 0.048690315425316;
  loc_TRI7[12][1] = 0.312865496004875;


  // weight of Gauss points --------------------------------------------------------------

  weight_TRI7[ 0] = -0.074785022233835;
  weight_TRI7[ 1] =  0.087807628716602;
  weight_TRI7[ 2] =  0.087807628716602;
  weight_TRI7[ 3] =  0.087807628716602;
  weight_TRI7[ 4] =  0.026673617804420;
  weight_TRI7[ 5] =  0.026673617804420;
  weight_TRI7[ 6] =  0.026673617804420;
  weight_TRI7[ 7] =  0.038556880445129;
  weight_TRI7[ 8] =  0.038556880445129;
  weight_TRI7[ 9] =  0.038556880445129;
  weight_TRI7[10] =  0.038556880445129;
  weight_TRI7[11] =  0.038556880445129;
  weight_TRI7[12] =  0.038556880445129;
}


// ---------------------------------------------------------------------------------------
// function to compute components of the jacobian transformation matrix
// for twodimensional elements; return value is the jacobian's determinant

double SHAPE::jacobi2D( int     nnd,
                        double* dfdx,
                        double* dfdy,
                        double* x,
                        double* y,
                        double  trafo[2][2] )
{
  double det;
  double jacobi[2][2];


  // initialize jacobian -----------------------------------------------------------------

  jacobi[0][0] = 0.0;
  jacobi[0][1] = 0.0;
  jacobi[1][0] = 0.0;
  jacobi[1][1] = 0.0;


  // compute components of jacobian ------------------------------------------------------

  for( int n=0; n<nnd; n++ )
  {
    jacobi[0][0] += dfdx[n] * x[n];
    jacobi[0][1] += dfdx[n] * y[n];

    jacobi[1][0] += dfdy[n] * x[n];
    jacobi[1][1] += dfdy[n] * y[n];
  }


  // compute the determinant -------------------------------------------------------------

  det = jacobi[0][0] * jacobi[1][1] - jacobi[0][1] * jacobi[1][0];


  // compute the inverse matrix ----------------------------------------------------------

  trafo[0][0] =  jacobi[1][1] / det;
  trafo[0][1] = -jacobi[0][1] / det;
  trafo[1][0] = -jacobi[1][0] / det;
  trafo[1][1] =  jacobi[0][0] / det;

  return (det);
}


void SHAPE::getCornerNodes( int mn, int* cn1, int* cn2 )
{
  *cn1 = *cn2 = -1;

  switch( ident )
  {
    case kLine:
      switch( mn )
      {
        case  2: *cn1 =  0;  *cn2 =  1;  break;
      }
      break;

    case kTriangle:
      switch( mn )
      {
        case  3: *cn1 =  0;  *cn2 =  1;  break;
        case  4: *cn1 =  1;  *cn2 =  2;  break;
        case  5: *cn1 =  2;  *cn2 =  0;  break;
      }
      break;

    case kSquare:
      switch( mn )
      {
        case  4: *cn1 =  0;  *cn2 =  1;  break;
        case  5: *cn1 =  1;  *cn2 =  2;  break;
        case  6: *cn1 =  2;  *cn2 =  3;  break;
        case  7: *cn1 =  3;  *cn2 =  0;  break;
      }
      break;
  }
}
