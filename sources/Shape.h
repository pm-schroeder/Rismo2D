// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// S H A P E
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Shape.h      : definition file of the class.
// Shape.cpp    : implementation file of the class.
//
// Line.cpp     : methods SHAPE::shapeOfLine()
//                        SHAPE::localLine()
//                        SHAPE::lLine()
//                        SHAPE::lLineDx()
//                        SHAPE::qLine()
//                        SHAPE::qLineDx()
// Square.cpp   : methods SHAPE::shapeOfSquare()
//                        SHAPE::localSquare()
//                        SHAPE::lSquare()
//                        SHAPE::lSquareDx()
//                        SHAPE::lSquareDy()
//                        SHAPE::qSquare()
//                        SHAPE::qSquareDx()
//                        SHAPE::qSquareDy()
// Triangle.cpp : methods SHAPE::shapeOfTriangle()
//                        SHAPE::localTriangle()
//                        SHAPE::lTriangle()
//                        SHAPE::lTriangleDx()
//                        SHAPE::lTriangleDy()
//                        SHAPE::bTriangle()
//                        SHAPE::bTriangleDx()
//                        SHAPE::bTriangleDy()
//                        SHAPE::qTriangle()
//                        SHAPE::qTriangleDx()
//                        SHAPE::qTriangleDy()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the shape functions on Finite Elements.
//
// -------------------------------------------------------------------------------------------------
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
// -------------------------------------------------------------------------------------------------
//
// HISTORY
//
//    date              changes
// ------------  ----  -----------------------------------------------------------------------------
//  01.01.1992    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHAPE_INCL
#define SHAPE_INCL 1

#include "Report.h"


#define kMaxGP2D    16         // maximum number of Gauss points, 2D
#define kMaxND       8         // maximum number of nodes

#define kLine        0         // element types
#define kTri         1
#define kTriangle    1         // will be removed in future
#define kQuad        2
#define kSquare      2         // will be removed in future


class SHAPE
{
  private:
    static int      nShape;
    static SHAPE*   pShape;

    static double loc_1D[4][4];         // location of 1D-GP (1, 3, 5, 7th-order)
    static double weight_1D[4][4];      // weight of Gauss points

    static double loc_TRI3[4][2];       // location of triangular GP (3rd-order)
    static double weight_TRI3[4];

    static double loc_TRI5[7][2];       // location of triangular GP (5th-order)
    static double weight_TRI5[7];

    static double loc_TRI7[13][2];      // location of triangular GP (7th-order)
    static double weight_TRI7[13];

    static void     gp();
    static void     types( int );

  public:
    int    dim;                         // dimension
    int    degree;                      // degree of GP integration
    int    ident;                       // identifier of shape function type
    int    order;                       // 1: linear, 2: quadratic; 3: bubble
    int    nnd;                         // number of nodes
    int    ngp;                         // number of Gauss points

    double f[kMaxGP2D][kMaxND];         // value of shape function at GP
    double dfdx[kMaxGP2D][kMaxND];      // derivatives of shape function at GP
    double dfdy[kMaxGP2D][kMaxND];
    double weight[kMaxGP2D];            // weight of Gauss points

    double dndx[kMaxND][kMaxND];        // derivatives of shape function at nodes
    double dndy[kMaxND][kMaxND];

  public:
   // Shape.cpp ----------------------------------------------------------------------------
    SHAPE();
    ~SHAPE();

    static void   init( int );
    static SHAPE* get( int, int );

    void   getCornerNodes( int, int*, int* );
    double jacobi2D( int,double *,double *,double *,double *,double [2][2] );

  private:
    // line.cpp ----------------------------------------------------------------------------
    void   shapeOfLine();
    void   localLine( int, double * );

    double lLine( int, double );                      // compute value of linear shape
    double lLineDx( int );                            // and its derivatives

    double qLine( int, double );                      // compute value of quadratic shape
    double qLineDx( int,double );                     // and its derivatives

    // triangle.cpp ------------------------------------------------------------------------
    void   shapeOfTriangle();
    int    localTriangle( int, double*, double* );    // compute local coordinates

    double lTriangle( int, double, double );          // linear shape
    double lTriangleDx( int );                        // and derivatives
    double lTriangleDy( int );

    double bTriangle( int, double, double );          // bubble shape (4 node triangle)
    double bTriangleDx( int, double, double );        // and derivatives
    double bTriangleDy( int, double, double );

    double qTriangle( int, double, double);           // quadratic shape
    double qTriangleDx( int, double, double);         // and derivatives
    double qTriangleDy( int, double, double);

    // square.cpp --------------------------------------------------------------------------
    void   shapeOfSquare();
    int    localSquare( int, double*, double* );      // compute local coordinates

    double lSquare( int, double, double );            // compute value of linear shape
    double lSquareDx( int, double );                  // and its derivatives
    double lSquareDy( int, double );

    double qSquare( int, double, double );            // compute value of quadratic shape
    double qSquareDx( int, double, double );          // and its derivatives
    double qSquareDy( int, double, double );
};

#endif
