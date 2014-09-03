// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class ELEM
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

#include "Defs.h"
#include "Shape.h"
#include "Node.h"
#include "Elem.h"


ELEM::ELEM()
{
  mark  = false;
  rigid = false;
  shape = 0;
  flag  = 0;
}


ELEM::~ELEM()
{
}


ELEM ELEM::operator =( const ELEM& e )
{
  no       = e.no;
  name     = e.name;
  link     = e.link;
  flag     = e.flag;
  shape    = e.shape;
  type     = e.type;
  areaFact = e.areaFact;
  U        = e.U;
  V        = e.V;
  P        = e.P;
  dz       = e.dz;
//  lShape   = e.lShape;
//  bShape   = e.bShape;
//  qShape   = e.qShape;

  for( int i=0; i<kMaxNodes2D; i++ )  nd[i] = e.nd[i];

  for( int i=0; i<kSimDF; i++ )
  {
    for( int j=0; j<kMaxNodes2D; j++ )  isLast[i][j] = e.isLast[i][j];
  }

  return *this;
}


void ELEM::Setshape( int id )
{
  switch( id )
  {
    case kLine:
    case kTriangle:
    case kSquare:
      shape = id;
      break;

    default:
      REPORT::rpt.Error( kParameterFault, "element type %d not supported %s",
                                          "(ELEM::setshape)", id );
  }
//  switch( id )
//  {
//    case kLine:
//                    lShape = SHAPE::get( id, 2 );
//                    bShape = NULL;
//                    qShape = SHAPE::get( id, 3 );   break;
//
//    case kTriangle: lShape = SHAPE::get( id, 3 );
//                    bShape = SHAPE::get( id, 4 );
//                    qShape = SHAPE::get( id, 6 );   break;
//
//    case kSquare:   lShape = SHAPE::get( id, 4 );
//                    bShape = SHAPE::get( id, 5 );
//                    qShape = SHAPE::get( id, 8 );   break;
//
//    default:
//      REPORT::rpt.Error( kParameterFault, "element type %d not supported %s",
//                                "(ELEM::setshape)", id );
//  }
}


SHAPE* ELEM::GetLShape()   { return SHAPE::get(shape,1); }
SHAPE* ELEM::GetQShape()   { return SHAPE::get(shape,2); }
SHAPE* ELEM::GetBShape()   { return SHAPE::get(shape,3); }
int    ELEM::GetshapeID()  { return GetLShape()->ident; }
int    ELEM::Getncn()      { return GetLShape()->nnd; }
int    ELEM::Getnnd()      { return GetQShape()->nnd; }


NODE* ELEM::Getnode( int i )
{
  return nd[i];
}


void ELEM::Center( double* xc, double* yc )
{
  int nnd = 0;

  *xc = *yc = 0.0;

  switch( GetshapeID() )
  {
    case kLine:       nnd = 2;  break;
    case kTriangle:   nnd = 3;  break;
    case kSquare:     nnd = 4;  break;
  }

  for( int i=0; i<nnd; i++ )
  {
    NODE* nd = Getnode(i);

    *xc += nd->x;
    *yc += nd->y;
  }

  *xc /= nnd;
  *yc /= nnd;
}


void ELEM::Center( double* xc, double* yc, double* zc )
{
  int nnd = 0;

  *xc = *yc = *zc = 0.0;

  switch( GetshapeID() )
  {
    case kLine:       nnd = 2;  break;
    case kTriangle:   nnd = 3;  break;
    case kSquare:     nnd = 4;  break;
  }

  for( int i=0; i<nnd; i++ )
  {
    NODE* nd = Getnode(i);

    *xc += nd->x;
    *yc += nd->y;
    *zc += nd->z;
  }

  *xc /= nnd;
  *yc /= nnd;
  *zc /= nnd;
}


double ELEM::area()
{
  int    nnd = 0;
  double A = 0.0;

  switch( GetshapeID() )
  {
    case kLine:       nnd = 0;  break;
    case kTriangle:   nnd = 3;  break;
    case kSquare:     nnd = 4;  break;
  }


  for( int i=0; i<nnd; i++ )
  {
    int j = (i + 1) % nnd;

    NODE* ndi = Getnode(i);
    NODE* ndj = Getnode(j);

    A += ( ndi->x - ndj->x ) * ( ndi->y + ndj->y );
  }

  return A / 2.0;
}


double ELEM::perimeter()
{
  int    nnd = 0;
  double P = 0.0;

  switch( GetshapeID() )
  {
    case kLine:       nnd = 0;  break;
    case kTriangle:   nnd = 3;  break;
    case kSquare:     nnd = 4;  break;
  }


  for( int i=0; i<nnd; i++ )
  {
    int   j   = (i + 1) % nnd;

    NODE* ndi = Getnode(i);
    NODE* ndj = Getnode(j);

    double dx = ndj->x - ndi->x;
    double dy = ndj->y - ndi->y;

    P += sqrt( dx*dx + dy*dy );
  }

  return P;
}
