// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class SECTION
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
#include "Report.h"

#include "Section.h"


SECTION::SECTION()
{
}


SECTION::~SECTION()
{
}


void SECTION::init( double xs, double ys, double xe, double ye, double z )
{
  this->elem = (ELEM*) 0;

  this->xs = xs;
  this->xe = xe;
  this->ys = ys;
  this->ye = ye;
  this->z  = z;

  tx = xe - xs;
  ty = ye - ys;

  double L  = sqrt( tx * tx  +  ty * ty );
  if( fabs(L) < 1.0e-6 )
    REPORT::rpt.Error( "zero length between section points - SECTION::init (1)" );

  tx /= L;
  ty /= L;

  nx = -ty;
  ny = tx;
}


double SECTION::getxs() { return xs; }
double SECTION::getys() { return ys; }
double SECTION::getxe() { return xe; }
double SECTION::getye() { return ye; }


double SECTION::length()
{
  return sqrt( (xe-xs) * (xe-xs)  +  (ye-ys) * (ye-ys) );
}


double SECTION::normalDistance( double x, double y )
{
  return (x-xs) * nx  +  (y-ys) * ny;
}


double SECTION::tangentialDistance( double x, double y )
{
  return (x-xs) * tx  +  (y-ys) * ty;
}


int SECTION::interpolate( SECTION* next,
                          double   x,
                          double   y,
                          double*  z )
{
  double epsilon = 1.0e-3;

  // check if coordinate (x,y) is inside of the quadrilateral
  // built of sections "this" and "next", and interpolate z

  SECTION left;
  SECTION rght;

  left.init( this->xs, this->ys, next->xs, next->ys, 0.0 );
  rght.init( this->xe, this->ye, next->xe, next->ye, 0.0 );

  double dThis = this->normalDistance( x, y );
  double dNext = next->normalDistance( x, y );

  double dLeft = left.normalDistance( x, y );
  if( left.normalDistance(this->xe, this->ye) > 0.0 ) dLeft *= -1.0;

  double dRght = rght.normalDistance( x, y );
  if( rght.normalDistance(this->xs, this->ys) < 0.0 ) dRght *= -1.0;

  if(     dThis > -epsilon  &&  dNext <  epsilon
      &&  dLeft <  epsilon  &&  dRght > -epsilon  )
  {
    *z = (-this->z * dNext  +  next->z * dThis) / (dThis - dNext);
    return 1;
  }

  return 0;
}


void SECTION::scaleFR( double lScale, double hScale )
{
  this->xs /= lScale;
  this->xe /= lScale;
  this->ys /= lScale;
  this->ye /= lScale;
  this->z  /= hScale;
}
