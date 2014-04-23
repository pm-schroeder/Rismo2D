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

#include "Node.h"
#include "Elem.h"
#include "Shape.h"

#include "Grid.h"


void GRID::AreaFactors()
{
  double x[4], y[4], z[4];

  for( int e=0; e<ne; e++ )
  {
    ELEM* el  = &elem[e];
    int   ncn = el->Getncn();

    if( isFS(el->flag,ELEM::kRegion)  &&  el->GetLShape()->dim == 2 )
    {
      // compute local coordinates ---------------------------------------------------------

      x[0] = el->nd[0]->x;
      y[0] = el->nd[0]->y;
      z[0] = el->nd[0]->z;

      for( int i=1; i<ncn; i++ )
      {
        x[i] = el->nd[i]->x - x[0];
        y[i] = el->nd[i]->y - y[0];
        z[i] = el->nd[i]->z - z[0];
      }


      // compute area and horizontal projected area of element -----------------------------

      double tmp1 = y[1] * z[2]  -  y[2] * z[1];
      double tmp2 = z[1] * x[2]  -  z[2] * x[1];
      double tmp3 = x[1] * y[2]  -  x[2] * y[1];

      double areaElem = sqrt( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 );
      double areaPro  = tmp3;


      // squares are treated as two triangles ----------------------------------------------

      if( ncn == 4 )
      {
        tmp1 = y[2] * z[3]  -  y[3] * z[2];
        tmp2 = z[2] * x[3]  -  z[3] * x[2];
        tmp3 = x[2] * y[3]  -  x[3] * y[2];

        areaElem += sqrt( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 );
        areaPro  += tmp3;
      }

      el->areaFact = areaElem / areaPro;
    }

    else
    {
      el->areaFact = 1.0;
    }
  }
}
