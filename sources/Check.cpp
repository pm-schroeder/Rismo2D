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

#include "Defs.h"
#include "Report.h"
#include "Shape.h"
#include "Node.h"
#include "Elem.h"

#include "Grid.h"


void GRID::Check()
{
  char   text[80];
  double x[kMaxNodes2D], y[kMaxNodes2D];
  double trafo2[2][2];

  long   errFlag = 0;
  double detj    = 0.0;
  double area    = 0.0;

  int    count[2];

  count[0] = count[1] = 0;


  /* -------------------------------------------------------------------------------- */

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = &elem[e];

    el->mark = false;

    int nnd = el->Getnnd();

    x[0] = el->nd[0]->x;
    y[0] = el->nd[0]->y;

    for( int i=1; i<nnd; i++ )
    {
      x[i] = el->nd[i]->x - *x;
      y[i] = el->nd[i]->y - *y;
    }

    x[0] = y[0] = 0.0;

    for( int i=0; i<el->GetQShape()->ngp; i++ )
    {
      switch( el->GetQShape()->dim )
      {
        case 2:
          detj = el->GetQShape()->jacobi2D ( nnd,
                                             el->GetQShape()->dfdx[i],
                                             el->GetQShape()->dfdy[i],
                                             x, y, trafo2 );
          if ( detj <= 0.0  )  SF( errFlag, 4l );
          break;
      }


      area += detj * el->GetQShape()->weight[i];
    }

    if( isFS(errFlag, 4l)  )
    {
      CF( errFlag, 4l );
      SF( errFlag, 1l );
      count[0]++;
      el->mark = true;
    }

    if( area <= 0.0 )
    {
      SF( errFlag, 2l );
      count[1]++;
      el->mark = true;
    }
  }


  if( errFlag & 1l )
  {
    sprintf( text, "\n (check)         warning: connectivity error in %d elements:\n",
                   count[0] );
    REPORT::rpt.Output( text );
  }

  if( errFlag & 2l )
  {
    sprintf( text, "\n (check)         warning: element area A <= 0 in %d elements:\n",
                   count[1] );
    REPORT::rpt.Output( text );
  }


  if( errFlag )
  {
    int i = 0;

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( el->mark )
      {
        sprintf( text, "  %5d", el->Getname() );
        REPORT::rpt.Output( text );

        i++;
        if( !( i % 10) )  REPORT::rpt.Output( "\n" );
      }
    }

    if( i % 10 )  REPORT::rpt.Output( "\n" );
  }
}
