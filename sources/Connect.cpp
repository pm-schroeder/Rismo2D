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


void GRID::Connection( long dry )
{
  static int    nElemBuf = 0;
  static ELEM** elemBuf  = (ELEM**) 0;


  // initializations ---------------------------------------------------------------------

  for( int n=0; n<np; n++ )  node[n].noel = 0;


  // count number of connected elements per node -----------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = &elem[e];

    if( isFS(el->flag, ELEM::kRegion)  &&  !isFS(el->flag, dry) )
    {
      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        int j = el->nd[i]->Getno();
        node[j].noel++;
      }
    }
  }


  // allcate dynamic memory for element connections --------------------------------------

  int elcnt = 0;
  for( int n=0; n<np; n++ ) elcnt += node[n].noel;


  if( elcnt > nElemBuf )
  {
    if( nElemBuf ) delete[] elemBuf;

    elemBuf = new ELEM* [elcnt];
    if( !elemBuf )
      REPORT::rpt.Error( kMemoryFault, "can not allocate memory (GRID::connection - 1)" );

    nElemBuf = elcnt;
  }

  elcnt = 0;
  for( int n=0; n<np; n++ )
  {
    if( node[n].noel > 0 )
    {
      node[n].el = elemBuf + elcnt;
      elcnt += node[n].noel;
    }

    node[n].noel = 0;
  }


  // set up connection -------------------------------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = &elem[e];

    if( isFS(el->flag,ELEM::kRegion) && !isFS(el->flag,dry) )
    {
      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        int j = el->nd[i]->Getno();
        int k = node[j].noel;

        node[j].el[k] = (ELEM *) el;       // set pointer to element
        node[j].noel++;                    // increase number of connected elements
      }
    }
  }


  REPORT::rpt.Output( "\n (GRID::connection)      element connectivity set up\n", 4 );
}
