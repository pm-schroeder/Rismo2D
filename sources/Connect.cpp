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
// Determine all elements connected to nodes.
// Only kRegion elements are taken into consideration.
//
// September 1992, M. Schroeder
// March     1993, M. Schroeder
// ---------------------------------------------------------------------------------------
//
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
