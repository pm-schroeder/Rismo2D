// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class MODEL
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
#include "Memory.h"
#include "Node.h"
#include "Elem.h"

#include "Model.h"


// ---------------------------------------------------------------------------------------
// determine last occurence of nodes during element assembling

void MODEL::LastNode()
{
  GRID* rg = region;


  // allocate temporary memory for counter vector and initialize -------------------------

  int* counter = (int*) MEMORY::memo.Array_nd( rg->Getnp() );

  for( int i=0; i<rg->Getnp(); i++ )  counter[i] = 0;


  // initialize the elem->isLast array ---------------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = elem[e];

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )  el->isLast[0][i] = 0;
  }


  // determine the number of elements associated to nodes ----------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = elem[e];

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      int no = el->nd[i]->Getno();

      counter[no]++;
    }
  }


  // determine the last occurence of nodes -----------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = elem[e];

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      int no = el->nd[i]->Getno();

      counter[no]--;

      if( counter[no] <  0 )
        REPORT::rpt.Error( kUnexpectedFault, "wrong count of nodes (MODEL::lastNode - 2)" );

      if( counter[no] == 0 )  el->isLast[0][i] = 1;
    }
  }

  MEMORY::memo.Detach( counter );


  char text[100];
  sprintf( text, "\n%-25s%s\n", " (MODEL::LastNode)",
                 "last occurence of nodes determined" );
  REPORT::rpt.Output( text, 2 );
}
