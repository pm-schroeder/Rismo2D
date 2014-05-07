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
