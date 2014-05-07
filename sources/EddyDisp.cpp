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
#include "Type.h"
#include "Node.h"
#include "Elem.h"

#include "Grid.h"

//define DEBUG 1


void GRID::EddyDisp()
{
  char text[80];


  // allocate temporary memory for counter vector ----------------------------------------

  int* counter = (int*) MEMORY::memo.Array_nd( np );


  // initialization ----------------------------------------------------------------------

  for( int n=0; n<np; n++ )
  {
    counter[n]   = 0;
    node[n].exx = 0.0;
    node[n].exy = 0.0;
    node[n].eyy = 0.0;
  }

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = &elem[e];

    // eddy difusivity (est) is computed from eddy viscosity (vt) with the
    // turbulent Schmidt number (proportional factor st)
    //                         est = vt / st
    // since a different parameter KDest was used for computing the eddy
    // viscosity, Exx and Eyy have to be adapted for multiplication with vt.

    TYPE* type = TYPE::Getid( el->type );

    double ExxInit = type->estx / type->st / type->KDest;
    double EyyInit = type->esty / type->st / type->KDest;

    NODE** nd  = el->nd;
    int    nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      nd[i]->exx += ExxInit;
      nd[i]->eyy += EyyInit;

      counter[nd[i]->Getno()]++;
    }
  }


  sprintf( text, "\n%-25s%s\n", " (GRID::EddyDisp)", "computing eddy diffusity" );
  REPORT::rpt.Output( text, 4 );


  // compute eddy diffusity from eddy viscosity-------------------------------------------

  for( int n=0; n<np; n++ )
  {
    if( counter[n] )
    {
      node[n].exx /= counter[n];
      node[n].eyy /= counter[n];
    }

//  changed: multiplication with eddy viscosity is now performed in coefsXYZ()
//  29. April 2004, SC
//  node[n].Exx *= node[n].vt;
//  node[n].Eyy *= node[n].vt;
  }


  // transform local Exx,Eyy into global Exx,Exy,Eyy -------------------------------------

  for( int n=0; n<np; n++ )
  {
    // angle between Ures and global x-Axis ----------------------------------------------

    double U     = node[n].v.U;
    double V     = node[n].v.V;
    double Ur    = sqrt( U*U + V*V );

    double sinus = 0.0;
    double cosin = 1.0;

    if( Ur > 1.0e-6 )
    {
      sinus = V / Ur;
      cosin = U / Ur;
    }


    // transformation ------------------------------------------------------------------------

    double Ex = node[n].exx;
    double Ey = node[n].eyy;

    node[n].exx =  Ex * cosin * cosin  +  Ey * sinus * sinus;
    node[n].exy = (Ex - Ey) * sinus * cosin;
    node[n].eyy =  Ex * sinus * sinus  +  Ey * cosin * cosin;
  }


  // free temporary used memory --------------------------------------------------------------

  MEMORY::memo.Detach( counter );


#ifdef DEBUG
{
  FILE *id = fopen( "eddyDisp.Report", "w" );

  if( id )
  {
    for( int n=0; n<np; n++ )
    {
      fprintf(id, "%5d   %12.4le  %12.4le  %12.4le\n",
                  node[n].Getname(),
                  node[n].exx,
                  node[n].exy,
                  node[n].eyy );
    }

    fclose( id );
  }
}
#endif
}
