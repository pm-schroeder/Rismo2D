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

#include "Grid.h"


double GRID::ReportCuPe( double dt, double vk )
{
  double maxCu_x, maxCu_y;
  double maxPe_x, maxPe_y;


  // report Courant- and Peclet-Number: Cu and Pe ------------------------------------------

  int*    counter = (int*)    MEMORY::memo.Array_nd( np );
  double* xnd     = (double*) MEMORY::memo.Array_nd( np );
  double* ynd     = (double*) MEMORY::memo.Array_nd( np );

  for( int i=0; i<np; i++ )
  {
    counter[i] = 0;
    xnd[i]     = 0.0;
    ynd[i]     = 0.0;
  }


  for( int e=0; e<ne; e++ )
  {
    ELEM* el = &elem[e];

    int    ncn, nnd;
    double x, y, xe, ye;

    if( isFS( el->flag, ELEM::kRegion ) )
    {
      ncn = el->Getncn();
      nnd = el->Getnnd();

      xe = ye = 0.0;

      for( int i=1; i<ncn; i++ )
      {
        x = fabs( el->nd[i]->x - el->nd[0]->x );
        y = fabs( el->nd[i]->y - el->nd[0]->y );

        if( x > xe ) xe = x;
        if( y > ye ) ye = y;
      }

      for( int i=0; i<nnd; i++ )
      {
        int no = el->nd[i]->Getno();

        if( no < 0  ||  no >= np )
          REPORT::rpt.Error( kRangeFault, "node number is out of range (GRID::reportCuPe - 1)" );

        counter[no]++;

        xnd[no] += xe;
        ynd[no] += ye;
      }
    }
  }


  maxCu_x = 0.0;      maxCu_y = 0.0;
  maxPe_x = 0.0;      maxPe_y = 0.0;

  for( int i=0; i<np; i++ )
  {
    double Cu, Pe, U, V;

    if( isFS(node[i].flag, NODE::kDry) )  continue;

    if( counter[i] > 0 )
    {
      xnd[i] /= counter[i];
      ynd[i] /= counter[i];
    }

    U = node[i].v.U;
    V = node[i].v.V;

    Cu = fabs( U * dt / xnd[i] );
    Pe = fabs( U * xnd[i] / (node[i].vt + vk) );

    if( Cu > maxCu_x ) maxCu_x = Cu;
    if( Pe > maxPe_x ) maxPe_x = Pe;

    Cu = fabs( V * dt / ynd[i] );
    Pe = fabs( V * ynd[i] / (node[i].vt + vk) );

    if( Cu > maxCu_y ) maxCu_y = Cu;
    if( Pe > maxPe_y ) maxPe_y = Pe;
  }


  REPORT::rpt.Message( 1, "\n\n%-25s%s\n",
                          " (GRID::ReportCuPe)", "maximum of       CU           Pe");

  REPORT::rpt.Message( 1, " %s %9.1le     %9.1le\n",
                          "                        x-direction: ", maxCu_x, maxPe_x );
  REPORT::rpt.Message( 1, " %s %9.1le     %9.1le\n",
                          "                        y-direction: ", maxCu_y, maxPe_y );


  MEMORY::memo.Detach( counter );
  MEMORY::memo.Detach( xnd );
  MEMORY::memo.Detach( ynd );


  if( maxCu_x > maxCu_y )  return maxCu_x;
  else                     return maxCu_y;
}
