// ======================================================================================
//
// Copyright (C) 1992-2010  by  P.M. SCHROEDER
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
#include "Bcon.h"
#include "Node.h"

#include "Grid.h"


void GRID::SetSlipFlow()
{
  char text[100];

  sprintf( text, "\n (GRID::SetSlipFlow)     %s\n",
                 "slip flow condition applied to boundary nodes" );
  REPORT::rpt.Output( text, 4 );


  /* --- rotate velocity vector to make perpendicular component zero ----- */

  for( int i=0; i<np; i++ )
  {
    if( isFS(node[i].flag, NODE::kRotat) )
    {
      double U = node[i].v.U;
      double V = node[i].v.V;

      double dUdt = node[i].v.dUdt;
      double dVdt = node[i].v.dVdt;


      // rotate velocity vector ----------------------------------------------------------

      double Ut = node[i].bc.Getrot(0,0)*U + node[i].bc.Getrot(1,0)*V;

      node[i].v.U = node[i].bc.Getrot(0,0) * Ut;
      node[i].v.V = node[i].bc.Getrot(1,0) * Ut;


      // rotate time derivatives ---------------------------------------------------------

      Ut = node[i].bc.Getrot(0,0)*dUdt + node[i].bc.Getrot(1,0)*dVdt;

      node[i].v.dUdt = node[i].bc.Getrot(0,0) * Ut;
      node[i].v.dVdt = node[i].bc.Getrot(1,0) * Ut;
    }
  }
}
