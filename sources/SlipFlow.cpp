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
