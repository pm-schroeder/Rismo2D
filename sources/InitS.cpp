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
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Project.h"
#include "Section.h"

#include "Grid.h"


void GRID::InitS( int nSection, SECTION* section )
{
  for( int i=0; i<np; i++ )
  {
    NODE* nd = &node[i];

    if( nSection == 1 )
    {
      nd->v.S = section[0].z;
    }

    else
    {
      for( int j=0; j<nSection-1; j++ )
      {
        double S;

        if( section[j].interpolate(&section[j+1], nd->x, nd->y, &S) )
        {
          nd->v.S = S;
          break;
        }
      }
    }
  }
}
