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
