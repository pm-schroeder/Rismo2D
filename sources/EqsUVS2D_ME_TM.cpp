// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_ME_TM
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
#include "Vars.h"
#include "Shape.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsUVS2D_ME_TM.h"


//#define kDebug


EQS_UVS2D_ME_TM::EQS_UVS2D_ME_TM() : EQS_UVS2D_ME()
{
  neq = 0;
}


EQS_UVS2D_ME_TM::~EQS_UVS2D_ME_TM()
{
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_ME_TM::Timegrad( PROJECT* project, double dt, double th )
{
  // call the base class ---------------------------------------------------------------------------
  EQS_UVS2D::Timegrad( project, dt, th );
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_ME_TM::Predict( PROJECT* project, int init, double dt, double th )
{
  // call the base class ---------------------------------------------------------------------------
  EQS_UVS2D::Predict( project, init, dt, th );

  // ---------------------------------------------------------------------------------------------
  // initialize U,V on previously dry elements to their current values

  for( int e=0; e<project->M2D->ne; e++ )
  {
    ELEM *el = project->M2D->elem[e];

    if( project->GetTimeStep() == project->timeint.firstTimeStep || isFS(el->flag, ELEM::kDryPrev) )
    {
      el->Uo = el->U;
      el->Vo = el->V;
    }
  }
  return;
}
