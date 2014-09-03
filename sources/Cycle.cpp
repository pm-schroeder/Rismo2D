// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class PROJECT
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

#include "Project.h"


int PROJECT::NextCycle( BCONSET* bconSet, int reset )
{
  static int actualCycle = 0;

  if( reset ) actualCycle = 0;


  // determine the next cycle ------------------------------------------------------------
  int next = 0;

  for( ;; )
  {
    next = bconSet->cycle[actualCycle];

    if( flowFlag )  next = abs(next);

    if( next >= 0 )
    {
      break;
    }

    else
    {
      if( actualCycle < kMaxCycles-1 )  actualCycle++;
      else                              bconSet->cycle[actualCycle] = 0;
    }
  }

  // determine the actual turbulence model -----------------------------------------------
  actualTurb = bconSet->turb[actualCycle];

  // determine the actual stationary settings --------------------------------------------
  actualStat = bconSet->stat[actualCycle];

  // determine the actual turbulence model -----------------------------------------------
  actualDisp = bconSet->disp[actualCycle];

  // determine the actual number of iterations per cycle ---------------------------------
  actualCycit = bconSet->cycit[actualCycle];
  if( actualCycit < 1 ) actualCycit = 1;

  // determine the actual solver ---------------------------------------------------------
  actualSolver = SOLVER::Getno( bconSet->solver[actualCycle] );

  // next == 0  ... stop iteration -------------------------------------------------------
  if( next == 0 )
  {
    return 0;
  }

  // else       ... get next cycle from cycle-array --------------------------------------
  else
  {
    if( actualCycle < kMaxCycles-1 )  actualCycle++;
    else                              bconSet->cycle[actualCycle] = 0;
  }

  theCycle = next;

  return theCycle;
}


void PROJECT::PrintTheCycle( int iter )
{
  char ltxt[120];
  char rtxt[120];

  if( REPORT::rpt.level < 1  &&  iter > 1 ) return;

  REPORT::rpt.Message( 1, "\n\n" );
  REPORT::rpt.Line1( 1 );

  sprintf( ltxt, "time step %d", iTM );
  sprintf( rtxt, "" );

  switch( theCycle )
  {
    case kUVSCyc:
      if( iter == 0 )
      {
        sprintf( rtxt, "prediction step in UVS cycle (Taylor-Hood)" );
      }
      else if( isFS(actualTurb, BCONSET::kVtAnisotrop) )
      {
        sprintf( rtxt, "%d. iteration in anisotrop UVS cycle (Taylor-Hood)", iter );
      }
      else
      {
        sprintf( rtxt, "%d. iteration in UVS cycle (Taylor-Hood)", iter );
      }
      break;

    case kUVS_TMCyc:
      if( iter == 0 )
      {
        sprintf( rtxt, "prediction step in unsteady UVS cycle (Taylor-Hood)" );
      }
      else if( isFS(actualTurb, BCONSET::kVtAnisotrop) )
      {
        sprintf( rtxt, "%d. iteration in unsteady anisotrop UVS cycle (Taylor-Hood)", iter );
      }
      else
      {
        sprintf( rtxt, "%d. iteration in unsteady UVS cycle (Taylor-Hood)", iter );
      }
      break;

    case kUVS_MECyc:
      sprintf( rtxt, "%d. iteration in UVS cycle (MINI)", iter );
      break;

    case kUVS_ME_TMCyc:
      if( iter == 0 )
      {
        sprintf( rtxt, "prediction step in unsteady UVS cycle (MINI)" );
      }
      else if( isFS(actualTurb, BCONSET::kVtAnisotrop) )
      {
        sprintf( rtxt, "%d. iteration in unsteady anisotrop UVS cycle (MINI)", iter );
      }
      else
      {
        sprintf( rtxt, "%d. iteration in unsteady UVS cycle (MINI)", iter );
      }
      break;

    case kDispCurv2D:
      sprintf( rtxt, "dispersion due to secondary flow (curvature)" );
      break;

    case kKDInitCyc:
      sprintf( rtxt, "KD initialization" );
      break;

    case kKLCyc:
      sprintf( rtxt, "%d. iteration in KL cycle", iter );
      break;

    case kKDCyc:
    case kKD_LCyc:
    case kKD_QCyc:
      sprintf( rtxt, "%d. iteration in KD cycle", iter );
      break;

    case kKCyc:
      sprintf( rtxt, "%d. iteration in K cycle", iter );
      break;

    case kDCyc:
      sprintf( rtxt, "%d. iteration in D cycle", iter );
      break;

    case kQbCyc:
      sprintf( rtxt, "%d. iteration in sediment transport cycle", iter );
      break;

    case kSLCyc:
      sprintf( rtxt, "%d. iteration in suspended load cycle", iter );
      break;

    case kBLCyc:
      sprintf( rtxt, "%d. iteration in bed load cycle", iter );
      break;

    case kBSLCyc:
      sprintf( rtxt, "%d. iteration in suspended + bed load cycle", iter );
      break;

    case kDiffBLCyc:
      sprintf( rtxt, "load difference cycle" );
      break;

    case kDivCyc:
      sprintf( rtxt, "%d. iteration in zero-div cycle", iter );
      break;

    case kDryRewet:
      sprintf( rtxt, "dry/rewet cycle");
      break;

    case kReOrderCyc:
      sprintf( rtxt, "reorder cycle");
      break;

    case kSurfaceCyc:
      sprintf( rtxt, "initialization of free surface");
      break;

    case kSmoothS:
      sprintf( rtxt, "smoothing cycle - free surface");
      break;

    case kCornerToVol:
      sprintf( rtxt, "free surface to volumes cycle");
      break;

    case kScaledOutputCyc:
      sprintf( rtxt, "scaled output cycle");
      break;

    case kOutputCyc:
      sprintf( rtxt, "output cycle");
      break;
  }

  REPORT::rpt.Message( 1, " %-26s  %50s\n", ltxt, rtxt );
}
