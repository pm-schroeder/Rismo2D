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
        sprintf( rtxt, "prediction step in UVS cycle" );
      }
      else if( isFS(actualTurb, BCONSET::kVtAnisotrop) )
      {
        sprintf( rtxt, "%d. iteration in anisotrop UVS cycle", iter );
      }
      else
      {
        sprintf( rtxt, "%d. iteration in UVS cycle", iter );
      }
      break;

    case kUVS_TMCyc:
      if( iter == 0 )
      {
        sprintf( rtxt, "prediction step in unsteady UVS cycle" );
      }
      else if( isFS(actualTurb, BCONSET::kVtAnisotrop) )
      {
        sprintf( rtxt, "%d. iteration in unsteady anisotrop UVS cycle", iter );
      }
      else
      {
        sprintf( rtxt, "%d. iteration in unsteady UVS cycle", iter );
      }
      break;

    case kUVS_LVCyc:
      sprintf( rtxt, "%d. iteration in UVS_LV cycle", iter );
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

    case kSurfaceToVol:
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
