// ======================================================================================
//
// Copyright (C) 1992-2013  by  P.M. SCHROEDER
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
#include "Asciifile.h"
#include "Report.h"
#include "Scale.h"
#include "Shape.h"
#include "Elem.h"
#include "Node.h"
#include "Grid.h"
#include "Model.h"
#include "Type.h"
#include "Memory.h"
#include "Section.h"
#include "Statist.h"
#include "Time.h"

#include "Project.h"

#include <sys/stat.h>
#include <sys/types.h>


void PROJECT::Compute()
{
  // read project data from input files --------------------------------------------------
  Input();                                      // common data

  if( name.timeStepFile[0] )                    // boundary conditions
  {
    timeint.Input( name.timeStepFile );
  }
  else
  {
    timeint.Input( name.inputFile );
  }

  // common initializations --------------------------------------------------------------
  REPORT::rpt.InitTheClock();

  // read data from geometry files -------------------------------------------------------
  REPORT::rpt.Message( 2, "\n" );
  REPORT::rpt.PrintTime( 2 );
  REPORT::rpt.Message( 3, "\n\n reading model data...\n\n" );

  M2D = new MODEL;
  if( !M2D )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (PROJECT::Compute #1)" );

  M2D->subdom = &subdom;

  int initial = M2D->Input( this, &name, &timeint.prevTime, &sed );

  GRID* R2D = M2D->region;

  // read statistical values -------------------------------------------------------------
  if( statistics  )
  {
    statist = new STATIST;

    statist->Init( R2D->Getnp() );
    statist->Read( R2D->Getnp(), name.inputStatistFile, &subdom );
  }

  // read nodal values of previous time step, if specified -------------------------------
  R2D->InitPrevious();

  // scale all variables -----------------------------------------------------------------
  switch( scale.gettype() )
  {
    case 1:
    case 2:
      scale.scale( this );
      break;
  }

  // check connectivity and element area -------------------------------------------------
  R2D->Check();

  // set list of gauges from boundary conditions -----------------------------------------
  int maxgauge = 10;
  gauge  = new int[maxgauge];
  ngauge = 0;

  for( int i=0; i<maxgauge; i++ ) gauge[i] = 0;

  for( int i=0; i<timeint.setsOfBcon; i++ )
  {
    BCONSET *bcs= &timeint.bconSet[i];
    for( int j=0; j<bcs->nofLines; j++ )
    {
      int gnd = bcs->bcLine[j].gct->nocg;
      if( gnd > 0 )
      {
        for( int k=0; k<ngauge; k++ )
        {
          if( gnd == gauge[k] )
          {
            gnd = 0;
            break;
          }
        }

        if( gnd )
        {
          gauge[ngauge] = gnd;
          ngauge++;

          if( ngauge >= maxgauge )
          {
            maxgauge *= 2;
            int *newgauge = new int[maxgauge];

            int k = 0;
            for( ; k<ngauge;   k++ ) newgauge[k] = gauge[k];
            for( ; k<maxgauge; k++ ) newgauge[k] = 0;

            delete[] gauge;
            gauge = newgauge;
          }
        }
      }
    }
  }

  gaugeS  = new double[maxgauge];
  gaugeSo = new double[maxgauge];

  for( int i=0; i<maxgauge; i++ )
  {
    gaugeS[i]  = 0.0;
    gaugeSo[i] = 0.0;
  }

  // initialize specified boundary conditions --------------------------------------------
  timeint.bconSet[0].InitBcon( this, &timeint.actualTime );

  // set location of boundary elements and nodes -----------------------------------------
  M2D->SetLocation();

  // continuity check of initial conditions (if specified) -------------------------------
  if( initial ) M2D->Continuity();

  // compute area factors to increase the friction coefficient ---------------------------
  // for inclined two-dimensional elements
  R2D->AreaFactors();


  // -------------------------------------------------------------------------------------
  // iteration loop for time integration

  int bcSetNo     = -1;
  int prevBcSetNo = -1;

  MEMORY::memo.PrintInfo();


  // store the report level into a local variable ----------------------------------------

  int report_level = REPORT::rpt.level;

  for( iTM=timeint.firstTimeStep; iTM<=timeint.lastTimeStep; iTM++ )
  {
    // -----------------------------------------------------------------------------------
    // check, if output should be written to the report file in this time step
    // if not, set the report level temporarily to 0

    if( REPORT::rpt.ntm > 0  &&  !REPORT::rpt.tmlist[iTM]
                             &&  !timeint.result[iTM - 1] )  REPORT::rpt.Setlevel(0);

    REPORT::rpt.Line2( 1 );


    // -----------------------------------------------------------------------------------
    // copy periodic boundary conditions for nodes

    for( int i=0; i<timeint.nPeriodicNode; i++ )
    {
      NODE* nd_src = M2D->region->Getnode(timeint.periodicNode[0][i]-1);
      NODE* nd_dst = M2D->region->Getnode(timeint.periodicNode[1][i]-1);

      if( isFS(nd_dst->bc.kind, BCON::kInlet) )
      {
        double H = nd_src->v.S - nd_src->z;

        if( H > 0.0 )
        {
          nd_dst->bc.val->U = -nd_src->v.U * H;
          nd_dst->bc.val->V = -nd_src->v.V * H;
        }
      }
    }


    // -----------------------------------------------------------------------------------
    // copy periodic boundary conditions for lines

    for( int i=0; i<timeint.nPeriodicLine; i++ )
    {
      static int    nct = 0;
      static NODE** src = NULL;
      static NODE** dst = NULL;

      // count number of edges in line [0] and line [2]
      int ns = 0;
      int nd = 0;

      for( int j=0; j<M2D->control->Getne(); j++ )
      {
        ELEM* ctrl = M2D->control->Getelem(j);

        if( ctrl->type == timeint.periodicLine[0][i] )  ns++;
        if( ctrl->type == timeint.periodicLine[2][i] )  nd++;

        ctrl->mark = 0;
      }

      if( ns > nct || nd > nct)
      {
        nct = (ns>nd)?(ns):(nd);

        if( src )  delete[] src;
        if( dst )  delete[] dst;

        src = new NODE* [3*nct];
        dst = new NODE* [3*nct];

        if( !src || !dst )
          REPORT::rpt.Error( kMemoryFault, "can not allocate memory (PROJECT::Compute #2)" );
      }

      src[0] = M2D->region->Getnode(timeint.periodicLine[1][i]-1);
      dst[0] = M2D->region->Getnode(timeint.periodicLine[3][i]-1);


      // find the two edges with node [1] and node [3] respectively

      int n = 0;

      for( ;; )
      {
        int found = 0;

        for( int j=0; j<M2D->control->Getne(); j++ )
        {
          ELEM* ctrl = M2D->control->Getelem(j);

          if( ctrl->type == timeint.periodicLine[0][i]  &&  !ctrl->mark )
          {
            if( ctrl->nd[0] == src[n] )
            {
              src[n+1] = ctrl->nd[2];
              src[n+2] = ctrl->nd[1];

              ctrl->mark = 1;
              found++;
            }

            else if( ctrl->nd[1] == src[n] )
            {
              src[n+1] = ctrl->nd[2];
              src[n+2] = ctrl->nd[0];

              ctrl->mark = 1;
              found++;
            }
          }

          if( ctrl->type == timeint.periodicLine[2][i]  &&  !ctrl->mark )
          {
            if( ctrl->nd[0] == dst[n] )
            {
              dst[n+1] = ctrl->nd[2];
              dst[n+2] = ctrl->nd[1];

              ctrl->mark = 1;
              found++;
            }
            
            else if( ctrl->nd[1] == dst[n] )
            {
              dst[n+1] = ctrl->nd[2];
              dst[n+2] = ctrl->nd[0];

              ctrl->mark = 1;
              found++;
            }
          }
        }

        if( found != 2 )  break;

        n += 2;
        if( n >= 3 * nct )  break;      // Indeed, this should not happen !
      }


      // copy values for U,V from outlet boundary to inlet boundary

      REPORT::rpt.Message( 1, " applying periodic boundary conditions (source / destination node ...\n" );

      for( int j=0; j<=n; j++ )
      {
        if( isFS(dst[j]->bc.kind, BCON::kInlet) )
        {
          double H = src[j]->v.S - src[j]->z;

          if( H > 0.0 )
          {
            dst[j]->bc.val->U = -src[j]->v.U * H;
            dst[j]->bc.val->V = -src[j]->v.V * H;

            REPORT::rpt.Message( 1, "  %6d > %6d\n", dst[j]->Getname(), src[j]->Getname() );
          }
        }
      }
    }


    // -----------------------------------------------------------------------------------
    // set times

    timeint.nextTime   = timeint.prevTime + timeint.deltaTime;
    timeint.incTime    = timeint.deltaTime;

    timeint.actualTime = timeint.prevTime;


    // -----------------------------------------------------------------------------------
    // choose number of boundary condition set

    int j;

    if( iTM > timeint.endLoop  &&  timeint.bcLoop > 0 )
    {
      int b = timeint.bcLoop;
      int n = 1 + timeint.endLoop - b;

      if( n == 0 )
      {
        b = timeint.endLoop;
        n = 1;
        REPORT::rpt.Warning( kParameterFault,
                   "Time step loop not specified; setting loop = %d", b );
      }
      j = b  +  (iTM - b) % n;
    }

    else
    {
      j = iTM;
    }

    bcSetNo = -1;

    for( int i=0; i<timeint.setsOfBcon; i++ )
    {
      if( timeint.bconSet[i].noTimeStep <= j )
      {
        bcSetNo = i;
      }

      else if( timeint.bconSet[i].noTimeStep > j )
      {
        break;
      }
    }

    if( bcSetNo < 0 )  REPORT::rpt.Error( "no boundary condition set found - main (3)" );

    timeint.actualBcSet = &timeint.bconSet[bcSetNo];


    // determine the next iteration cycle ------------------------------------------------

    NextCycle( &timeint.bconSet[bcSetNo], true );


    // prepare for time integration ------------------------------------------------------

    REPORT::rpt.Message( 1, " computation for %5d. time step; %s: %s\n%s %d\n\n", iTM,
                            "                     time", timeint.actualTime.Get(),
                            " number of boundary condition set is", bcSetNo + 1 );


    if( bcSetNo != prevBcSetNo )
    {
      prevBcSetNo = bcSetNo;

      // set specified boundary conditions -----------------------------------------------

      timeint.bconSet[bcSetNo].InitBcon( this, &timeint.actualTime );

      M2D->SetLocation();
      M2D->SetNormal();
      M2D->SetRotation();
    }

    // -----------------------------------------------------------------------------------
    // reset Sums in Statistics
    if( timeint.reset_statist )
    {
      if( timeint.reset_statist[iTM - 1] )  statist->Reset( M2D );
    }

    // -----------------------------------------------------------------------------------
    // solve equations

    for( ;; )
    {
      switch( theCycle )
      {
        // -------------------------------------------------------------------------------
        // Navier-Stokes cycle

        case kUVSCyc:                                                   // 2D: UVS coupled
          if( isFS(actualTurb, BCONSET::kVtAnisotrop) )
          {
            eqs_uvs2d_ai.Execute( this, actualStat );
          }
          else
          {
            eqs_uvs2d.Execute( this, actualStat );
          }
          break;


        case kUVS_TMCyc:                                                // 2D: UVS coupled
          if( isFS(actualTurb, BCONSET::kVtAnisotrop) )
          {
            eqs_uvs2d_tmai.Execute( this, actualStat );
          }
          else
          {
            eqs_uvs2d_tm.Execute( this, actualStat );
          }
          break;


        case kUVS_LVCyc:                           // 2D: UVS coupled, linear UV, const. S
          eqs_uvs2d_lv.Execute( this, actualStat );
          break;

/*
        case kUVS_LVXCyc:
          eqs_uvs2d_lvx.Execute( this );
          break;
*/

        case kDispCurv2D:
          eqs_disp.Execute( this );
          break;

        // -------------------------------------------------------------------------------
        // algebraic eddy viscosity model cycles

        case kKDInitCyc:
          // print information on actual iteration
          PrintTheCycle( 1 );
          REPORT::rpt.PrintTime( 1 );

          // compute friction coefficients
          M2D->DoFriction( this );

          // initialize Reynolds stresses and eddy viscosity
          R2D->Turbulence( this );

          // initialize K and D with algebraic model
          R2D->InitKD( this );
          M2D->SetBoundKD( this );
//        R2D->SmoothKD( smoothPassesKD );

          // initialize Reynolds stresses and eddy viscosity (once again)
          R2D->Turbulence( this );
          break;


        // -------------------------------------------------------------------------------
        // one equation turbulence model cycle

        case kKLCyc:
#         ifdef _MPI_
          if( subdom.npr > 1 )
            REPORT::rpt.Error( kParameterFault, "currently no cycle %d in MPI-Version", theCycle );
#         endif
          eqs_kl2d.Execute( this, actualStat );
//        R2D->smoothKD( smoothPassesKD );
          break;


        // -------------------------------------------------------------------------------
        // k-epsilon cycle (two equation turbulence cycle)

        case kKDCyc:
          eqs_kd2d.Execute( this, actualStat, 0 );
          break;

        case kKD_LCyc:
          eqs_kd2d.Execute( this, actualStat, 1 );
          break;

        case kKD_QCyc:
          eqs_kd2d.Execute( this, actualStat, 2 );
          break;

        // -------------------------------------------------------------------------------

        case kKCyc:
          eqs_k2d.Execute( this, actualStat, 0 );
          break;

        // -------------------------------------------------------------------------------

        case kDCyc:
          eqs_d2d.Execute( this, actualStat, 0 );
          break;

        // -------------------------------------------------------------------------------
        // suspended and bed load cycles

        // sediment transport cycle ------------------------------------------------------
        case kQbCyc:
          // initialize sediment parameters as qbe, Ls, sx, sy, ... ----------------------
          sed.Initialize( this );

          eqs_bl2d.Execute( this, EQS_BL2D::kQuadratic );

          // detach memory for parameters ------------------------------------------------
          sed.Detach();
          break;

        // suspended load ----------------------------------------------------------------
        case kSLCyc:
          // initialize sediment parameters as qbe, Ls, sx, sy, ... ----------------------
          sed.Initialize( this );

          eqs_sl2d.Execute( this );
          eqs_dz.Execute( this, EQS_BL2D::kQuadratic );

          // detach memory for parameters ------------------------------------------------
          sed.Detach();
          break;

        // bed load cycle ----------------------------------------------------------------
        case kBLCyc:
          // initialize sediment parameters as qbe, Ls, sx, sy, ... ----------------------
          sed.Initialize( this );

          eqs_bl2d.Execute( this, EQS_BL2D::kQuadratic );
          eqs_dz.Execute( this, EQS_BL2D::kQuadratic );

          // detach memory for parameters ------------------------------------------------
          sed.Detach();
          break;


        // bed and suspended load cycle --------------------------------------------------
        case kBSLCyc:
          // initialize sediment parameters as qbe, Ls, sx, sy, ... ----------------------
          sed.Initialize( this );

          eqs_bl2d.Execute( this, EQS_BL2D::kQuadratic );
          eqs_dz.Execute( this, EQS_BL2D::kQuadratic );
          eqs_sl2d.Execute( this );

          // detach memory for parameters ------------------------------------------------
          sed.Detach();
          break;


        // difference in bed load used for bottom evolution ------------------------------
        case kDiffBLCyc:
          // initialize sediment parameters as qbe, Ls, sx, sy, ... ----------------------
          sed.Initialize( this );

          eqs_bl2d.Execute( this, EQS_BL2D::kQuadratic );
          eqs_dz.QbDiff( this, EQS_BL2D::kQuadratic );

          // detach memory for parameters ------------------------------------------------
          sed.Detach();
          break;


        // -------------------------------------------------------------------------------
        // divergence free flow field cycle

        case kDivCyc:
          eqs_ppe2d.Execute( this );
          break;
/*
        case kDivCyc_LV:
          // print information on actual iteration
          PrintTheCycle( 1 );
          printTime( 1 );

          eqs_ppe2d_lv.Execute( this, false );
          break;

        case kDivCyc_LV_LM:
          // print information on actual iteration
          PrintTheCycle( 1 );
          printTime( 1 );

          eqs_ppe2d_lv.Execute( this, true );
          break;
*/

        // -------------------------------------------------------------------------------
        // other useful cycles

        case kDryRewet:
          // print information on actual iteration
          PrintTheCycle( 1 );
          REPORT::rpt.PrintTime( 1 );

          // dry and rewet algorithm
          M2D->DoDryRewet( this );
          break;


        case kReOrderCyc:                                // Reorder Elements
#         ifdef _MPI_
          if( subdom.npr > 1 )
            REPORT::rpt.Error( kParameterFault, "currently no cycle %d in MPI-Version", theCycle );
#         endif
          // print information on actual iteration
          PrintTheCycle( 1 );
          REPORT::rpt.PrintTime( 1 );

          R2D->Connection( 0l );

          // reorder Elements
          M2D->list = M2D->ReorderElem( nSection, section );

          M2D->Initialize();
          break;


        case kSurfaceCyc:
          // print information on actual iteration
          PrintTheCycle( 1 );
          REPORT::rpt.PrintTime( 1 );

          R2D->InitS( nSection, section );

          // compute friction coefficients
          M2D->DoFriction( this );

          // initialize Reynolds stresses and eddy viscosity
          R2D->Turbulence( this );
          break;

        case kSurfaceToVol:
          // print information on actual iteration
          PrintTheCycle( 1 );
          REPORT::rpt.PrintTime( 1 );

          {
            for( int e=0; e<R2D->Getne(); e++ )
            {
              ELEM* el = R2D->Getelem(e);

              el->P = 0.0;

              if( !isFS(el->flag, ELEM::kDry) )
              {
                int ncn = el->Getncn();

                for( int i=0; i<ncn; i++ )
                {
                  el->U += el->nd[i]->v.U;
                  el->V += el->nd[i]->v.V;
                  el->P += el->nd[i]->v.S;
                }

                el->P /= ncn;
                el->U /= ncn;
                el->V /= ncn;
              }
            }
          }
          break;

        case kOutputCyc:
          // print information on actual iteration
          PrintTheCycle( 1 );
          REPORT::rpt.PrintTime( 1 );

          // compute friction coefficients
          M2D->DoFriction( this );

          // initialize Reynolds stresses and eddy viscosity
          R2D->Turbulence( this );

          // write output files
          R2D->OutputData( this, iTM, timeint.actualTime.Get() );

          M2D->Output( this, iTM );
          M2D->DetachOutput( this );

          if( *name.geometryFile )
          {
            R2D->OutputGeom( this, iTM );
          }
          break;


        default:
          // compute friction coefficients
          M2D->DoFriction( this );

          // initialize Reynolds stresses and eddy viscosity
          R2D->Turbulence( this );
          break;
      }

      MEMORY::memo.PrintInfo();

      if( errLevel & kErr_interrupt )  break;

      // determine the next iteration cycle ----------------------------------------------
      theCycle = NextCycle( &timeint.bconSet[bcSetNo], false );

      if( !theCycle )
      {
        // proceed to next time step -----------------------------------------------------
        for( int i=0; i<R2D->Getnp(); i++ )
        {
          NODE* nd = R2D->Getnode(i);

          if( isFS(nd->flag, NODE::kDry) )  SF( nd->flag, NODE::kDryPrev );
          else                              CF( nd->flag, NODE::kDryPrev );

          if( isFS(nd->flag, NODE::kMarsh) )  SF( nd->flag, NODE::kMarshPrev );
          else                                CF( nd->flag, NODE::kMarshPrev );

          nd->vo.U = nd->v.U;
          nd->vo.V = nd->v.V;
          nd->vo.S = nd->v.S;

          nd->vo.K = nd->v.K;
          nd->vo.D = nd->v.D;

          nd->vo.C = nd->v.C;

          nd->vo.dUdt = nd->v.dUdt;
          nd->vo.dVdt = nd->v.dVdt;
          nd->vo.dSdt = nd->v.dSdt;
        }

        // -------------------------------------------------------------------------------
        timeint.actualTime += timeint.incTime;

        if( timeint.actualTime < timeint.nextTime )
        {
          REPORT::rpt.Message( 1, "\n actual time at end of time step: %s\n",
                                  timeint.actualTime.Get() );

          timeint.incTime = timeint.nextTime - timeint.actualTime;

          theCycle = NextCycle( &timeint.bconSet[bcSetNo], true );
        }
        else
        {
          timeint.prevTime = timeint.nextTime;
          break;
        }
      }
    }

    // do statistics ---------------------------------------------------------------------
    if( statistics )  statist->Sum( M2D );


    // write results file ----------------------------------------------------------------
    if( timeint.result[iTM - 1]  ||  (errLevel & kErr_interrupt) )
    {
      // create directories when indicated
      char pathname[500];
      ReplaceMacro( outputPath, pathname, iTM, subdom.pid+1 );

      // search for last occurence of "/" in pathname
      char *s = pathname + strlen(pathname)-1;
      while( *s != *pathname )
      {
        if( *s == '/' )
        {
          *s = '\0';
          break;
        }
        s--;
      }

#     if defined(LINUX)
      mkdir( pathname, 0777 );

#     elif defined(MS_WIN)
      mkdir( pathname );

#     else
      REPORT::rpt.Warning( kMethodFault, "No directories for output created! \"mkdir()\" not supported" );
#     endif

      // write output to files
      R2D->OutputData( this, iTM, timeint.actualTime.Get() );

      M2D->Output( this, iTM );

      if( statistics )
      {
        char fileName[80];

        if( name.rgAvsStatistFile[0] )
          sprintf( fileName, name.rgAvsStatistFile, iTM );
        else
          fileName[0] = '\0';

        statist->Write( M2D, release, name.outputStatistFile, fileName, iTM, this );
      }

      if( name.geometryFile[0] )
      {
        R2D->OutputGeom( this, iTM );
      }
    }

    // write time series -----------------------------------------------------------------
    for( int its=0; its<ntmser; its++ )
    {
      if( tmser[its].ntm == 0  ||  (tmser[its].ntm > iTM && tmser[its].tmlist[iTM]) )
      {
        M2D->OutputSeries( this, iTM, &tmser[its] );
      }
    }

    M2D->DetachOutput( this );

    // -----------------------------------------------------------------------------------
    REPORT::rpt.Setlevel( report_level );

    // -----------------------------------------------------------------------------------
    if( errLevel & kErr_interrupt )  break;
  }


  REPORT::rpt.Line2( 1 );

  if( errLevel & kErr_interrupt )
    REPORT::rpt.Message( 0, "= program execution interrupted due to errors\n" );

  else if( errLevel & kErr_some_errors )
    REPORT::rpt.Message( 0, "= normal program termination with some errors\n" );

  else
    REPORT::rpt.Message( 0, "= normal program termination\n" );

  if( errLevel & kErr_no_conv_cg )
    REPORT::rpt.Message( 0, "= cg-solver did not converge in each cycle\n" );

  if( errLevel & kErr_no_conv_nr )
    REPORT::rpt.Message( 0, "= Newton-Raphson did not converge in each time step\n" );

  REPORT::rpt.Line2( 0 );
  REPORT::rpt.PrintTime( 0 );


  // close report file -------------------------------------------------------------------
  REPORT::rpt.Close();
}
