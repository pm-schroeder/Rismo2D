// ======================================================================================
//
// Copyright (C) 1992-2012  by  P.M. SCHROEDER
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
#include "Bcon.h"
#include "Report.h"
#include "Project.h"

#include "Timeint.h"


TIMEINT::TIMEINT()
{
  setsOfBcon = 0;

  firstTimeStep = 1;
  lastTimeStep  = 1;
  bcLoop        = 1;

  thetaFlow = 0.5;
  thetaTurb = 0.5;
  thetaSedi = 0.5;

  result = NULL;
  reset_statist = NULL;

  actualBcSet = NULL;

  // -------------------------------------------------------------------------------------

  static DATKEY dk[] =
  {
    kTM_STEPS,           "TM_STEPS",            //  1
    kTM_INTERVAL,        "TM_INTERVAL",         //  2
    kTM_WEIGHT,          "TM_WEIGHT",           //  3
    kTM_OUTPUT,          "TM_OUTPUT",           //  4
    kTM_SETTIME,         "TM_SETTIME",          //  5

    kTM_STEP_NO,         "TM_STEP_NO",          //  6
    kTM_CYCLE,           "TM_CYCLE",            //  7
    kTM_STATIONARY,      "TM_STATIONARY",       //  8
    kTM_TURBULENCE,      "TM_TURBULENCE",       //  9
    kTM_DISPERSION,      "TM_DISPERSION",       // 10
    kTM_MAXITER,         "TM_MAXITER",          // 11
    kTM_SOLVER,          "TM_SOLVER",           // 12

    kTM_BOUND_LINE,      "TM_BOUND_LINE",       // 13
    kTM_BOUND_NODE,      "TM_BOUND_NODE",       // 14

    kTM_LINE,            "TM_LINE",             // 15
    kTM_NODE,            "TM_NODE",             // 16

    kTM_PERIODIC_LINE,   "TM_PERIODIC_LINE",    // 17
    kTM_PERIODIC_NODE,   "TM_PERIODIC_NODE",    // 18

    kTM_RESET_STATIST,   "TM_RESET_STATIST"     // 19
  };

  datkey = dk;
  nkey   = 19;

  set = false;
  startTime.Set( "0" );

  // -------------------------------------------------------------------------------------

  static RELOC rl[] =
  {
    RELOC::kInflow,    BCON::kInlet,       //   1: inlet boundary
    RELOC::kOutflow,   BCON::kOutlet,      //   2: outlet: water elevation specified
    RELOC::kOpenflow,  BCON::kOpenBnd,     //   3: open boundary; no boundary conditions
    RELOC::kSlip,      BCON::kSlip,        //   4: slip velocity; flow direction
    RELOC::kFixVelo,   BCON::kSetUV,       //   5: fixed velocity (U,V)
    RELOC::kFixS,      BCON::kSetS,        //   6: fixed value for flow depth h
    RELOC::kFixKD,     BCON::kSetKD,       //   7: fixed values for K and D
    RELOC::kFlowC,     BCON::kRateC,       //   8: sediment transport rate [m3/s]
    RELOC::kQInflow,   BCON::kQInlet,      //   9: discharge along control line (inlet)
    RELOC::kQTInflow,  BCON::kQTInlet,     //  10: Q(t)-relation for inlet
    RELOC::kSQOutflow, BCON::kSQOutlet,    //  11: H(Q)-relation for outlet
    RELOC::kSTOutflow, BCON::kSTOutlet,    //  12: S(t)-relation for outlet
    RELOC::kFixC,      BCON::kSetC,        //  13: fixed values for concentration [kg/m3]
    RELOC::kQSource,   BCON::kSource       //  14: diffuse sink/source Q
  };

  reloc  = rl;
  nreloc = 14;
}


TIMEINT::~TIMEINT()
{
  if( result )  delete[] result;
  if( reset_statist )  delete[] reset_statist;
}


void TIMEINT::Input( char* fileName )
{
  // open input file for boundary conditions ---------------------------------------------

  ASCIIFILE* file = new ASCIIFILE( fileName, "r" );

  if( !file  ||  !file->getid() )
    REPORT::rpt.Error( kOpenFileFault, "%s %s (TIMEINT::Input #1)",
                       "can not open time step file", fileName );


  // check version of input file ---------------------------------------------------------

  release = 0;

  char* textLine = file->next();
  sscanf( textLine, " $RISMO2D %d", &release );
  REPORT::rpt.Screen( 2, "\n ### release of time step file:          %d ###\n\n", release );

       if( release >= 40100 )  Input_40100( file );
  else if( release >= 30900 )  Input_30900( file );
  else                         Input_00000( file );

  delete file;
}


void TIMEINT::Input_00000( ASCIIFILE* file )
{
  char* textLine;
  char  text[400];

  file->rewind();


  // -------------------------------------------------------------------------------------
  // read parameters for time integration

  double deltaTm, relaxTmFlow, dummy, relaxTmTurb;

  textLine = file->nextLine();
  sscanf( textLine, " %d %d %d %d %lf %lf %lf %lf",
          &setsOfBcon,
          &firstTimeStep,
          &lastTimeStep,
          &bcLoop,
          &deltaTm,
          &relaxTmFlow,
          &dummy,
          &relaxTmTurb );

  if( relaxTmFlow <= 0.0 || relaxTmFlow > deltaTm )  relaxTmFlow = deltaTm;
  if( relaxTmTurb <= 0.0 || relaxTmTurb > deltaTm )  relaxTmTurb = deltaTm;

  deltaTime.Setsec( deltaTm );
  relaxTimeFlow.Setsec( relaxTmFlow );
  relaxTimeTurb.Setsec( relaxTmTurb );

  REPORT::rpt.Output( "\n", 2 );
  REPORT::rpt.OutputLine2( 2 );

  sprintf( text, "\n  %30s  %-2d\n  %30s  %-4d\n  %30s  %-4d\n",
          "sets of boundary conditions:", setsOfBcon,
          "first time step:",             firstTimeStep,
          "last time step:",              lastTimeStep );
  REPORT::rpt.Output( text, 2 );

  sprintf( text, "  %30s  %-4d\n  %30s  %-s\n",
          "repetition from:",             bcLoop,
          "time interval:",               deltaTime.Get() );
  REPORT::rpt.Output( text, 2 );

  sprintf( text, "  %30s  %-s\n  %30s  %-s\n",
          "relaxed time for flow:",       relaxTimeFlow.Get(),
          "relaxed time for turbulence:", relaxTimeTurb.Get() );
  REPORT::rpt.Output( text, 2 );


  if( firstTimeStep <= 0  ||  lastTimeStep <= 0 )
    REPORT::rpt.Error( "numbers for time steps must be greater than 0 (TIMEINT::Input #2)" );


  REPORT::rpt.Output( "\n", 2 );
  REPORT::rpt.OutputLine2( 2 );

  textLine = file->nextLine();
  sscanf( textLine, " %lf", &thetaFlow );
  thetaTurb = thetaFlow;
  thetaSedi = thetaFlow;

  sprintf( text, "\n  %30s  %-4.2lf\n",
           "time weighting THETA          :", thetaFlow );
  REPORT::rpt.Output( text, 2 );


  // read time step numbers for output ---------------------------------------------------

  result = new int [lastTimeStep];
  if( !result ) REPORT::rpt.Error( "can not allocate memory - TIMEINT::Input #3" );

  for( int i=0; i<lastTimeStep; i++ )  result[i] = false;


  int nLines;
  textLine = file->nextLine();
  sscanf( textLine, " %d", &(nLines) );


  for( int i=0; i<nLines; i++ )
  {
    int start, stop, step;

    textLine = file->nextLine();

    int j = sscanf( textLine, " %d to %d step %d", &start, &stop, &step );

    if( j == 3 )
    {
      for( int k=start; k<=stop && k<=lastTimeStep; k+=step )
      {
        result[k-1] = true;
      }
    }
    else
    {
      int stepno[10];

      j = sscanf( textLine, " %d %d %d %d %d %d %d %d %d %d",
                            &(stepno[0]), &(stepno[1]),
                            &(stepno[2]), &(stepno[3]),
                            &(stepno[4]), &(stepno[5]),
                            &(stepno[6]), &(stepno[7]),
                            &(stepno[8]), &(stepno[9]) );

      for( int k=0; k<10; k++ )
      {
        if( stepno[k] <= lastTimeStep )  result[stepno[k]-1] = true;
      }
    }
  }


  REPORT::rpt.Output( "\n",2 );
  REPORT::rpt.OutputLine2( 2 );

  sprintf( text, "\n output will be written for time steps ...\n" );
  REPORT::rpt.Output( text, 2 );

  for( int i=0, j=0; i<lastTimeStep; i++ )
  {
    if( result[i] )
    {
      sprintf( text, "  %4d", i + 1 );
      REPORT::rpt.Output( text, 2 );
      j++;

      if ( j % 10 == 0 )  REPORT::rpt.Output( "\n", 2 );
    }
  }


  // allocate memory for different sets of boundary conditions ---------------------------

  bconSet = new BCONSET [setsOfBcon];

  if( !bconSet )
    REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input #4)" );


  for( int i=0; i<setsOfBcon; i++ )
  {
    int ncyc;


    // read specifications for ith time step ---------------------------------------------
    // bc's at first time step are for stationary solution

    textLine = file->nextLine();
    sscanf( textLine, " %d %d %d", &bconSet[i].noTimeStep,
                                   &bconSet[i].nofLines,
                                   &bconSet[i].nofNodes );

    endLoop = bconSet[i].noTimeStep;

    if( i == 0 )
    {
      if( bconSet[i].noTimeStep != 1 )
        REPORT::rpt.Error( "initial time step number must be 1 (TIMEINT::Input #5)" );
    }

    REPORT::rpt.Output( "\n\n", 2 );
    REPORT::rpt.OutputLine1( 2 );

    sprintf( text, "  %s %3d. %s\n\n",
                   "specifications for",
                   bconSet[i].noTimeStep, "time step" );
    REPORT::rpt.Output( text, 2 );


    // read cycle sequence ---------------------------------------------------------------

    for( int j=0; j<kMaxCycles; j++ )  bconSet[i].cycle[j] = 0;

    textLine = file->nextLine();
    ncyc = sscanf( textLine, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                             &(bconSet[i].cycle[ 0]), &(bconSet[i].cycle[ 1]),
                             &(bconSet[i].cycle[ 2]), &(bconSet[i].cycle[ 3]),
                             &(bconSet[i].cycle[ 4]), &(bconSet[i].cycle[ 5]),
                             &(bconSet[i].cycle[ 6]), &(bconSet[i].cycle[ 7]),
                             &(bconSet[i].cycle[ 8]), &(bconSet[i].cycle[ 9]),
                             &(bconSet[i].cycle[10]), &(bconSet[i].cycle[11]),
                             &(bconSet[i].cycle[12]), &(bconSet[i].cycle[13]),
                             &(bconSet[i].cycle[14]), &(bconSet[i].cycle[15]),
                             &(bconSet[i].cycle[16]), &(bconSet[i].cycle[17]),
                             &(bconSet[i].cycle[18]), &(bconSet[i].cycle[19]),
                             &(bconSet[i].cycle[20]), &(bconSet[i].cycle[21]),
                             &(bconSet[i].cycle[22]), &(bconSet[i].cycle[23]),
                             &(bconSet[i].cycle[24]), &(bconSet[i].cycle[25]),
                             &(bconSet[i].cycle[26]), &(bconSet[i].cycle[27]),
                             &(bconSet[i].cycle[28]), &(bconSet[i].cycle[29]),
                             &(bconSet[i].cycle[30]), &(bconSet[i].cycle[31]),
                             &(bconSet[i].cycle[32]), &(bconSet[i].cycle[33]),
                             &(bconSet[i].cycle[34]), &(bconSet[i].cycle[35]) );

    REPORT::rpt.Output( "\n  sequence of cycles ...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 6 )  REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "  %8d", bconSet[i].cycle[j] );
      REPORT::rpt.Output( text, 2 );
    }


    // read turbulence model according to cycle switch -----------------------------------

    for( int j=0; j<kMaxCycles; j++ )  bconSet[i].turb[j] = 0;

    textLine = file->nextLine();
    sscanf( textLine, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                      &(bconSet[i].turb[ 0]), &(bconSet[i].turb[ 1]),
                      &(bconSet[i].turb[ 2]), &(bconSet[i].turb[ 3]),
                      &(bconSet[i].turb[ 4]), &(bconSet[i].turb[ 5]),
                      &(bconSet[i].turb[ 6]), &(bconSet[i].turb[ 7]),
                      &(bconSet[i].turb[ 8]), &(bconSet[i].turb[ 9]),
                      &(bconSet[i].turb[10]), &(bconSet[i].turb[11]),
                      &(bconSet[i].turb[12]), &(bconSet[i].turb[13]),
                      &(bconSet[i].turb[14]), &(bconSet[i].turb[15]),
                      &(bconSet[i].turb[16]), &(bconSet[i].turb[17]),
                      &(bconSet[i].turb[18]), &(bconSet[i].turb[19]),
                      &(bconSet[i].turb[20]), &(bconSet[i].turb[21]),
                      &(bconSet[i].turb[22]), &(bconSet[i].turb[23]),
                      &(bconSet[i].turb[24]), &(bconSet[i].turb[25]),
                      &(bconSet[i].turb[26]), &(bconSet[i].turb[27]),
                      &(bconSet[i].turb[28]), &(bconSet[i].turb[29]),
                      &(bconSet[i].turb[30]), &(bconSet[i].turb[31]),
                      &(bconSet[i].turb[32]), &(bconSet[i].turb[33]),
                      &(bconSet[i].turb[34]), &(bconSet[i].turb[35]) );

    REPORT::rpt.Output( "\n\n  turbulence model ...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 6 ) REPORT::rpt.Output( "\n", 2 );

      int turb = bconSet[i].turb[j];

      sprintf( text, "  %8d", turb );
      REPORT::rpt.Output( text, 2 );

      bconSet[i].turb[j] = 0;

      switch( turb % 10 )
      {
        case 1:   bconSet[i].turb[j] |= BCONSET::kVtConstant;   break;
        case 2:   bconSet[i].turb[j] |= BCONSET::kVtPrandtlKol; break;
        case 3:   bconSet[i].turb[j] |= BCONSET::kVtPrandtlKol
                                      | BCONSET::kVtAnisotrop;  break;
      }

      switch( (turb / 10) % 10 )
      {
        case 1:   bconSet[i].turb[j] |= BCONSET::kVtMin;        break;
      }

      switch( (turb / 100) % 10 )
      {
        case 1:   bconSet[i].turb[j] |= BCONSET::kVtIterat;     break;
      }
    }


    // read maximum number of iterations per cycle ---------------------------------------

    for( int j=0; j<kMaxCycles; j++ )  bconSet[i].cycit[j] = 1;

    textLine = file->nextLine();
    sscanf( textLine, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                      &(bconSet[i].cycit[ 0]), &(bconSet[i].cycit[ 1]),
                      &(bconSet[i].cycit[ 2]), &(bconSet[i].cycit[ 3]),
                      &(bconSet[i].cycit[ 4]), &(bconSet[i].cycit[ 5]),
                      &(bconSet[i].cycit[ 6]), &(bconSet[i].cycit[ 7]),
                      &(bconSet[i].cycit[ 8]), &(bconSet[i].cycit[ 9]),
                      &(bconSet[i].cycit[10]), &(bconSet[i].cycit[11]),
                      &(bconSet[i].cycit[12]), &(bconSet[i].cycit[13]),
                      &(bconSet[i].cycit[14]), &(bconSet[i].cycit[15]),
                      &(bconSet[i].cycit[16]), &(bconSet[i].cycit[17]),
                      &(bconSet[i].cycit[18]), &(bconSet[i].cycit[19]),
                      &(bconSet[i].cycit[20]), &(bconSet[i].cycit[21]),
                      &(bconSet[i].cycit[22]), &(bconSet[i].cycit[23]),
                      &(bconSet[i].cycit[24]), &(bconSet[i].cycit[25]),
                      &(bconSet[i].cycit[26]), &(bconSet[i].cycit[27]),
                      &(bconSet[i].cycit[28]), &(bconSet[i].cycit[29]),
                      &(bconSet[i].cycit[30]), &(bconSet[i].cycit[31]),
                      &(bconSet[i].cycit[32]), &(bconSet[i].cycit[33]),
                      &(bconSet[i].cycit[34]), &(bconSet[i].cycit[35]) );

    REPORT::rpt.Output( "\n\n  iterations per cycle...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 6 )  REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "  %8d", bconSet[i].cycit[j] );
      REPORT::rpt.Output( text, 2 );
    }


    // read maximum number of iterations per cycle ---------------------------------------

    for( int j=0; j<kMaxCycles; j++ )  bconSet[i].solver[j] = 1;

    textLine = file->nextLine();
    sscanf( textLine, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                      &(bconSet[i].solver[ 0]), &(bconSet[i].solver[ 1]),
                      &(bconSet[i].solver[ 2]), &(bconSet[i].solver[ 3]),
                      &(bconSet[i].solver[ 4]), &(bconSet[i].solver[ 5]),
                      &(bconSet[i].solver[ 6]), &(bconSet[i].solver[ 7]),
                      &(bconSet[i].solver[ 8]), &(bconSet[i].solver[ 9]),
                      &(bconSet[i].solver[10]), &(bconSet[i].solver[11]),
                      &(bconSet[i].solver[12]), &(bconSet[i].solver[13]),
                      &(bconSet[i].solver[14]), &(bconSet[i].solver[15]),
                      &(bconSet[i].solver[16]), &(bconSet[i].solver[17]),
                      &(bconSet[i].solver[18]), &(bconSet[i].solver[19]),
                      &(bconSet[i].solver[20]), &(bconSet[i].solver[21]),
                      &(bconSet[i].solver[22]), &(bconSet[i].solver[23]),
                      &(bconSet[i].solver[24]), &(bconSet[i].solver[25]),
                      &(bconSet[i].solver[26]), &(bconSet[i].solver[27]),
                      &(bconSet[i].solver[28]), &(bconSet[i].solver[29]),
                      &(bconSet[i].solver[30]), &(bconSet[i].solver[31]),
                      &(bconSet[i].solver[32]), &(bconSet[i].solver[33]),
                      &(bconSet[i].solver[34]), &(bconSet[i].solver[35]) );

    REPORT::rpt.Output( "\n\n  solver for cycle...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 6 ) REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "  %8d", bconSet[i].solver[j] );
      REPORT::rpt.Output( text, 2 );
    }

    REPORT::rpt.Output( "\n\n", 2 );
    REPORT::rpt.OutputLine1( 2 );


    // read specifications for lines -----------------------------------------------------

    if( bconSet[i].nofLines > 0 )
    {
      bconSet[i].bcLine = new BCONLINE [bconSet[i].nofLines];

      if( !bconSet[i].bcLine )
        REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input #6)" );

      BCONLINE* line = bconSet[i].bcLine;


      for( int j=0; j<bconSet[i].nofLines; j++ )
      {
        long fix;

        textLine = file->nextLine();
        sscanf( textLine, " %d %ld %d", &line[j].no,
                                        &fix,
                                        &line[j].ndat );

        line[j].kind = reloc->RelocBcon( fix, nreloc, reloc );

        line[j].x  = new double [line[j].ndat];
        line[j].y  = new double [line[j].ndat];
        line[j].U  = new double [line[j].ndat];
        line[j].V  = new double [line[j].ndat];
        line[j].S  = new double [line[j].ndat];
        line[j].K  = new double [line[j].ndat];
        line[j].D  = new double [line[j].ndat];
        line[j].C  = new double [line[j].ndat];
        line[j].Qb = new double [line[j].ndat];

        if(     !line[j].x  ||  !line[j].y
            ||  !line[j].U  ||  !line[j].V
            ||  !line[j].S  ||  !line[j].K  ||  !line[j].D
            ||  !line[j].C  ||  !line[j].Qb )
            REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input #7)" );

        for( int k=0; k<line[j].ndat; k++ )
        {
          textLine = file->nextLine();
          sscanf( textLine, " %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                            &line[j].x[k],
                            &line[j].y[k],
                            &line[j].U[k],
                            &line[j].V[k],
                            &line[j].S[k],
                            &line[j].K[k],
                            &line[j].D[k],
                            &line[j].C[k],
                            &line[j].Qb[k] );
        }
      }
    }


    // read specifications for nodes -----------------------------------------------------

    if( bconSet[i].nofNodes > 0 )
    {
      bconSet[i].bcNode = new BCON [ bconSet[i].nofNodes ];
      BCVAL* bcval      = new BCVAL[ bconSet[i].nofNodes ];

      if( !bconSet[i].bcNode || !bcval )
        REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input #8)" );


      for( int j=0; j<bconSet[i].nofNodes; j++ )
      {
        int    no;
        long   type;
        double USpec, VSpec, SSpec, KSpec, DSpec, CSpec, qbSpec;

        textLine = file->nextLine();
        sscanf( textLine, " %d %ld %lf %lf %lf %lf %lf %lf %lf",
                          &no, &type, &USpec, &VSpec, &SSpec,
                          &KSpec, &DSpec, &CSpec, &qbSpec );


        bconSet[i].bcNode[j].no         = no - 1;
        bconSet[i].bcNode[j].kind       = reloc->RelocBcon( type, nreloc, reloc );

        bconSet[i].bcNode[j].val        = &bcval[j];

        bconSet[i].bcNode[j].val->U     = USpec;
        bconSet[i].bcNode[j].val->V     = VSpec;
        bconSet[i].bcNode[j].val->S     = SSpec;
        bconSet[i].bcNode[j].val->K     = KSpec;
        bconSet[i].bcNode[j].val->D     = DSpec;
        bconSet[i].bcNode[j].val->C[0]  = CSpec;
        bconSet[i].bcNode[j].val->Qb[0] = qbSpec;
      }
    }


    // report boundary conditions --------------------------------------------------------

    REPORT::rpt.Output( "\n", 2 );
    REPORT::rpt.OutputLine1( 2 );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////

void TIMEINT::Input_30900( ASCIIFILE* file )
{
  // -------------------------------------------------------------------------------------

  char* textLine;
  char  text[400];

  // -------------------------------------------------------------------------------------
  // the file is read twice

  // first reading of file:
  //   1. determine number of boundary sets
  //   2. determine number of last time step
  //   3. determine maximum number of boundary lines and nodes
  //   4. determine number of nodes & lines with periodic boundary condition

  int cntOutput = 0;
  int cntbcNode = 0;
  int maxbcNode = 0;

  int nbcLine = 100;

  int* bcLine    = new int[nbcLine];
  int* cntbcLine = new int[nbcLine];
  int* maxbcLine = new int[nbcLine];
  if( !bcLine || !cntbcLine || !maxbcLine )
    REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #1)" );

  memset( bcLine,    0, nbcLine );
  memset( cntbcLine, 0, nbcLine );
  memset( maxbcLine, 0, nbcLine );

  nPeriodicNode = 0;
  nPeriodicLine = 0;

  while( !feof(file->getid()) )
  {
    if( !(textLine = file->nextLine()) )  break;

    char key[50];
    int  ikey = -99;

    // remove leading blanks
    while( *textLine == ' ' )  textLine++;

    if( textLine[0] == '$' )
    {
      sscanf( textLine, "$%s", key );

      for( int i=0; i<nkey; i++ )
      {
        if( strcmp(key,datkey[i].name) == 0 )
        {
          ikey = datkey[i].id;
          break;
        }
      }
    }

    switch( ikey )
    {
      case kTM_OUTPUT:
        cntOutput++;
        break;

      // read number of first and last time step and bcLoop ------------------------------
      case kTM_STEPS:
        sscanf( textLine, "$TM_STEPS %d %d %d",
                          &firstTimeStep,
                          &lastTimeStep,
                          &bcLoop );

        if( firstTimeStep <= 0  ||  lastTimeStep <= 0 )
          REPORT::rpt.Error( "time step numbers > 0 expected (TIMEINT::Input_30900 #2)" );

        result = new int[lastTimeStep];
        if( !result )  REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #3)" );

        memset( result, 0, lastTimeStep );
        result[lastTimeStep-1] = true;
        break;

      // read start time -----------------------------------------------------------------
      case kTM_SETTIME:
        {
          char tmstr[30];
          sscanf( textLine, "$TM_SETTIME %s", tmstr );
          set = true;
          startTime.Set( tmstr );
        }
        break;

      // read length of time step in seconds ---------------------------------------------
      case kTM_INTERVAL:
        {
          double deltaTm, relaxTmFlow, relaxTmTurb, dummy;

          if( release >= 40070 )
          {
            sscanf( textLine, "$TM_INTERVAL %lf %lf %lf",
                              &deltaTm,
                              &relaxTmFlow,
                              &relaxTmTurb );
          }
          else
          {
            sscanf( textLine, "$TM_INTERVAL %lf %lf %lf %lf",
                              &deltaTm,
                              &relaxTmFlow,
                              &dummy,
                              &relaxTmTurb );
          }

          if( relaxTmFlow <= 0.0 || relaxTmFlow > deltaTm )  relaxTmFlow = deltaTm;
          if( relaxTmTurb <= 0.0 || relaxTmTurb > deltaTm )  relaxTmTurb = deltaTm;

          deltaTime.Setsec( deltaTm );
          relaxTimeFlow.Setsec( relaxTmFlow );
          relaxTimeTurb.Setsec( relaxTmTurb );
        }
        break;

      // read weighting of time step -----------------------------------------------------
      case kTM_WEIGHT:
        sscanf( textLine, "$TM_WEIGHT %lf %lf %lf", &thetaFlow, &thetaTurb, &thetaSedi );
        break;

      // determine number of boundary sets -----------------------------------------------
      case kTM_STEP_NO:
        setsOfBcon++;

        for( int i=0; i<nbcLine; i++ )
        {
          if( bcLine[i] == 0 )  break;
          if( cntbcLine[i] > maxbcLine[i] )  maxbcLine[i] = cntbcLine[i];
        }
        memset( cntbcLine, 0, nbcLine );

        if( cntbcNode > maxbcNode )  maxbcNode = cntbcNode;
        break;

      // determine number of lines with boundary conditions ------------------------------
      case kTM_BOUND_LINE:
        {
          int no;
          sscanf( textLine, "$TM_BOUND_LINE %d", &no );

          if( no > 0 )
          {
            for( int i=0; i<nbcLine; i++ )
            {
              if( bcLine[i] == 0  ||  bcLine[i] == no )
              {
                bcLine[i] = no;
                cntbcLine[i]++;
                no = -1;
                break;
              }
            }

            // increase nbcLine (=100) and the corresponding arrays
            // if the number of boundary lines is larger
            if( no > 0 )
            {
              int* new_bcLine    = new int[2*nbcLine];
              int* new_cntbcLine = new int[2*nbcLine];
              int* new_maxbcLine = new int[2*nbcLine];

              if( !new_bcLine || !new_cntbcLine || !new_maxbcLine )
                REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #4)" );

              memset( new_bcLine,    0, 2*nbcLine );
              memset( new_cntbcLine, 0, 2*nbcLine );
              memset( new_maxbcLine, 0, 2*nbcLine );

              memcpy( new_bcLine,    bcLine,    nbcLine );
              memcpy( new_cntbcLine, cntbcLine, nbcLine );
              memcpy( new_maxbcLine, maxbcLine, nbcLine );

              delete[] bcLine;
              delete[] cntbcLine;
              delete[] maxbcLine;

              bcLine    = new_bcLine;
              cntbcLine = new_cntbcLine;
              maxbcLine = new_maxbcLine;

              bcLine[nbcLine] = no;
              cntbcLine[nbcLine]++;

              nbcLine *= 2;
            }
          }
        }
        break;

      // determine number of nodes with boundary conditions ------------------------------
      case kTM_BOUND_NODE:
        cntbcNode++;
        break;

      // determine number of nodes with periodic boundary condition ----------------------
      case kTM_PERIODIC_NODE:
        nPeriodicNode++;
        break;

      // determine number of lines with periodic boundary condition ----------------------
      case kTM_PERIODIC_LINE:
        nPeriodicLine++;
        break;
    }
  }

  for( int i=0; i<nbcLine; i++ )
  {
    if( bcLine[i] == 0 )  break;
    if( cntbcLine[i] > maxbcLine[i] )  maxbcLine[i] = cntbcLine[i];
  }

  if( cntbcNode > maxbcNode )  maxbcNode = cntbcNode;


  // -------------------------------------------------------------------------------------
  // allocate memory for some arrays

  int set  = -1;
  int ncyc =  0;

  if( setsOfBcon <= 0 )
  {
    endLoop    = 1;
    setsOfBcon = 1;
    set = 0;
  }

  bconSet = new BCONSET[setsOfBcon+1];
  if( !bconSet )
    REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #5)" );

  bconSet[0].noTimeStep = 1;


  // -------------------------------------------------------------------------------------

  for( int i=0; i<nbcLine; i++ )
  {
    if( bcLine[i] == 0 )
    {
      nbcLine = i;
      break;
    }
  }

  for( int i=0; i<setsOfBcon; i++ )
  {
    bconSet[i].bcNode = new BCON [maxbcNode];
    BCVAL* bcval      = new BCVAL[maxbcNode];

    bconSet[i].bcLine = new BCONLINE[nbcLine];

    if( !bconSet[i].bcNode || !bcval || !bconSet[i].bcLine )
      REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #6)" );

    for( int j=0; j<maxbcNode; j++ )  bconSet[i].bcNode[j].val = &bcval[j];

    for( int j=0; j<nbcLine; j++ )
    {
      bconSet[i].bcLine[j].x  = new double [maxbcLine[j]];
      bconSet[i].bcLine[j].y  = new double [maxbcLine[j]];
      bconSet[i].bcLine[j].U  = new double [maxbcLine[j]];
      bconSet[i].bcLine[j].V  = new double [maxbcLine[j]];
      bconSet[i].bcLine[j].S  = new double [maxbcLine[j]];
      bconSet[i].bcLine[j].K  = new double [maxbcLine[j]];
      bconSet[i].bcLine[j].D  = new double [maxbcLine[j]];
      bconSet[i].bcLine[j].C  = new double [maxbcLine[j]];
      bconSet[i].bcLine[j].Qb = new double [maxbcLine[j]];

      if(     !bconSet[i].bcLine[j].x  ||  !bconSet[i].bcLine[j].y
          ||  !bconSet[i].bcLine[j].U  ||  !bconSet[i].bcLine[j].V
          ||  !bconSet[i].bcLine[j].S  ||  !bconSet[i].bcLine[j].K
          ||  !bconSet[i].bcLine[j].D  ||  !bconSet[i].bcLine[j].C
          ||  !bconSet[i].bcLine[j].Qb )
          REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #7)" );
    }
  }

  delete[] bcLine;
  delete[] cntbcLine;
  delete[] maxbcLine;


  // -------------------------------------------------------------------------------------

  if( !result )
  {
    result = new int[lastTimeStep];
    if( !result )  REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #8)" );
  }

  if( !cntOutput )
  {
    for( int i=0; i<lastTimeStep; i++ )  result[i] = true;
  }
  else
  {
    for( int i=0; i<lastTimeStep; i++ )  result[i] = false;
  }


  // -------------------------------------------------------------------------------------

  if( nPeriodicNode > 0 )
  {
    periodicNode[0] = new int[nPeriodicNode];
    periodicNode[1] = new int[nPeriodicNode];
    if( !periodicNode[0] || !periodicNode[1] )
      REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #9)" );
  }

  if( nPeriodicLine > 0 )
  {
    periodicLine[0] = new int[nPeriodicLine];
    periodicLine[1] = new int[nPeriodicLine];
    periodicLine[2] = new int[nPeriodicLine];
    periodicLine[3] = new int[nPeriodicLine];
    if( !periodicLine[0] || !periodicLine[1] || !periodicLine[2] || !periodicLine[3] )
      REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #10)" );
  }

  nPeriodicNode = 0;
  nPeriodicLine = 0;


  // -------------------------------------------------------------------------------------
  // second reading of file
  // -------------------------------------------------------------------------------------

  file->rewind();

  while( !feof(file->getid()) )
  {
    if( !(textLine = file->nextLine()) )  break;

    char key[50];
    int  ikey = -99;

    if( textLine[0] == '$' )
    {
      sscanf( textLine, "$%s", key );

      for( int i=0; i<nkey; i++ )
      {
        if( strcmp(key,datkey[i].name) == 0 )
        {
          ikey = datkey[i].id;
          break;
        }
      }
    }

    switch( ikey )
    {
      // read time step numbers for output -----------------------------------------------
      case kTM_OUTPUT:
        {
          int start, stop, step;
          int n = sscanf( textLine, "$TM_OUTPUT %d to %d step %d", &start, &stop, &step );

          if( n == 3 )
          {
            for( int i=start; i<=stop && i<=lastTimeStep; i+=step )  result[i-1] = true;
          }
          else
          {
            char  list[500];
            char  seps[] = " ,\t\n\r";
            char* token;

            strcpy( list, textLine+10 );
            token = strtok( list, seps );

            while( token != NULL )
            {
              int stepno;
              sscanf( token, "%d", &stepno );
              ncyc++;

              if( stepno <= lastTimeStep )  result[stepno-1] = true;

              token = strtok( NULL, seps );
            }
          }
        }
        break;

      case kTM_STEP_NO:
        {
          set++;
          sscanf( textLine, "$TM_STEP_NO %d", &bconSet[set].noTimeStep );

          endLoop = bconSet[set].noTimeStep;

          if( set == 0 )
          {
            if( bconSet[set].noTimeStep != 1 )
              REPORT::rpt.Error( "initial time step number must be 1 (TIMEINT::Input_30900 #11)" );
          }

          ncyc = 0;
        }
        break;

      case kTM_CYCLE:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          strcpy( list, textLine+9 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].cycle[ncyc] );
            ncyc++;

            token = strtok( NULL, seps );
          }
        }
        break;

      case kTM_TURBULENCE:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+14 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].turb[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].turb[i] = bconSet[set].turb[n-1];
        }
        break;

      case kTM_DISPERSION:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+14 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].disp[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].disp[i] = bconSet[set].disp[n-1];
        }
        break;

      case kTM_MAXITER:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+11 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].cycit[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].cycit[i] = bconSet[set].cycit[n-1];
        }
        break;

      case kTM_SOLVER:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+10 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].solver[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].solver[i] = bconSet[set].solver[n-1];
        }
        break;

      case kTM_BOUND_LINE:
        {
          int  no;
          int  id   = -1;
          long kind =  0;

          sscanf( textLine, "$TM_BOUND_LINE %d %ld", &no, &kind );

          if( no > 0 )
          {
            for( int i=0; i<bconSet[set].nofLines; i++ )
            {
              if( bconSet[set].bcLine[i].no == no )
              {
                id = i;
                break;
              }
            }

            if( id == -1 )
            {
              id = bconSet[set].nofLines;
              bconSet[set].nofLines++;
            }

            bconSet[set].bcLine[id].no   = no;
            bconSet[set].bcLine[id].kind = reloc->RelocBcon( kind, nreloc, reloc );

            int k = bconSet[set].bcLine[id].ndat;

            bconSet[set].bcLine[id].x[k]  = 0.0;
            bconSet[set].bcLine[id].y[k]  = 0.0;
            bconSet[set].bcLine[id].U[k]  = 0.0;
            bconSet[set].bcLine[id].V[k]  = 0.0;
            bconSet[set].bcLine[id].S[k]  = 0.0;
            bconSet[set].bcLine[id].K[k]  = 0.0;
            bconSet[set].bcLine[id].D[k]  = 0.0;
            bconSet[set].bcLine[id].C[k]  = 0.0;
            bconSet[set].bcLine[id].Qb[k] = 0.0;

            sscanf( textLine, "$TM_BOUND_LINE %d %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                              &no, &kind,
                              &bconSet[set].bcLine[id].x[k],
                              &bconSet[set].bcLine[id].y[k],
                              &bconSet[set].bcLine[id].U[k],
                              &bconSet[set].bcLine[id].V[k],
                              &bconSet[set].bcLine[id].S[k],
                              &bconSet[set].bcLine[id].K[k],
                              &bconSet[set].bcLine[id].D[k],
                              &bconSet[set].bcLine[id].C[k],
                              &bconSet[set].bcLine[id].Qb[k] );

            bconSet[set].bcLine[id].ndat++;
          }
        }
        break;

      case kTM_BOUND_NODE:
        {
          int  no;
          long kind = 0;
          double USpec  = 0.0,
                 VSpec  = 0.0,
                 SSpec  = 0.0,
                 KSpec  = 0.0,
                 DSpec  = 0.0,
                 CSpec  = 0.0,
                 qbSpec = 0.0;

          int j = bconSet[set].nofNodes;

          sscanf( textLine, "$TM_BOUND_NODE %d %ld %lf %lf %lf %lf %lf %lf %lf",
                            &no, &kind,
                            &USpec, &VSpec, &SSpec,
                            &KSpec, &DSpec,
                            &CSpec, &qbSpec );

          bconSet[set].bcNode[j].no         = no - 1;
          bconSet[set].bcNode[j].kind       = reloc->RelocBcon( kind, nreloc, reloc );

          bconSet[set].bcNode[j].val->U     = USpec;
          bconSet[set].bcNode[j].val->V     = VSpec;
          bconSet[set].bcNode[j].val->S     = SSpec;
          bconSet[set].bcNode[j].val->K     = KSpec;
          bconSet[set].bcNode[j].val->D     = DSpec;
          bconSet[set].bcNode[j].val->C[0]  = CSpec;
          bconSet[set].bcNode[j].val->Qb[0] = qbSpec;

          bconSet[set].nofNodes++;
        }
        break;


      // read periodic boundary conditions for nodes -------------------------------------
      case kTM_PERIODIC_NODE:
        sscanf( textLine, "$TM_PERIODIC_NODE %d %d",
                &periodicNode[0][nPeriodicNode], &periodicNode[1][nPeriodicNode] );
        nPeriodicNode++;
        break;


      // read periodic boundary conditions for nodes -------------------------------------
      case kTM_PERIODIC_LINE:
        sscanf( textLine, "$TM_PERIODIC_LINE %d %d %d %d",
                &periodicLine[0][nPeriodicLine], &periodicLine[1][nPeriodicLine],
                &periodicLine[2][nPeriodicLine], &periodicLine[3][nPeriodicLine] );
        nPeriodicLine++;
        break;
    }
  }


  // -------------------------------------------------------------------------------------
  // report parameters for time integration

  REPORT::rpt.Output( "\n", 2 );
  REPORT::rpt.OutputLine2( 2 );

  sprintf( text, "\n  %30s  %-2d\n  %30s  %-4d\n  %30s  %-4d\n",
          "sets of boundary conditions:", setsOfBcon,
          "first time step:",             firstTimeStep,
          "last time step:",              lastTimeStep );
  REPORT::rpt.Output( text, 2 );

  sprintf( text, "  %30s  %-4d\n  %30s  %-s\n",
          "repetition from:",             bcLoop,
          "time interval:",               deltaTime.Get() );
  REPORT::rpt.Output( text, 2 );

  sprintf( text, "  %30s  %-s\n  %30s  %-s\n",
          "relaxed time for flow:",       relaxTimeFlow.Get(),
          "relaxed time for turbulence:", relaxTimeTurb.Get() );
  REPORT::rpt.Output( text, 2 );

  REPORT::rpt.Output( "\n", 2 );
  REPORT::rpt.OutputLine1( 2 );

  sprintf( text, "\n  %30s  %-4.2lf\n  %30s  %-4.2lf\n  %30s  %-4.2lf\n",
          "time weighting THETA,  flow:", thetaFlow,
          "                 turbulence:", thetaTurb,
          "                   sediment:", thetaSedi );
  REPORT::rpt.Output( text, 2 );


  // report time step numbers for output -------------------------------------------------

  REPORT::rpt.Output( "\n",2 );
  REPORT::rpt.OutputLine1( 2 );

  sprintf( text, "\n   output will be written for time steps ...\n" );
  REPORT::rpt.Output( text, 2 );

  for( int i=0, j=0; i<lastTimeStep; i++ )
  {
    if( this->result[i] )
    {
      sprintf( text, "   %5d", i + 1 );
      REPORT::rpt.Output( text, 2 );
      j++;

      if ( j % 10 == 0 )  REPORT::rpt.Output( "\n", 2 );
    }
  }


  REPORT::rpt.Output( "\n\n",2 );
  REPORT::rpt.OutputLine2( 2 );


  // report specifications for boundary conditions ---------------------------------------

  for( int i=0; i<setsOfBcon; i++ )
  {
    sprintf( text, "\n   %24s  %-d\n   %24s  %-d\n   %24s  %-d\n",
                   "specifications for set:", bconSet[i].noTimeStep,
                   "       number of lines:", bconSet[i].nofLines,
                   "       number of nodes:", bconSet[i].nofNodes );
    REPORT::rpt.Output( text, 2 );


    // report cycle sequence -------------------------------------------------------------

    REPORT::rpt.Output( "\n   sequence of cycles ...\n\n", 2 );

    for( int j=0; j<kMaxCycles; j++ )
    {
      ncyc = j;
      if( bconSet[i].cycle[j] == 0 )  break;

      if( j == 8 )  REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "   %6d", bconSet[i].cycle[j] );
      REPORT::rpt.Output( text, 2 );
    }


    // report turbulence model according to cycle switch ---------------------------------

    REPORT::rpt.Output( "\n\n   turbulence model ...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 8 ) REPORT::rpt.Output( "\n", 2 );

      int turb = bconSet[i].turb[j];

      sprintf( text, "   %6d", turb );
      REPORT::rpt.Output( text, 2 );

      bconSet[i].turb[j] = 0;

      switch( turb%10 )
      {
        case 1:   bconSet[i].turb[j] |= BCONSET::kVtConstant;                     break;
        case 2:   bconSet[i].turb[j] |= BCONSET::kVtAlgebraic;                    break;
        case 3:   bconSet[i].turb[j] |= BCONSET::kVtAlgShear;                     break;
        case 4:   bconSet[i].turb[j] |= BCONSET::kVtMixingLength;                 break;
        case 5:   bconSet[i].turb[j] |= BCONSET::kVtLES;                          break;
        case 6:   bconSet[i].turb[j] |= BCONSET::kVtPrandtlKol;                   break;
        case 7:   bconSet[i].turb[j] |= BCONSET::kVtLES | BCONSET::kVtAlgebraic;  break;
      }

      if( (turb/10)%10 == 1 )
      {
        bconSet[i].turb[j] |= BCONSET::kVtAnisotrop;
      }

      if( (turb/100)%10 & 1 )
      {
        bconSet[i].turb[j] |= BCONSET::kVtMin;
      }

      if( (turb/100)%10 & 2 )
      {
        bconSet[i].turb[j] |= BCONSET::kVtIterat;
      }
    }


    // report dispersion model according to cycle switch ---------------------------------

    REPORT::rpt.Output( "\n\n   dispersion model ...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 8 )  REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "   %6d", bconSet[i].disp[j] );
      REPORT::rpt.Output( text, 2 );
    }


    // report maximum number of iterations per cycle -------------------------------------

    REPORT::rpt.Output( "\n\n   iterations per cycle...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 8 )  REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "   %6d", bconSet[i].cycit[j] );
      REPORT::rpt.Output( text, 2 );
    }


    // report solver type according to cycle switch --------------------------------------

    REPORT::rpt.Output( "\n\n   solver for cycle...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 8 ) REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "   %6d", bconSet[i].solver[j] );
      REPORT::rpt.Output( text, 2 );
    }

    REPORT::rpt.Output( "\n\n", 2 );
    REPORT::rpt.OutputLine1( 2 );


    // report specifications for lines ---------------------------------------------------

    for( int j=0; j<bconSet[i].nofLines; j++ )
    {
      sprintf( text, "\n   boundary conditions at line %d of type %d\n\n",
                     bconSet[i].bcLine[j].no, bconSet[i].bcLine[j].kind );
      REPORT::rpt.Output( text, 2 );

      sprintf( text, "   %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n",
                     "x", "y", "U", "V", "S", "K", "D", "C", "qb" );
      REPORT::rpt.Output( text, 2 );

      for( int k=0; k<bconSet[i].bcLine[j].ndat; k++ )
      {
        BCONLINE* bcl = &bconSet[i].bcLine[j];

        sprintf( text, "   %-12.4le %-12.4le %-12.4le %-12.4le %-12.4le %-12.4le %-12.4le %-12.4le %-12.4le\n",
                       bcl->x[k], bcl->y[k],
                       bcl->U[k], bcl->V[k], bcl->S[k],
                       bcl->K[k], bcl->D[k],
                       bcl->C[k], bcl->Qb[k] );
        REPORT::rpt.Output( text, 2 );
      }
    }

    REPORT::rpt.Output( "\n", 2 );
    REPORT::rpt.OutputLine1( 2 );


    // report specifications for nodes ---------------------------------------------------

    sprintf( text, "\n   %d boundary conditions at nodes\n\n", bconSet[i].nofNodes );

    if( bconSet[i].nofNodes > 0 )
    {
      sprintf( text, "   %-4s %-6s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n",
                     "no", "type", "U", "V", "S", "K", "D", "C", "qb" );
      REPORT::rpt.Output( text, 2 );

      for( int j=0; j<bconSet[i].nofNodes; j++ )
      {
        BCON* bcn = &bconSet[i].bcNode[j];

        sprintf( text, "   %-4d %-6d %-12.4le %-12.4le %-12.4le %-12.4le %-12.4le %-12.4le %-12.4le\n",
                       bcn->no, bcn->kind,
                       bcn->val->U, bcn->val->V, bcn->val->S,
                       bcn->val->K, bcn->val->D,
                       bcn->val->C[0], bcn->val->Qb[0] );
        REPORT::rpt.Output( text, 2 );
      }
    }

    REPORT::rpt.Output( "\n", 2 );
    REPORT::rpt.OutputLine2( 2 );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////

void TIMEINT::Input_40100( ASCIIFILE *file )
{
  char *textLine;
  char  text[400];


  // -------------------------------------------------------------------------------------
  // the file is read twice

  // first reading of file:
  //   1. determine number of boundary sets
  //   2. determine number of last time step
  //   3. determine maximum number of boundary lines and nodes
  //   4. determine number of nodes & lines with periodic boundary condition

  int cntOutput = 0;

  int cntbcNode = 0;
  int maxbcNode = 0;

  int nbcline = 100;

  int* bclno     = new int[nbcline];
  int* bclkind   = new int[nbcline];
  int* cntbcline = new int[nbcline];
  int* maxbcline = new int[nbcline];
  if( !bclno || !bclkind || !cntbcline || !maxbcline )
    REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_40100 #1)" );

  memset( bclno,     0, nbcline );
  memset( bclkind,   0, nbcline );
  memset( cntbcline, 0, nbcline );
  memset( maxbcline, 0, nbcline );

  nPeriodicNode = 0;
  nPeriodicLine = 0;

  while( !feof(file->getid()) )
  {
    if( !(textLine = file->nextLine()) )  break;

    char key[50];
    int  ikey = -99;

    // remove leading blanks
    while( *textLine == ' ' )  textLine++;

    if( textLine[0] == '$' )
    {
      sscanf( textLine, "$%s", key );

      for( int i=0; i<nkey; i++ )
      {
        if( strcmp(key,datkey[i].name) == 0 )
        {
          ikey = datkey[i].id;
          break;
        }
      }
    }

    switch( ikey )
    {
      case kTM_OUTPUT:
        cntOutput++;
        break;

      // read number of first and last time step and bcLoop ------------------------------
      case kTM_STEPS:
        sscanf( textLine, "$TM_STEPS %d %d %d",
                          &firstTimeStep,
                          &lastTimeStep,
                          &bcLoop );

        if( firstTimeStep <= 0  ||  lastTimeStep <= 0 )
          REPORT::rpt.Error( "time step numbers > 0 expected (TIMEINT::Input_40100 #2)" );

        result = new int[lastTimeStep];
        if( !result )  REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_40100 #3)" );

        memset( result, 0, lastTimeStep );
        result[lastTimeStep-1] = true;
        break;

      // read start time -----------------------------------------------------------------
      case kTM_SETTIME:
        {
          char tmstr[30];
          sscanf( textLine, "$TM_SETTIME %s", tmstr );
          set = true;
          startTime.Set( tmstr );
        }
        break;

      // read length of time step in seconds ---------------------------------------------
      case kTM_INTERVAL:
        {
          double deltaTm, relaxTmFlow, relaxTmTurb, dummy;

          if( release >= 40070 )
          {
            sscanf( textLine, "$TM_INTERVAL %lf %lf %lf",
                              &deltaTm,
                              &relaxTmFlow,
                              &relaxTmTurb );
          }
          else
          {
            sscanf( textLine, "$TM_INTERVAL %lf %lf %lf %lf",
                              &deltaTm,
                              &relaxTmFlow,
                              &dummy,
                              &relaxTmTurb );
          }

          if( relaxTmFlow <= 0.0 || relaxTmFlow > deltaTm )  relaxTmFlow = deltaTm;
          if( relaxTmTurb <= 0.0 || relaxTmTurb > deltaTm )  relaxTmTurb = deltaTm;

          deltaTime.Setsec( deltaTm );
          relaxTimeFlow.Setsec( relaxTmFlow );
          relaxTimeTurb.Setsec( relaxTmTurb );
        }
        break;

      // read weighting of time step -----------------------------------------------------
      case kTM_WEIGHT:
        sscanf( textLine, "$TM_WEIGHT %lf %lf %lf", &thetaFlow, &thetaTurb, &thetaSedi );
        break;

      // determine number of boundary sets -----------------------------------------------
      case kTM_STEP_NO:
        setsOfBcon++;

        for( int i=0; i<nbcline; i++ )
        {
          if( bclno[i] == 0 )  break;
          if( cntbcline[i] > maxbcline[i] )  maxbcline[i] = cntbcline[i];
        }
        memset( cntbcline, 0, nbcline );

        if( cntbcNode > maxbcNode )  maxbcNode = cntbcNode;
        break;

      // determine number of lines with boundary conditions ------------------------------
      case kTM_BOUND_LINE:
        {
          int no;
          sscanf( textLine, "$TM_BOUND_LINE %d", &no );

          if( no > 0 )
          {
            for( int i=0; i<nbcline; i++ )
            {
              if( bclno[i] == 0 )
              {
                bclno[i] = no;
                cntbcline[i]++;
                no = -1;
                break;
              }
            }

            // increase nbcLine (=100) and the corresponding arrays
            // if the number of boundary lines is larger
            if( no > 0 )
            {
              int* new_bclno     = new int[2*nbcline];
              int* new_cntbcline = new int[2*nbcline];
              int* new_maxbcline = new int[2*nbcline];

              if( !new_bclno || !new_cntbcline || !new_maxbcline )
                REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_40100 #4)" );

              memset( new_bclno,     0, 2*nbcline );
              memset( new_cntbcline, 0, 2*nbcline );
              memset( new_maxbcline, 0, 2*nbcline );

              memcpy( new_bclno,     bclno,     nbcline );
              memcpy( new_cntbcline, cntbcline, nbcline );
              memcpy( new_maxbcline, maxbcline, nbcline );

              delete[] bclno;
              delete[] cntbcline;
              delete[] maxbcline;

              bclno     = new_bclno;
              cntbcline = new_cntbcline;
              maxbcline = new_maxbcline;

              bclno[nbcline] = no;
              cntbcline[nbcline]++;

              nbcline *= 2;
            }
          }
        }
        break;

      case kTM_LINE:
        {
          int no   = 0;
          int kind = 0;
          char key[20];

          sscanf( textLine, "$TM_LINE %s %d", key, &no );

          // search for keyword to determine the kind of boundary condition --------------

          for( int i=0; i<BCON::nkey; i++ )
          {
            if( strcmp(key,BCON::datkey[i].name) == 0 )
            {
              kind = BCON::datkey[i].id;
              break;
            }
          }

          if( no > 0 )
          {
            for( int i=0; i<nbcline; i++ )
            {
              if( bclno[i] == no  &&  bclkind[i] == kind )
              {
                cntbcline[i]++;
                no = -1;
                break;
              }
              else if( bclno[i] == 0 )
              {
                bclno[i]   = no;
                bclkind[i] = kind;
                cntbcline[i]++;
                no = -1;
                break;
              }
            }

            // increase nbcline (=100) and the corresponding arrays
            // if the number of boundary lines gets larger than 100
            if( no > 0 )
            {
              int* new_bclno     = new int[2*nbcline];
              int* new_bclkind   = new int[2*nbcline];
              int* new_cntbcline = new int[2*nbcline];
              int* new_maxbcline = new int[2*nbcline];

              if( !new_bclno || !new_bclkind || !new_cntbcline || !new_maxbcline )
                REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_40100 #4)" );

              memset( new_bclno,     0, 2*nbcline );
              memset( new_bclkind,   0, 2*nbcline );
              memset( new_cntbcline, 0, 2*nbcline );
              memset( new_maxbcline, 0, 2*nbcline );

              memcpy( new_bclno,     bclno,     nbcline );
              memcpy( new_bclkind,   bclkind,   nbcline );
              memcpy( new_cntbcline, cntbcline, nbcline );
              memcpy( new_maxbcline, maxbcline, nbcline );

              delete[] bclno;
              delete[] bclkind;
              delete[] cntbcline;
              delete[] maxbcline;

              bclno     = new_bclno;
              bclkind   = new_bclkind;
              cntbcline = new_cntbcline;
              maxbcline = new_maxbcline;

              bclno[nbcline]   = no;
              bclkind[nbcline] = kind;
              cntbcline[nbcline]++;

              nbcline *= 2;
            }
          }
        }
        break;

      // determine number of nodes with boundary conditions ------------------------------
      case kTM_BOUND_NODE:
        cntbcNode++;
        break;

      // determine number of nodes with boundary conditions ------------------------------
      case kTM_NODE:
        cntbcNode++;
        break;

      // determine number of nodes with periodic boundary condition ----------------------
      case kTM_PERIODIC_NODE:
        nPeriodicNode++;
        break;

      // determine number of lines with periodic boundary condition ----------------------
      case kTM_PERIODIC_LINE:
        nPeriodicLine++;
        break;
    }
  }

  for( int i=0; i<nbcline; i++ )
  {
    if( bclno[i] == 0 )  break;
    if( cntbcline[i] > maxbcline[i] )  maxbcline[i] = cntbcline[i];
  }

  if( cntbcNode > maxbcNode )  maxbcNode = cntbcNode;


  // -------------------------------------------------------------------------------------
  // allocate memory for some arrays

  int set  = -1;
  int ncyc =  0;

  if( setsOfBcon <= 0 )
  {
    endLoop    = 1;
    setsOfBcon = 1;
    set = 0;
  }

  bconSet = new BCONSET[setsOfBcon+1];
  if( !bconSet )
    REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #5)" );

  bconSet[0].noTimeStep = 1;


  // -------------------------------------------------------------------------------------

  for( int i=0; i<nbcline; i++ )
  {
    if( bclno[i] == 0 )
    {
      nbcline = i;
      break;
    }
  }

  for( int i=0; i<setsOfBcon; i++ )
  {
    bconSet[i].bcNode = new BCON [maxbcNode];
    BCVAL* bcval      = new BCVAL[maxbcNode];

    bconSet[i].bcLine = new BCONLINE[nbcline];

    if( !bconSet[i].bcNode || !bcval || !bconSet[i].bcLine )
      REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_40100 #6)" );

    for( int j=0; j<maxbcNode; j++ ) bconSet[i].bcNode[j].val = &bcval[j];

    for( int j=0; j<nbcline; j++ )
    {
      bconSet[i].bcLine[j].x   = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].y   = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].U   = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].V   = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].S   = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].K   = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].D   = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].C   = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].Qb  = new double         [maxbcline[j]];
      bconSet[i].bcLine[j].gct = new BCVAL::GAUGECT [maxbcline[j]];

      if(     !bconSet[i].bcLine[j].x  ||  !bconSet[i].bcLine[j].y
          ||  !bconSet[i].bcLine[j].U  ||  !bconSet[i].bcLine[j].V
          ||  !bconSet[i].bcLine[j].S  ||  !bconSet[i].bcLine[j].K
          ||  !bconSet[i].bcLine[j].D  ||  !bconSet[i].bcLine[j].C
          ||  !bconSet[i].bcLine[j].Qb ||  !bconSet[i].bcLine[j].gct )
          REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_40100 #7)" );
    }
  }

  delete[] bclno;
  delete[] bclkind;
  delete[] cntbcline;
  delete[] maxbcline;


  // -------------------------------------------------------------------------------------

  if( !result )
  {
    result = new int[lastTimeStep];
    if( !result )  REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #8)" );
  }

  if( !cntOutput )
  {
    for( int i=0; i<lastTimeStep; i++ )  result[i] = true;
  }
  else
  {
    for( int i=0; i<lastTimeStep; i++ )  result[i] = false;
  }

  reset_statist = new int[lastTimeStep];
  if( !reset_statist )  REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_40100 #3)" );
  for( int i=0; i<lastTimeStep; i++ )  reset_statist[i] = false;

  // -------------------------------------------------------------------------------------

  if( nPeriodicNode > 0 )
  {
    periodicNode[0] = new int[nPeriodicNode];
    periodicNode[1] = new int[nPeriodicNode];
    if( !periodicNode[0] || !periodicNode[1] )
      REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #9)" );
  }

  if( nPeriodicLine > 0 )
  {
    periodicLine[0] = new int[nPeriodicLine];
    periodicLine[1] = new int[nPeriodicLine];
    periodicLine[2] = new int[nPeriodicLine];
    periodicLine[3] = new int[nPeriodicLine];
    if( !periodicLine[0] || !periodicLine[1] || !periodicLine[2] || !periodicLine[3] )
      REPORT::rpt.Error( "can not allocate memory (TIMEINT::Input_30900 #10)" );
  }

  nPeriodicNode = 0;
  nPeriodicLine = 0;


  // -------------------------------------------------------------------------------------
  // second reading of file
  // -------------------------------------------------------------------------------------

  file->rewind();

  while( !feof(file->getid()) )
  {
    if( !(textLine = file->nextLine()) )  break;

    char key[50];
    int  ikey = -99;

    if( textLine[0] == '$' )
    {
      sscanf( textLine, "$%s", key );

      for( int i=0; i<nkey; i++ )
      {
        if( strcmp(key,datkey[i].name) == 0 )
        {
          ikey = datkey[i].id;
          break;
        }
      }
    }

    switch( ikey )
    {
      // read time step numbers for output -----------------------------------------------
      case kTM_OUTPUT:
        {
          int start, stop, step;
          int n = sscanf( textLine, "$TM_OUTPUT %d to %d step %d", &start, &stop, &step );

          if( n == 3 )
          {
            for( int i=start; i<=stop && i<=lastTimeStep; i+=step )  result[i-1] = true;
          }
          else
          {
            char  list[500];
            char  seps[] = " ,\t\n\r";
            char* token;

            strcpy( list, textLine+10 );
            token = strtok( list, seps );

            while( token != NULL )
            {
              int stepno;
              sscanf( token, "%d", &stepno );
              ncyc++;

              if( stepno <= lastTimeStep )  result[stepno-1] = true;

              token = strtok( NULL, seps );
            }
          }
        }
        break;

      case kTM_RESET_STATIST:
        {
          int start, stop, step;
          int n = sscanf( textLine, "$TM_RESET_STATIST %d to %d step %d", &start, &stop, &step );

          if( n == 3 )
          {
            for( int i=start; i<=stop && i<=lastTimeStep; i+=step )  reset_statist[i-1] = true;
          }
          else
          {
            char  list[500];
            char  seps[] = " ,\t\n\r";
            char* token;

            strcpy( list, textLine+10 );
            token = strtok( list, seps );

            while( token != NULL )
            {
              int stepno;
              sscanf( token, "%d", &stepno );
              ncyc++;

              if( stepno <= lastTimeStep )  reset_statist[stepno-1] = true;

              token = strtok( NULL, seps );
            }
          }
        }
        break;

      case kTM_STEP_NO:
        {
          set++;
          sscanf( textLine, "$TM_STEP_NO %d", &bconSet[set].noTimeStep );

          endLoop = bconSet[set].noTimeStep;

          if( set == 0 )
          {
            if( bconSet[set].noTimeStep != 1 )
              REPORT::rpt.Error( "initial time step number must be 1 (TIMEINT::Input_30900 #11)" );
          }

          ncyc = 0;
        }
        break;

      case kTM_CYCLE:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          strcpy( list, textLine+9 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].cycle[ncyc] );
            ncyc++;

            token = strtok( NULL, seps );
          }
        }
        break;

      case kTM_STATIONARY:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+14 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].stat[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].stat[i] = bconSet[set].stat[n-1];
        }
        break;

      case kTM_TURBULENCE:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+14 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].turb[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].turb[i] = bconSet[set].turb[n-1];
        }
        break;

      case kTM_DISPERSION:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+14 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].disp[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].disp[i] = bconSet[set].disp[n-1];
        }
        break;

      case kTM_MAXITER:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+11 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].cycit[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].cycit[i] = bconSet[set].cycit[n-1];
        }
        break;

      case kTM_SOLVER:
        {
          char  list[500];
          char  seps[] = " ,\t\n\r";
          char* token;

          int   n = 0;

          strcpy( list, textLine+10 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            sscanf( token, "%d", &bconSet[set].solver[n] );
            n++;

            token = strtok( NULL, seps );
          }

          for( int i=n; i<ncyc; i++ )  bconSet[set].solver[i] = bconSet[set].solver[n-1];
        }
        break;

      case kTM_BOUND_LINE:
        {
          int  no;
          int  id   = -1;
          long kind =  0;

          sscanf( textLine, "$TM_BOUND_LINE %d %ld", &no, &kind );

          if( no > 0 )
          {
            for( int i=0; i<bconSet[set].nofLines; i++ )
            {
              if( bconSet[set].bcLine[i].no == no )
              {
                id = i;
                break;
              }
            }

            if( id == -1 )
            {
              id = bconSet[set].nofLines;
              bconSet[set].nofLines++;
            }

            bconSet[set].bcLine[id].no   = no;
            bconSet[set].bcLine[id].kind = kind;

            int k = bconSet[set].bcLine[id].ndat;

            bconSet[set].bcLine[id].x[k]  = 0.0;
            bconSet[set].bcLine[id].y[k]  = 0.0;
            bconSet[set].bcLine[id].U[k]  = 0.0;
            bconSet[set].bcLine[id].V[k]  = 0.0;
            bconSet[set].bcLine[id].S[k]  = 0.0;
            bconSet[set].bcLine[id].K[k]  = 0.0;
            bconSet[set].bcLine[id].D[k]  = 0.0;
            bconSet[set].bcLine[id].C[k]  = 0.0;
            bconSet[set].bcLine[id].Qb[k] = 0.0;

            sscanf( textLine, "$TM_BOUND_LINE %d %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                              &no, &kind,
                              &bconSet[set].bcLine[id].x[k],
                              &bconSet[set].bcLine[id].y[k],
                              &bconSet[set].bcLine[id].U[k],
                              &bconSet[set].bcLine[id].V[k],
                              &bconSet[set].bcLine[id].S[k],
                              &bconSet[set].bcLine[id].K[k],
                              &bconSet[set].bcLine[id].D[k],
                              &bconSet[set].bcLine[id].C[k],
                              &bconSet[set].bcLine[id].Qb[k] );

            bconSet[set].bcLine[id].ndat++;
          }
        }
        break;

      case kTM_BOUND_NODE:
        {
          int  no;
          long kind = 0;
          double USpec  = 0.0,
                 VSpec  = 0.0,
                 SSpec  = 0.0,
                 KSpec  = 0.0,
                 DSpec  = 0.0,
                 CSpec  = 0.0,
                 QbSpec = 0.0;

          int j = bconSet[set].nofNodes;

          sscanf( textLine, "$TM_BOUND_NODE %d %ld %lf %lf %lf %lf %lf %lf %lf",
                            &no, &kind,
                            &USpec, &VSpec, &SSpec,
                            &KSpec, &DSpec,
                            &CSpec, &QbSpec );

          bconSet[set].bcNode[j].no         = no - 1;
          bconSet[set].bcNode[j].kind       = kind;

          bconSet[set].bcNode[j].val->U     = USpec;
          bconSet[set].bcNode[j].val->V     = VSpec;
          bconSet[set].bcNode[j].val->S     = SSpec;
          bconSet[set].bcNode[j].val->K     = KSpec;
          bconSet[set].bcNode[j].val->D     = DSpec;
          bconSet[set].bcNode[j].val->C[0]  = CSpec;
          bconSet[set].bcNode[j].val->Qb[0] = QbSpec;

          bconSet[set].nofNodes++;
        }
        break;

      // read new type of boundary conditions for lines (Rismo Rev. 4.01.00) -------------
      case kTM_LINE:
        {
          int  no;
          char key[20];

          sscanf( textLine, "$TM_LINE %s %d", key, &no );

          if( no > 0 )
          {
            int kind = 0;

            // search for keyword to determine the kind of boundary condition ----------

            for( int i=0; i<BCON::nkey; i++ )
            {
              if( strcmp(key,BCON::datkey[i].name) == 0 )
              {
                kind = BCON::datkey[i].id;
                break;
              }
            }

            if( kind <= 0 ) break;

            // search for already specified boundary values of same kind at line no ------

            BCONSET*  bcset = &bconSet[set];
            BCONLINE* bcl   = NULL;

            for( int i=0; i<bcset->nofLines; i++ )
            {
              if( bcset->bcLine[i].no == no  &&  bcset->bcLine[i].kind == kind )
              {
                bcl = &bcset->bcLine[i];
                break;
              }
            }

            if( !bcl )
            {
              bcl = &bcset->bcLine[bcset->nofLines];
              bcset->nofLines++;

              bcl->no    = no;
              bcl->kind  = kind;
              bcl->ndat = 0;
            }

            // read boundary values depending on the specified keyword -------------------

            int k = bcl->ndat;

            switch( kind )
            {
              // read flow q (=U,V) and coordinates (x,y)
              case BCON::kInlet:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf %lf",
                        key, &no,
                        &bcl->U[k], &bcl->V[k], &bcl->x[k], &bcl->y[k] );
                bcl->ndat++;
                break;

              // read water elevation (S) and coordinates (x,y)
              case BCON::kOutlet:
                bcl->gct[k].nocg = 0;
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf %d %lf",
                        key, &no,
                        &bcl->S[k], &bcl->x[k], &bcl->y[k],
                        &bcl->gct[k].nocg, &bcl->gct[k].So );
                bcl->ndat++;
                break;

              // no further values expected
              case BCON::kOpenBnd:
                bcl->ndat = 1;
                break;

              // read discharge Q(=U) and water elevation (S); only 1 value per line is
              // accepted and further readings will overwrite the first values
              case BCON::kQInlet:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf",
                        key, &no,
                        &bcl->U[k], &bcl->S[k] );
                bcl->ndat = 1;
                break;

              // read a hydrograph: time T(=x), discharge Q(=U) and water elevation (S)
              case BCON::kQTInlet:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf",
                        key, &no,
                       &bcl->U[k], &bcl->S[k],  &bcl->x[k] );
                bcl->ndat++;
                break;

              // read outlet condition: discharge Q(=U) and water elevation (S)
              case BCON::kSQOutlet:
                bcl->gct[k].nocg = 0;
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %d %lf",
                        key, &no,
                        &bcl->U[k], &bcl->S[k],
                        &bcl->gct[k].nocg, &bcl->gct[k].So );
                bcl->ndat++;
                break;

              // read outlet condition: time T(=x) and water elevation (S)
              case BCON::kSTOutlet:
                bcl->gct[k].nocg = 0;
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %d %lf",
                        key, &no,
                        &bcl->S[k], &bcl->x[k],
                        &bcl->gct[k].nocg, &bcl->gct[k].So );
                bcl->ndat++;
                break;

              // read specified flow direction (U,V) and coordinates (x,y)
              case BCON::kSlip:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf %lf",
                        key, &no,
                        &bcl->U[k], &bcl->V[k], &bcl->x[k], &bcl->y[k] );
                bcl->ndat++;
                break;

              // read specified log law parameter (available only on boundaries)
              case BCON::kLoglaw:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf",
                        key, &no, &bcl->dw, &bcl->kw );
                bcl->ndat++;
                break;

              // read specified velocites (U,V) and coordinates (x,y)
              case BCON::kSetUV:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf %lf",
                        key, &no,
                        &bcl->U[k], &bcl->V[k], &bcl->x[k], &bcl->y[k] );
                bcl->ndat++;
                break;

              // read specified water elevation (S) and coordinates (x,y)
              case BCON::kSetS:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf",
                        key, &no,
                        &bcl->S[k], &bcl->x[k], &bcl->y[k] );
                bcl->ndat++;
                break;

              // read specified turbulence k&eps (K,D) and coordinates (x,y)
              case BCON::kSetKD:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf %lf",
                        key, &no,
                        &bcl->K[k], &bcl->D[k], &bcl->x[k], &bcl->y[k] );
                bcl->ndat++;
                break;

              // read specified sediment concentration (C) and coordinates (x,y)
              case BCON::kSetC:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf",
                        key, &no,
                        &bcl->C[k], &bcl->x[k], &bcl->y[k] );
                bcl->ndat++;
                break;

              // read specified sediment rate (qb) and coordinates (x,y)
              case BCON::kRateC:
                sscanf( textLine, "$TM_LINE %s %d %lf %lf %lf",
                        key, &no,
                        &bcl->Qb[k], &bcl->x[k], &bcl->y[k] );
                bcl->ndat++;
                break;
            }
          }
        }
        break;

      // read new type of boundary conditions for nodes (Rismo Rev. 4.01.00) -------------
      case kTM_NODE:
        {
          int  no;
          int  kind;
          char key[20];

          double USpec  = 0.0;
          double VSpec  = 0.0;
          double SSpec  = 0.0;
          double KSpec  = 0.0;
          double DSpec  = 0.0;
          double CSpec  = 0.0;
          double QbSpec = 0.0;

          sscanf( textLine, "$TM_NODE %s %d", key, &no );

          // search for keyword to determine the kind of boundary condition --------------

          for( int i=0; i<BCON::nkey; i++ )
          {
            if( strcmp(key,BCON::datkey[i].name) == 0 )
            {
              kind = BCON::datkey[i].id;
              break;
            }
          }

          int j = bconSet[set].nofNodes;

          bconSet[set].nofNodes++;
          bconSet[set].bcNode[j].no   = no - 1;
          bconSet[set].bcNode[j].kind = kind;

          switch( kind )
          {
            case BCON::kInlet:
              sscanf( textLine, "$TM_NODE %s %d %lf %lf",
                                key, &no, &USpec, &VSpec );
              bconSet[set].bcNode[j].val->U = USpec;
              bconSet[set].bcNode[j].val->V = VSpec;
              break;

            case BCON::kOutlet:
              bconSet[set].bcNode[j].val->gct.nocg = 0;
              sscanf( textLine, "$TM_NODE %s %d %lf %d %lf",
                                key, &no, &SSpec,
                                &bconSet[set].bcNode[j].val->gct.nocg,
                                &bconSet[set].bcNode[j].val->gct.So);
              bconSet[set].bcNode[j].val->S = SSpec;
              break;

            case BCON::kOpenBnd:
              break;

            case BCON::kSlip:
              sscanf( textLine, "$TM_NODE %s %d %lf %lf",
                                key, &no, &USpec, &VSpec );
              bconSet[set].bcNode[j].val->U  = USpec;
              bconSet[set].bcNode[j].val->V  = VSpec;
              break;

            case BCON::kSetUV:
              sscanf( textLine, "$TM_NODE %s %d %lf %lf",
                                key, &no, &USpec, &VSpec );
              bconSet[set].bcNode[j].val->U = USpec;
              bconSet[set].bcNode[j].val->V = VSpec;
              break;

            case BCON::kSource:
              sscanf( textLine, "$TM_NODE %s %d %lf",
                                key, &no, &USpec );
              bconSet[set].bcNode[j].val->Q = USpec;
              USpec = 0.0;
              break;

            case BCON::kSetS:
              sscanf( textLine, "$TM_NODE %s %d %lf",
                                key, &no, &SSpec );
              bconSet[set].bcNode[j].val->S = SSpec;
              break;

            case BCON::kSetKD:
              sscanf( textLine, "$TM_NODE %s %d %lf %lf",
                                key, &no, &KSpec, &DSpec );
              bconSet[set].bcNode[j].val->K = KSpec;
              bconSet[set].bcNode[j].val->D = DSpec;
              break;

            case BCON::kSetC:
              sscanf( textLine, "$TM_NODE %s %d %lf",
                                key, &no, &CSpec );
              bconSet[set].bcNode[j].val->C[0] = CSpec;
              break;

            case BCON::kRateC:
              sscanf( textLine, "$TM_NODE %s %d %lf",
                                key, &no, &QbSpec );
              bconSet[set].bcNode[j].val->Qb[0] = QbSpec;
              break;
          }
        }
        break;

      // read periodic boundary conditions for nodes -------------------------------------
      case kTM_PERIODIC_NODE:
        sscanf( textLine, "$TM_PERIODIC_NODE %d %d",
                &periodicNode[0][nPeriodicNode], &periodicNode[1][nPeriodicNode] );
        nPeriodicNode++;
        break;


      // read periodic boundary conditions for nodes -------------------------------------
      case kTM_PERIODIC_LINE:
        sscanf( textLine, "$TM_PERIODIC_LINE %d %d %d %d",
                &periodicLine[0][nPeriodicLine], &periodicLine[1][nPeriodicLine],
                &periodicLine[2][nPeriodicLine], &periodicLine[3][nPeriodicLine] );
        nPeriodicLine++;
        break;
    }
  }


  // -------------------------------------------------------------------------------------
  // report parameters for time integration

  REPORT::rpt.Output( "\n", 2 );
  REPORT::rpt.OutputLine2( 2 );

  sprintf( text, "\n  %30s  %-2d\n  %30s  %-4d\n  %30s  %-4d\n",
          "sets of boundary conditions:", setsOfBcon,
          "first time step:",             firstTimeStep,
          "last time step:",              lastTimeStep );
  REPORT::rpt.Output( text, 2 );

  sprintf( text, "  %30s  %-4d\n  %30s  %-s\n",
          "repetition from:",             bcLoop,
          "time interval:",               deltaTime.Get() );
  REPORT::rpt.Output( text, 2 );

  sprintf( text, "  %30s  %-s\n  %30s  %-s\n",
          "relaxed time for flow:",       relaxTimeFlow.Get(),
          "relaxed time for turbulence:", relaxTimeTurb.Get() );
  REPORT::rpt.Output( text, 2 );

  REPORT::rpt.Output( "\n", 2 );
  REPORT::rpt.OutputLine1( 2 );

  sprintf( text, "\n  %30s  %-4.2lf\n  %30s  %-4.2lf\n  %30s  %-4.2lf\n",
          "time weighting THETA,  flow:", thetaFlow,
          "                 turbulence:", thetaTurb,
          "                   sediment:", thetaSedi );
  REPORT::rpt.Output( text, 2 );


  // report time step numbers for output -------------------------------------------------

  REPORT::rpt.Output( "\n",2 );
  REPORT::rpt.OutputLine1( 2 );

  sprintf( text, "\n   output will be written for time steps ...\n" );
  REPORT::rpt.Output( text, 2 );

  for( int i=0, j=0; i<lastTimeStep; i++ )
  {
    if( this->result[i] )
    {
      sprintf( text, "   %5d", i + 1 );
      REPORT::rpt.Output( text, 2 );
      j++;

      if ( j % 10 == 0 )  REPORT::rpt.Output( "\n", 2 );
    }
  }

  // report time step numbers for reset of statistics -------------------------------------------------
  if (reset_statist)
  {
    REPORT::rpt.Output( "\n",2 );
    REPORT::rpt.OutputLine1( 2 );

    sprintf( text, "\n   statistics will be reseted for time steps ...\n" );
    REPORT::rpt.Output( text, 2 );

    for( int i=0, j=0; i<lastTimeStep; i++ )
    {
      if( this->reset_statist[i] )
      {
        sprintf( text, "   %5d", i + 1 );
        REPORT::rpt.Output( text, 2 );
        j++;

        if ( j % 10 == 0 )  REPORT::rpt.Output( "\n", 2 );
      }
    }
  }


  REPORT::rpt.Output( "\n\n",2 );
  REPORT::rpt.OutputLine2( 2 );


  // report specifications for boundary conditions ---------------------------------------

  for( int i=0; i<setsOfBcon; i++ )
  {
    sprintf( text, "\n   %24s  %-d\n   %24s  %-d\n   %24s  %-d\n",
                   "specifications for set:", bconSet[i].noTimeStep,
                   "       number of lines:", bconSet[i].nofLines,
                   "       number of nodes:", bconSet[i].nofNodes );
    REPORT::rpt.Output( text, 2 );


    // report cycle sequence -------------------------------------------------------------

    REPORT::rpt.Output( "\n   sequence of cycles ...\n\n", 2 );

    for( int j=0; j<kMaxCycles; j++ )
    {
      ncyc = j;
      if( bconSet[i].cycle[j] == 0 )  break;

      if( j == 8 )  REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "   %6d", bconSet[i].cycle[j] );
      REPORT::rpt.Output( text, 2 );
    }


    // report turbulence model according to cycle switch ---------------------------------

    REPORT::rpt.Output( "\n\n   turbulence model ...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 8 ) REPORT::rpt.Output( "\n", 2 );

      int turb = bconSet[i].turb[j];

      sprintf( text, "   %6d", turb );
      REPORT::rpt.Output( text, 2 );

      bconSet[i].turb[j] = 0;

      switch( turb%10 )
      {
        case 1:   bconSet[i].turb[j] |= BCONSET::kVtConstant;                     break;
        case 2:   bconSet[i].turb[j] |= BCONSET::kVtAlgebraic;                    break;
        case 3:   bconSet[i].turb[j] |= BCONSET::kVtAlgShear;                     break;
        case 4:   bconSet[i].turb[j] |= BCONSET::kVtMixingLength;                 break;
        case 5:   bconSet[i].turb[j] |= BCONSET::kVtLES;                          break;
        case 6:   bconSet[i].turb[j] |= BCONSET::kVtPrandtlKol;                   break;
        case 7:   bconSet[i].turb[j] |= BCONSET::kVtLES | BCONSET::kVtAlgebraic;  break;
      }

      if( (turb/10)%10 == 1 )
      {
        bconSet[i].turb[j] |= BCONSET::kVtAnisotrop;
      }

      if( (turb/100)%10 & 1 )
      {
        bconSet[i].turb[j] |= BCONSET::kVtMin;
      }

      if( (turb/100)%10 & 2 )
      {
        bconSet[i].turb[j] |= BCONSET::kVtIterat;
      }
    }


    // report dispersion model according to cycle switch ---------------------------------

    REPORT::rpt.Output( "\n\n   dispersion model ...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 8 )  REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "   %6d", bconSet[i].disp[j] );
      REPORT::rpt.Output( text, 2 );
    }


    // report maximum number of iterations per cycle -------------------------------------

    REPORT::rpt.Output( "\n\n   iterations per cycle...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 8 )  REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "   %6d", bconSet[i].cycit[j] );
      REPORT::rpt.Output( text, 2 );
    }


    // report solver type according to cycle switch --------------------------------------

    REPORT::rpt.Output( "\n\n   solver for cycle...\n\n", 2 );

    for( int j=0; j<ncyc; j++ )
    {
      if( j == 8 ) REPORT::rpt.Output( "\n", 2 );

      sprintf( text, "   %6d", bconSet[i].solver[j] );
      REPORT::rpt.Output( text, 2 );
    }

    REPORT::rpt.Output( "\n\n", 2 );
    REPORT::rpt.OutputLine1( 2 );


    // report specifications for lines ---------------------------------------------------

    for( int j=0; j<bconSet[i].nofLines; j++ )
    {
      BCONLINE* bcl = &bconSet[i].bcLine[j];

      sprintf( text, "\n   boundary conditions at line %d (type %d)\n\n",
               bcl->no, bcl->kind );
      REPORT::rpt.Output( text, 2 );

      for( int k=0; k<bconSet[i].bcLine[j].ndat; k++ )
      {
        switch( bcl->kind )
        {
          case BCON::kInlet:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le %-12.4le\n",
                     "INLET", bcl->no, "U,V,x,y", bcl->U[k], bcl->V[k], bcl->x[k], bcl->y[k] );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kOutlet:
            if( bcl->gct[k].nocg > 0 )
            {
              sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le %-10d %-12.4le\n",
                       "OUTLET", bcl->no, "S,x,y,G,So", bcl->S[k], bcl->x[k], bcl->y[k],
                                 bcl->gct[k].nocg, bcl->gct[k].So );
            }
            else
            {
              sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le\n",
                       "OUTLET", bcl->no, "S,x,y", bcl->S[k], bcl->x[k], bcl->y[k] );
            }
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kOpenBnd:
            sprintf( text, "   %-10s %3d\n",
                     "OPEN", bcl->no  );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kQInlet:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le\n",
                     "QINLET", bcl->no, "Q,S", bcl->U[k], bcl->S[k] );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kQTInlet:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le\n",
                     "QTINLET", bcl->no, "t,Q,S", bcl->x[k], bcl->U[k], bcl->S[k] );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kSQOutlet:
            if( bcl->gct[k].nocg > 0 )
            {
              sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-10d %-12.4le\n",
                       "SQOUTLET", bcl->no, "Q,S,G,So", bcl->U[k], bcl->S[k],
                                   bcl->gct[k].nocg, bcl->gct[k].So );
            }
            else
            {
              sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le\n",
                       "SQOUTLET", bcl->no, "Q,S", bcl->U[k], bcl->S[k] );
            }
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kSTOutlet:
            if( bcl->gct[k].nocg > 0 )
            {
              sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-10d %-12.4le\n",
                       "STOUTLET", bcl->no, "t,S,G,So", bcl->x[k], bcl->S[k],
                                   bcl->gct[k].nocg, bcl->gct[k].So );
            }
            else
            {
              sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le\n",
                       "STOUTLET", bcl->no, "t,S", bcl->x[k], bcl->S[k] );
            }
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kSlip:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le %-12.4le\n",
                     "SLIP", bcl->no, "U,V,x,y", bcl->U[k], bcl->V[k], bcl->x[k], bcl->y[k] );
            REPORT::rpt.Output( text, 2 );
            break;

          // read specified velocites (U,V) and coordinates (x,y)
          case BCON::kSetUV:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le %-12.4le\n",
                     "SET_UV", bcl->no, "U,V,x,y", bcl->U[k], bcl->V[k], bcl->x[k], bcl->y[k] );
            REPORT::rpt.Output( text, 2 );
            break;

          // read specified water elevation (S) and coordinates (x,y)
          case BCON::kSetS:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le\n",
                     "SET_S", bcl->no, "S,x,y", bcl->S[k], bcl->x[k], bcl->y[k] );
            REPORT::rpt.Output( text, 2 );
            break;

          // read specified turbulence k&eps (K,D) and coordinates (x,y)
          case BCON::kSetKD:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le %-12.4le\n",
                     "SET_KD", bcl->no, "K,D,x,y", bcl->K[k], bcl->D[k], bcl->x[k], bcl->y[k] );
            REPORT::rpt.Output( text, 2 );
            break;

          // read specified sediment concentration (C) and coordinates (x,y)
          case BCON::kSetC:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le\n",
                     "SET_C", bcl->no, "C,x,y", bcl->C[k], bcl->x[k], bcl->y[k] );
            REPORT::rpt.Output( text, 2 );
            break;

          // read specified sediment rate (qb) and coordinates (x,y)
          case BCON::kRateC:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le %-12.4le\n",
                     "SET_Qb", bcl->no, "Qb,x,y", bcl->Qb[k], bcl->x[k], bcl->y[k] );
            REPORT::rpt.Output( text, 2 );
            break;
        }
      }
    }

    REPORT::rpt.Output( "\n", 2 );
    REPORT::rpt.OutputLine1( 2 );


    // report specifications for nodes ---------------------------------------------------

    sprintf( text, "\n   %d boundary conditions at nodes\n\n", bconSet[i].nofNodes );
    REPORT::rpt.Output( text, 2 );

    if( bconSet[i].nofNodes > 0 )
    {
      for( int j=0; j<bconSet[i].nofNodes; j++ )
      {
        BCON* bcn = &bconSet[i].bcNode[j];

        switch( bcn->kind )
        {
          case BCON::kInlet:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le\n",
                     "INLET", bcn->no+1, "U,V", bcn->val->U, bcn->val->V );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kOutlet:
            if( bcn->val->gct.nocg > 0 )
            {
              sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-10d %-12.4le\n",
                       "OUTLET", bcn->no+1, "S,G,So", bcn->val->S,
                                 bcn->val->gct.nocg, bcn->val->gct.So );
            }
            else
            {
              sprintf( text, "   %-10s %3d   %-15s  %-12.4le\n",
                       "OUTLET", bcn->no+1, "S", bcn->val->S );
            }
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kOpenBnd:
            sprintf( text, "   %-10s %3d\n",
                     "OPEN", bcn->no+1  );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kSlip:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le\n",
                     "SLIP", bcn->no+1, "U,V", bcn->val->U, bcn->val->V );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kSetUV:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le\n",
                     "SET_UV", bcn->no+1, "U,V", bcn->val->U, bcn->val->V );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kSetS:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le\n",
                     "SET_S", bcn->no+1, "S", bcn->val->S );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kSetKD:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le %-12.4le\n",
                     "SET_KD", bcn->no+1, "K,D", bcn->val->K, bcn->val->D );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kSetC:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le\n",
                     "SET_C", bcn->no+1, "C", bcn->val->C[0] );
            REPORT::rpt.Output( text, 2 );
            break;

          case BCON::kRateC:
            sprintf( text, "   %-10s %3d   %-15s  %-12.4le\n",
                     "SET_Qb", bcn->no+1, "Qb", bcn->val->Qb[0] );
            REPORT::rpt.Output( text, 2 );
            break;
        }
      }
    }

    REPORT::rpt.Output( "\n", 2 );
    REPORT::rpt.OutputLine2( 2 );
  }
}
