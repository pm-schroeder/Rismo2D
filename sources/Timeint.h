// ======================================================================================
//                                   T I M E I N T
// ======================================================================================
// This class implements the interface to the time step file.
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
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.1992     sc     first implementation / first concept
// 29.03.2010     sc     Rismo-Version 4.01.00, new keywords: kTM_NODE, kTM_LINE
//                       class RELOC to ensure backward compatibility
// 13.10.2012     sc     Rismo-Version 4.03.00, new keyword: kTM_STATIONARY
//
// ======================================================================================

// read this header file only once ...

#ifndef TIMEINT_INCL
#define TIMEINT_INCL

#include "Datkey.h"
#include "Bcon.h"
#include "Time.h"

class ASCIIFILE;

// -------------------------------------------------------------------------------------
// class RELOC relocates old and new flags of boundary conditions

class RELOC
{
  public:
    enum
    {
      kInflow    =     1l,    //   1: inlet boundary
      kOutflow   =     2l,    //   2: outlet: water elevation specified
      kOpenflow  =     4l,    //   3: open boundary; no boundary conditions
      kSlip      =     8l,    //   4: slip velocity; flow direction specified
      kFixVelo   =    16l,    //   5: fixed velocity (U,V)
      kFixS      =    32l,    //   6: fixed value for flow depth h
      kFixKD     =    64l,    //   7: fixed values for K and D
      kFlowC     =   128l,    //   8: sediment transport rate [m3/s]
      kQInflow   =   256l,    //   9: discharge along control line (inlet)
      kQTInflow  =   512l,    //  10: Q(t)-relation for inlet
      kSQOutflow =  1024l,    //  11: H(Q)-relation for outlet
      kSTOutflow =  2048l,    //  12: S(t)-relation for outlet
      kFixC      =  4096l,    //  13: fixed values for concentration [kg/m3]
      kQSource   =  8192l     //  14: diffuse sink/source Q
    };

    int old_BCON;
    int new_BCON;

  public:
    long RelocBcon( long o, int nreloc, RELOC* reloc )
    {
      long n = 0l;

      for( int i=0; i<nreloc; i++ )
      {
        if( isFS(o, reloc[i].old_BCON) )  SF( n, reloc[i].new_BCON );
      }

      return n;
    }
};


class TIMEINT
{
  private:
    int     nkey;
    DATKEY* datkey;

    enum
    {
      kTM_STEPS,         kTM_INTERVAL,      kTM_WEIGHT,      kTM_OUTPUT,
      kTM_SETTIME,       kTM_STEP_NO,       kTM_CYCLE,       kTM_STATIONARY,
      kTM_TURBULENCE,    kTM_DISPERSION,    kTM_MAXITER,     kTM_SOLVER,
      kTM_BOUND_NODE,    kTM_BOUND_LINE,    kTM_NODE,        kTM_LINE,
      kTM_PERIODIC_NODE, kTM_PERIODIC_LINE, kTM_RESET_STATIST
    };

    int      release;

    int      nreloc;
    RELOC*   reloc;

  public:
    int      setsOfBcon;         // number of sets of boundary conditions
    BCONSET* bconSet;
    BCONSET* actualBcSet;

    int      firstTimeStep;      // number of first time step
    int      lastTimeStep;       // number of last time step
    int      bcLoop;             // time step number to start repetition
    int      endLoop;            // internal: last time step with bc

    TIME     deltaTime;
    TIME     relaxTimeFlow;
    TIME     relaxTimeTurb;

    TIME     prevTime;
    TIME     nextTime;
    TIME     actualTime;

    TIME     incTime;

    int      set;
    TIME     startTime;

    double   thetaFlow;            // time weighting, flow equation
    double   thetaTurb;            // time weighting, turbulence equation
    double   thetaSedi;            // time weighting, sediment equation

    int*     result;               // array of time steps for output

    int      nPeriodicNode;        // number of node pairs with periodic boundary condition
    int*     periodicNode[2];      // list of node pairs

    int      nPeriodicLine;        // number of line pairs with periodic boundary condition
    int*     periodicLine[4];      // list of line & node pairs, node = first point of line
                                   // [0] = source line number      / [1] = first node of [0]
                                   // [2] = destination line number / [3] = first node of [2]

    int*     reset_statist;        // array of time steps on which stitisic values will be resetted

  public:
    TIMEINT();
    ~TIMEINT();

    void Input( char *fileName );
    void Input_00000( ASCIIFILE *file );
    void Input_30900( ASCIIFILE *file );
    void Input_40100( ASCIIFILE *file );
};

#endif
