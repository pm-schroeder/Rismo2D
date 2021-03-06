// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// T I M E I N T
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Timeint.h   : definition file of the class.
// Timeint.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the interface to the time step file.
//
// -------------------------------------------------------------------------------------------------
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
// -------------------------------------------------------------------------------------------------
//
// HISTORY
//
//    date              changes
// ------------  ----  -----------------------------------------------------------------------------
//  01.01.1992    sc    first implementation / first concept
//  29.03.2010    sc    Rismo-Version 4.01.00, new keywords: kTM_NODE, kTM_LINE
//                      class RELOC to ensure backward compatibility
//  13.10.2012    sc    Rismo-Version 4.03.00, new keyword: kTM_STATIONARY
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TIMEINT_INCL
#define TIMEINT_INCL

#include "Datkey.h"
#include "Bcon.h"
#include "Times.h"

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
