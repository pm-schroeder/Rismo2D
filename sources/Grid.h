// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// D R Y R E W
// G R I D
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Grid.h         : definition file of the class.
// Grid.cpp       : implementation file of the class.
//
// Arfact.cpp     : method  GRID::AreaFactors()
// Check.cpp      : method  GRID::Check()
// Connect.cpp    : method  GRID::Connection()
// Courant.cpp    : method  GRID::ReportCuPe()
// Dispersion.cpp : method  GRID::Dispersion()
// DryRewet.cpp   : methods GRID::Dry()
//                          GRID::Rewet()
//                          GRID::ReportDry()
//                          GRID::DryRewet()
//                          GRID::RewetDry()
// EddyDisp.cpp   : method  GRID::EddyDisp()
// Init.cpp       : method  GRID::InitKD()
// InitS.cpp      : method  GRID::InitS()
// Lumped.cpp     : method  GRID::LumpedMassMatrix()
// SlipFlow.cpp   : method  GRID::SetSlipFlow()
// Smooth.cpp     : methods GRID::SmoothS()
//                          GRID::SmoothKD()
// Turbulence.cpp : method  GRID::Turbulence()
// VeloGrad.cpp   : method  GRID::VeloGrad()
//
// Interpol.cpp   : method  DRYREW::interpolate()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a finite element grid and the included class DRYREW to manage
// the flooding an drying of parts of the grid.
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
//  16.02.2013    sc    rewetting of nodes in DryRewet() and RewetDry() adapted
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GRID_INCL
#define GRID_INCL

#include "Defs.h"

class ELEM;
class NODE;
class PROJECT;
class SECTION;
class SUBDOM;
class TIME;
class TYPE;


// parameters for dry/rewet algorithm --------------------------------------------------------------

class DRYREW
{
  public:
    int    method;             // method for drying and rewetting
    int    dryRewFreq;         // dry/rewet frequency during iterations
    int    rewetPasses;        // number of passes to rewet elements
    int    countDown;          // countDown for elements to be rewetted
    double rewetLimit;         // limit for flow depth to rewet nodes
    double dryLimit;           // limit for flow depth to eliminate nodes

  public:
    DRYREW()
    {
      method       = 2;
      dryRewFreq   = 10;
      dryLimit     = 0.001;
      rewetLimit   = 0.005;
      rewetPasses  = 100;
      countDown    = 5;
    }

    double interpolate( NODE*, int, int );
};


class GRID
{
  private:
    int    np;
    NODE*  node;

    int    ne;
    ELEM*  elem;

  public:
    static DRYREW dryRew;
    int    firstDryRew;

  public:
    // Grid.cpp ------------------------------------------------------------------------------------
    GRID();
    ~GRID();

    void   KillNode();
    void   KillElem();

    int    Getnp();
    NODE*  Getnode( int n );
    void   Setnode( int np, NODE* node );

    int    Getne();
    ELEM*  Getelem( int e );
    void   Setelem( int ne, ELEM* elem );

    void   Alloc( int np, int ne );
    void   Free();

    void   InputRegion( char* fileName, SUBDOM* subdom );
    void   InputControl( char* fileName, GRID* grid, SUBDOM* subdom );
    void   InputInitial( int isAscii, char* name, TIME* time, SUBDOM* subdom, int zb_init );
    void   InitPrevious();
    void   ScanElement( ELEM*, char*, char*, int*, int*, int* );

    void   OutputData( PROJECT* project, int timeStep, char* time );
    void   OutputGeom( PROJECT* project, int timeStep );

    // ArFact.cpp ----------------------------------------------------------------------------------
    void   AreaFactors();

    // Check.cpp -----------------------------------------------------------------------------------
    void   Check();

    // Connect.cpp ---------------------------------------------------------------------------------
    void   Connection( long );

    // Courant.cpp ---------------------------------------------------------------------------------
    double ReportCuPe( double, double );

    // Dispersion.cpp ------------------------------------------------------------------------------
    void   Dispersion( PROJECT* project, double* Duu, double* Dvv,
                       double* Duv, double* Vsec );

    // DryRewet.cpp --------------------------------------------------------------------------------
    int    Dry( double, int );
    int    Rewet( double, int, PROJECT* );
    void   ReportDry( PROJECT*, double, int );
    void   DryRewet(double dryLimit, double rewetLimit, int countDown, int *dried, int *wetted );
    void   RewetDry(double dryLimit, double rewetLimit, int countDown, int *dried, int *wetted );

    // EddyDisp.cpp --------------------------------------------------------------------------------
    void   EddyDisp();

    // Turbulence.cpp ------------------------------------------------------------------------------
    void   Turbulence( PROJECT* );

    // Init.cpp ------------------------------------------------------------------------------------
    void   InitKD( PROJECT* );

    // InitS.cpp -----------------------------------------------------------------------------------
    void   InitS( int, SECTION* );

    // Lumped.cpp ----------------------------------------------------------------------------------
    void   LumpedMassMatrix( double** );

    // SlipFlow.cpp --------------------------------------------------------------------------------
    void   SetSlipFlow();

    // Smooth.cpp ----------------------------------------------------------------------------------
    void   SmoothS( int );
    void   SmoothKD( int );

    // VeloGrad.cpp --------------------------------------------------------------------------------
    void VeloGrad( PROJECT* p, double* dUdx, double* dUdy, double* dVdx, double* dVdy );
};

#endif
