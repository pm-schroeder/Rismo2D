// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// N O D E
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Node.h   : definition file of the class.
// Node.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a Finite Element Node.
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
//  01.09.2004    sc    class SUB implemented for parallel computing
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef NODE_INCL
#define NODE_INCL

#include "Defs.h"
#include "Vars.h"
#include "Bcon.h"

class  ELEM;

class SUB
{
  public:
    int no;
    int dry;

    SUB* next;

  public:
    SUB()
    {
      no   = 0;
      dry  = false;
      next = NULL;
    };

    ~SUB()
    { };
};

class NODE
{
  private:
    int no;                             // node number
    int name;                           // node name

  public:
    enum {
      kInlet        =          1,       //  0: inflow boundary
      kOutlet       =          2,       //  1: outflow; water elevation specified
      kOpenBnd      =          4,       //  2: open boundary; no boundary conditions
      kGridBnd      =          8,       //  3: boundary node of the whole grid
      kBound        =         16,       //  4: boundary node of the wet grid
      kRotat        =         32,       //  5: apply rotation to boundary node

      kDry          =         64,       //  6: dry flag
      kDryPrev      =        128,       //  7:
      kMarsh        =        256,       //  8:
      kMarshPrev    =        512,       //  9:

      kCornNode     =       1024,       // 10:
      kMidsNode     =       2048,       // 11:
      kNoMoment     =       4096,       // 12:

      kInface       =       8192,       // 13: interface
      kInface_UP    =      16384,       // 14: upstream interface
      kInface_DN    =      32768        // 15: downstream interface
    };

  public:
    unsigned int flag         : 16;     // general flags
    unsigned int noel         : 12;     // number of connected elements
    unsigned int mark         :  1;     // mark elements/nodes (temporary used)
    unsigned int countDown    : 12;     // countDown for rewetting
    unsigned int fixqb        :  1;     //

    SUB*   sub;                         // pointer to list of connected subdomains
    double x, y;                        // horizontal coordinates
    double z;                           // elevation, may differ to zor at dry nodes
    double zor;                         // elevation

    double cf;                          // bottom friction coefficient
    double cfw;                         // wall friction coefficient
    double vt;                          // eddy viscosity
    double exx, exy, eyy;               // eddy diffusivity
    double uu, uv, vv;                  // Reynolds stresses
    double Dxx, Dxy, Dyy;               // dispersion
    double Vsec;                        // secondary flow velocity at bottom

    double dz;                          // accumulated change of bottom elevation
    double zero;                        // elevation of unerodible bed ( zero < zor )
    double qbo;                         // transport rate in zero state

    VARS   v;                           // variables at actual time step
    VARS   vo;                          // variables at previous time step
    BCON   bc;                          // boundary conditions
    ELEM** el;                          // list of connected elements

  public:
    NODE();
    ~NODE();

    NODE operator =( const NODE& );

    void Setno( int no )     { this->no = no; };
    int  Getno()             { return this->no; };

    void Setname( int name ) { this->name = name; };
    int  Getname()           { return this->name; };
};

#endif
