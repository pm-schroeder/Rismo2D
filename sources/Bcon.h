// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// B C O N
// B C O N L I N E
// B C O N S E T
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Bcon.h       : definition file of the 3 classes BCON, BCONLINE and BCONSET.
// Bcon.cpp     : implementation of the class BCON
// Bconline.cpp : implementation of the class BCONLINE
// Bconset.cpp  : implementation of the class BCONSET
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// These 3 classes implement boundary conditions.
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
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.1992     sc     first implementation
// 25.12.2009     ko     implementation of diffuse Q (Source/Sink)
// 27.03.2010     sc     class BCVAL appended
// 29.03.2010     sc     Rismo-Version 4.01.00: flags for boundary conditions:
//                       changed order, new keywords introduced in TIMINT::InputXXXXXX()
// 12.12.2013     sc     Rismo-Version 4.05.17: new boundary conditions for nodes
//                       TARGET_S and TARGET_ST to control the outlet water elevation
//                       from a different node
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef BCON_INCL
#define BCON_INCL

#include "Datkey.h"
#include "Vars.h"

class NODE;
class ELEM;
class MODEL;
class PROJECT;
class SUBDOM;
class TIME;

class BCVAL
{
  public:
    //                        -- specified boundary values
    double  U;                // flow velocity in x-direction
    double  V;                // flow velocity in y-direction
    double  Q;                // discharge
    double  S;                // water surface elevation

    double  K;                // turbulent kinetic energy
    double  D;                // turbulent kinetic energy dissipation

    double  dw;               // wall distance (log law boundaries)
    double  kw;               // roughness coefficient in log law (slip boundaries)

    unsigned int fracs;       // number of sediment fractions (not yet implemented)
    double* C;                // sediment concentration
    double* Qb;               // sediment transportation rate

    double  A;                // area of connected wet elements at a node
                              // parameter internally used for sink/source implementation

    struct GAUGECT            // gauge date to control outlet condition
    {
      public:
        int    node;          // name of controlled node on outlet boundary
        int    nocg;          // name of controlling gauge
        double So;            // targeted water level
        double dS;            // iterated difference to the specified outlet BCVAL::S
    };

    GAUGECT gct;

  public:
    BCVAL()
    {
      U     = 0.0;
      V     = 0.0;
      S     = 0.0;
      Q     = 0.0;
      K     = 0.0;
      D     = 0.0;

      dw    = 0.0;
      kw    = 0.0;

      fracs = 1;
      C     = new double[1];
      Qb    = new double[1];
      C[0]  = 0.0;
      Qb[0] = 0.0;

      A     = 0.0;
    }

    ~BCVAL()
    {
      if( C )  delete[] C;
      if( Qb ) delete[] Qb;
    }
};

class BCON
{
  public:
    static int    nkey;
    static DATKEY datkey[];

    enum
    {                              // --- user flags for boundary conditions
      kNone       =        0l,

      kInlet      =        1l,     //   0: inlet boundary
      kOutlet     =        2l,     //   1: outlet: water elevation specified
      kOpenBnd    =        4l,     //   2: open boundary; no boundary conditions

      kQInlet     =        8l,     //   3: discharge along control line (inlet)
      kQTInlet    =       16l,     //   4: Q(t)-relation for inlet

      kSQOutlet   =       32l,     //   5: H(Q)-relation for outlet
      kSTOutlet   =       64l,     //   6: S(t)-relation for outlet

      kSlip       =      128l,     //   7: slip velocity; flow direction
      kLoglaw     =      256l,     //   8: slip velocity; log law

      kSetUV      =      512l,     //   9: Dirichlet velocity (U,V)
      kSetS       =     1024l,     //  10: Dirichlet flow depth h

      kSource     =     2048l,     //  11: sink/source Q

      kSetKD      =     4096l,     //  12: Dirichlet K and D

      kSetC       =     8192l,     //  13: Dirichlet sediment concentration C [kg/m3]
      kRateC      =    16384l,     //  14: Neumann; sediment transport rate [m3/s]

                                   // --- internal flags for boundary conditions
      kFixU       =    32768l,     //  15: fixed value for velocity U
      kFixV       =    65536l,     //  16: fixed value for velocity V
      kFixK       =   131072l,     //  17: fixed value for turbulent kinetic energy
      kFixD       =   262144l,     //  18: fixed value for dissipation of K
      kFix_1      =   524288l,     //  19:
      kFix_2      =  1048576l,     //  20:
      kFix_3      =  2097152l,     //  21:
      kAutoSlip   =  4194304l,     //  22: slip velocity at boundary elements
      kAutoKD     =  8388608l      //  23: apply boundary model for K and D
    };


    // Switches to determine the rotation matrix from the kind of boundary.
    // The value is stored in the 3-bit unsigned int BCBAL::rot.
    enum
    {
      kTangentFlowRot = 1,         // tangential flow on closed boundaries
      kNormalFlowRot  = 2,         // normal flow on open boundaries
      kSlipFlowRot    = 3          // specified flow direction
    };

    unsigned int kind   : 28;
    unsigned int mark   :  1;
    unsigned int rot    :  3;      // flag to determine rotation matrix

  public:
    int     no;                    // node number (3D-mesh: surface node)
    double  nx, ny;                // normal direction on closed boundaries
    double  niox, nioy;            // normal direction on open boundaries (in-/outlet)
    BCVAL*  val;                   // pointer to boundary values if present

  public:
    BCON();
    ~BCON();

    int    Consistent( int kind );      // check for consistence of boundary conditions
    double Getrot( int r, int c );      // determine element (r,c) of rotation matrix
};


class BCONLINE : public BCON
{
  public:
    int             ndat;     // number of data values
    double         *x;        // specified data along line
    double         *y;
    double         *U;        // velocity
    double         *V;
    double         *S;        // water surface
    double         *K;        // turbulence
    double         *D;
    double         *C;        // sediment concentration
    double         *Qb;       // bed load transport rate
    BCVAL::GAUGECT *gct;

    double          dw;       // wall distance in log law
    double          kw;       // roughness height in log law

  public:
    BCONLINE();
    ~BCONLINE();

    int   GenBcon( MODEL *model, TIME *at, int b, BCON *bcon, double *preQ );
    int   FindInterval( double x, double y, int n, double *xar, double *yar, double *l );
    BCON* FindBcon( int no, int b, BCON *bcon );
    int   Inlet( MODEL *model, TIME *at, int b, BCON *bcon );
    int   Outlet( MODEL *model, TIME *at, int b, BCON *bcon, double *preQ );
    int   Further( MODEL *model, TIME *at, int b, BCON *bcon, long *ki );
};


class BCONSET
{
  public:
    enum
    {                      // --- type of turbulence implementation
      kVtConstant     =   1l,       // constant eddy viscosity
      kVtAlgebraic    =   2l,       // algebraic model for eddy viscosity
      kVtAlgShear     =   4l,       // algebraic model for eddy viscosity with production
      kVtMixingLength =   8l,       // Prandtl's mixing length model
      kVtLES          =  16l,       // large eddy model
      kVtPrandtlKol   =  32l,       // Prandtl-Kolmogorov equation for eddy viscosity

      kVtAnisotrop    =  64l,       // anisotrop eddy viscosity (Elder's model)

                                // --- lower limit
      kVtMin          = 128l,       // specified minimum
      kVtIterat       = 256l        // adapt eddy viscosity in iterations
    };

    int       noTimeStep;         // actual time step number
    double    time;               // time of time step

    int       nbcbuf;             // size of array bc
    int       nbc;                // number of specified boundary conditions
    BCON*     bc;                 // pointer to array of boundary conditions
    BCVAL*    bcval;

    int       nofLines;           // number of boundary-conditions at lines
    BCONLINE* bcLine;             // specified boundary conditions at lines

    int       nofNodes;           // number of boundary-conditions at nodes
    BCON*     bcNode;             // specified boundary conditions at nodes

    int       cycle [kMaxCycles]; // sequence of cycles
    int       stat  [kMaxCycles]; // solve stationary equations according to cycles
    int       turb  [kMaxCycles]; // turbulence model according to cycles
    int       disp  [kMaxCycles]; // dispersion model according to cycles
    int       cycit [kMaxCycles]; // iterations per cycle
    int       solver[kMaxCycles]; // solver for cycle

    int             ngct;         // number of gauge controlled nodes on outlet
    BCVAL::GAUGECT *gct;          // list of gauge controlled nodes on outlet

  public:
    BCONSET();
    ~BCONSET();

    // Bconset.cpp -----------------------------------------------------------------------
    BCON*  GetBcon( int );
    void   InitBcon( PROJECT* project, TIME* actualTime, double* =NULL );
    double Loglaw( double Us, double dw, double kw, double ka, double vk, double g );
};

#endif
