// ======================================================================================
//                                        S E D
// ======================================================================================
// This class implements the sediment class.
// ======================================================================================
//
// Copyright (C) 2005-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software.
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 31.10.2005     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once ...

#ifndef SED_INCL
#define SED_INCL

#include "Elem.h"

class SED
{
  private:
    int     isinit;

  public:
    enum
    {
      kLinear = 1, kQuadratic = 2
    };

    double* qbc;              // array to hold transport capacity (equilibrium)
    double* PLs;              // array to hold nonequilibrium parameters Ls, Pls
    double* sx;               // array to hold flow direction
    double* sy;
    double* dzds;             // bed slope in flow direction (+/-)
    double* dzmx;             // maximum of bed slope in flow direction (+)
    double* dhds;             // change of flow depth in flow direction (+)

  public:
    int    nfrac;             // number of sediment fractions (place holder)

    double rhob;              // density of grain

    double M;                 // proportional factor for rate of erosion
    double tauc;              // critical friction for erosion
    double taus;              // critical friction for sedimentation
    double us;                // sinking velocity

    double d50;               // characteristic grain diameter (=d50)
    double d90;               // grain diameter for roughness height (=d90)
    double por;               // porosity
    double phir;              // friction angle; internal use: tan(phir)

    double minQb;             // lower limit for qb (qb = 0 for qb < minQb)
    double convQb;            // convergence limit for qb

    double maxDz;             // maximum change of bottom elevation
                              // to determine the morphological time step
    int    zb_init;

    int    exnerEq;           // exner equation to solve for bed evolution
                              // 1 = algebraic solve: (qb-qbe)/Ls
                              // 2 = numerical solve: div(qb) by finite volumes
                              // 3 = numerical solve: div(qb) by finite elements

    int    loadeq;            // equation for equilibrium bed load

    int    lsType;            // equation for non-equilibrium parameter Ls
    double minLs;             // minimum value for loading parameter Ls
    double factLs;            // factor: Ls >= d50 * lsfact
    double alfaLs;            // constant in Phillips-Sutherland law for Ls

    int    slope;             // equation to incorporate sloping effects
    double maxSlope;          // maximum slope of bed to allow erosion

    double alfaSlope;         // parameters to account for bed slope
    double betaSlope;
    double gammaSlope;
    double deltaSlope;

    enum
    {
      kSLOPE_Max=1,   kSLOPE_Shields=2,  kSLOPE_Qbe=4,
      kSLOPE_Tau=8,   kSLOPE_Angle=16,   kSLOPE_Dm=32
    };

    enum
    {
      kVanRijn=1,  kMeyerPM=2
    };

  public:
    SED()
    {
      isinit     = false;

      qbc        = NULL;
      PLs        = NULL;
      sx         = NULL;
      sy         = NULL;
      dzds       = NULL;
      dzmx       = NULL;
      dhds       = NULL;

      rhob       = 2500.0;

      M          = 0.0;
      tauc       = 0.0;
      taus       = 0.0;
      us         = 0.0;

      d50        = 0.0005;
      d90        = 0.0015;
      por        = 0.35;
      phir       = tan( PI/6.0 );        // tan(30�)

      loadeq     = 1;

      lsType     = 1;
      minLs      = 1.00;
      factLs     = 100.0;
      alfaLs     = 6000.0;

      slope      = 0;
      maxSlope   = tan( PI/4.0 );        // tan(45�)

      minQb      = 1.0e-9;
      convQb     = 1.0e-9;

      maxDz      = 0.01;

      zb_init    = false;

      exnerEq    = 3;

      alfaSlope  = 1.00;
      betaSlope  = 1.00;
      gammaSlope = 0.45;
      deltaSlope = 0.25;
    };

    int Getinit()  { return isinit; };

    //////////////////////////////////////////////////////////////////////////////////////
    // methods in module Sed.cpp

    void   Initialize( PROJECT* p );
    void   Detach();

    double Erode( NODE* nd );

    void   Bcon( PROJECT* p, MODEL* m );

    void   Capacity( PROJECT* p, MODEL* m, double* qbc, double* PLs,
                     double* dzds, double* dhds );
    void   Slope( PROJECT* p, MODEL* m, double* dzds, double* dzdn,
                  double* dzmx, double* dhds, double* dzdx, double* dzdy );
    void   Direction( PROJECT* p, MODEL* m, double* sx, double* sy );

    // -----------------------------------------------------------------------------------
    // functions to obtain transport capacity (equilibrium bed load)

    double Equilib( double Us, double H, double dzds, double dhds, PROJECT* p );
    double Meyerpm( double Us, double H, double dzds, double dhds, PROJECT* p );
    double Vanrijn( double Us, double H, double dzds, double dhds, PROJECT* p );

    // -----------------------------------------------------------------------------------
    // determine values for nonequilibrium parameter Ls
    double GetLs( double Us, double H, PROJECT* p );
    double Ls_vanrijn( double Us, double H, PROJECT* p );
    double Ls_phillips( double Us, double H, PROJECT* p );

    // -----------------------------------------------------------------------------------
    // determine the critical value of Shields parameter
    double Shields_crit( double Dst, double dzdx );

    // -----------------------------------------------------------------------------------
    double GetUtau( double Us, double H, double dzds, PROJECT* p );

    // -----------------------------------------------------------------------------------
    double Balance( ELEM* elem, int shape, double* etaQb, double V[kMaxNodes2D],
                    double* A, double* Vmax );
    double Balance( ELEM* elem, int shape, double* etaQb, double V[kMaxNodes2D] );
};

#endif

