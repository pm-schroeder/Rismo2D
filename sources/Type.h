// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// T Y P E
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Type.h   : definition file of the class.
// Type.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class computes the friction coefficient cf from different friction laws.
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
//  24.02.2000    sc    First implementation of C++ class.
//  03.03.2004    sc    Complete documentation of the source code (Type.h).
//  27.02.2005    sc    Some work on the import-method.
//  21.05.2006    sc    Sediment parameter removed.
//  01.08.2006    sc    Dispersion coefficient betaSf.
//  11.02.2008    sc    No further use of Default-Manning-Value nDef. If the specified
//                      roughness height k is less than flow depth, a Manning coefficient
//                      is computed internally for consistent transition.
//                      In later versions nDef will be removed later from the input files.
//  02.05.2008    sc    Implementation of van Velzen's assumption for flow over roughness.
//  09.05.2008    sc    Implementation of form roughness acoording to van Rijn (1993)
//  13.10.2012    sc    roughness of walls (log law) is no longer set by material zones
//                      thus rtype[2] and rcoef[2] are now rtype and rcoef
//                      wall roughnees may be applied by slip flow boundary conditions
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TYPE_INCL
#define TYPE_INCL

// ======================================================================================
// includes
// ======================================================================================

#include "Report.h"


// ======================================================================================

class ASCIIFILE;


class TYPE
{
  // ====================================================================================
  //                               V A R I A B L E N
  // ====================================================================================

  // ------------------------------------------------------------------------------------
  // private static variables (static = global inside of class context)
  private:
    static int   nType;        // number of different roughness zones
                               // = size of array pType
    static TYPE* pType;        // pointer to array of roughness zones

    static int   cWRav;
    static int   cWRmax;
    static int   cWRcount;
    static int   aNLav;
    static int   aNLmax;
    static int   aNLcount;
    static int   itErr;

    enum        { kMaxUndef = 100 };
    static int    nUndef;
    static int    list_of_undefined[kMaxUndef];

    static int    release;      // release of roughness-file

  public:
    enum        { kCOWI = 2,    // Colebrook-White's law
                  kNIKU = 3,    // Nikuradse's law
                  kCHEZ = 4,    // Chezy's law
                  kMANN = 5 };  // Manning's law

    static double _cb;          // fricition coefficient of bottom roughness
    static double _cp;          // fricition coefficient of non-submerged vegetation
    static double _cWR;         // drag coefficient of vegetation
    static double _aNL;         // wake length
    static double _Ust;         // ratio of flow velocity
    static double _Hst;         // ratio of flow depth

    static double _Hd;          // height and
    static double _Ld;          // length of dunes
    static double _Hr;          // height and
    static double _Lr;          // length of ripples

  // ------------------------------------------------------------------------------------
  // public static variables
  public:
    static TYPE  deflt;        // default values

  // ------------------------------------------------------------------------------------
  // private variables
  private:
    int    _id;                // ident number (incremented list)
    int    _no;                // material number

  // ------------------------------------------------------------------------------------
  // public variables
  public:
                           // --- $ROUGHNESS: bottom and wall roughness
    int    rtype;              // type of roughness coefficient
    double rcoef;              // value of roughness coefficient
    int    kslaw;              // which ks-parameter to use
                               // kslaw = 0:  ks = coef
                               // kslaw = 1:  ks = ksfact * d90
                               // kslaw = 2:  ks = ksfact * SED::d90
    double ksfact;             // see kslaw

                           // --- $ROUGHNESS: vegetation roughness
    double hp;                 // height of vegetation
    double dp[2];              // averaged diameter of roughness [m]
    double sp[2];              // averaged spacing of roughness [m]
                               // dp[0] / sp[0]      = non-submerged vegetation
                               // dp[1] / sp[1] / hp = submerged vegetation
    double d90;                // grain roughness ks = sediment parameter d90
    double d50;                // median of sediment
    int    form;               // account for form roughness, if form > 0
    double duneCoef;

                           // --- $TURBULENCE
    double vt;                 // constant eddy viscosity
    double st;                 // turbulent PRANDTL/SCHMIDT number
    double KDest;              // dimensionless depth-averaged diffusivity
    double lm;                 // mixing length constant

    double estx;               // dimensionless depth-averaged diffusivity
    double esty;

    double betaSf;             // dispersion coefficient beta (secondary flow)

  // ====================================================================================
  //                               M E T H O D E N
  // ====================================================================================

  public:
    // ----------------------------------------------------------------------------------
    // constructor: ...initialize all members of the class
    TYPE();

    // ----------------------------------------------------------------------------------
    // destructor:  ...delete allocated memory
    ~TYPE();

    // ----------------------------------------------------------------------------------
    // copy operator
    TYPE& operator =( const TYPE& );

    // ----------------------------------------------------------------------------------
    // import method to read the roughness map
    static void Import( const char* filename );
    static void Import_00000( ASCIIFILE* tbl );
    static void Import_30900( ASCIIFILE* tbl );
    static void Import_40000( ASCIIFILE* tbl );

    // ----------------------------------------------------------------------------------
    // static methods for screen output
/*
    static void err( const char* format, ... )
    {
      prnt( "\n=== ERROR: " );
      prnt( format );
      prnt( "\n" );

#     ifdef _MPI_
      MPI_Finalize();
#     endif

      exit( 0 );
    }

    static void prnt( const char* format, ... )
    {
      char text[200];

      va_list argPtr;
      va_start( argPtr, format );

      vsprintf( text, format, argPtr );
      REPORT::rpt.Output( text, 5 );

      va_end( argPtr );
    }
*/

    // ----------------------------------------------------------------------------------
    // static method to initialize the class
    // <n> is the number of roughness/material zones
    static void Init( int n )
    {
      nType = n;
      pType = new TYPE [nType];

      if( !pType )
        REPORT::rpt.Error( kMemoryFault, "can not allocate memory (TYPE::init - 1)" );

      for( int i=0; i<nType; i++ )
      {
        pType[i]._id = i + 1;
        pType[i]._no = 0;

        pType[i] = deflt;
      }
    }

    // ----------------------------------------------------------------------------------
    // static method Get() returns the number of roughness zones

    static int Get()
    {
      return nType;
    }

    // ----------------------------------------------------------------------------------
    // static method Getid() returns a pointer to the roughness zone with index <i>

    static TYPE* Getid( int i )
    {
      if( i <= 0 )      return &deflt;
      if( i <= nType )  return &pType[i-1];

      REPORT::rpt.Error( kParameterFault,
                         "element type %d out of range %s", i, "(TYPE::getid #1)" );
      return NULL;
    }

    // ----------------------------------------------------------------------------------
    // static method Getno() returns a pointer to the roughness zone <no>

    static TYPE* Getno( int no )
    {
      for( int i=0; i<nType; i++ )
      {
        if( pType[i]._no == no )  return &pType[i];
      }

      for( int i=0; i<nUndef; i++ )
      {
        if( list_of_undefined[i] == no )  return &deflt;
      }

      if( nUndef < kMaxUndef )
      {
        REPORT::rpt.Warning( kParameterFault,
                             "undefined element type %d: %s (%s)",
                             no, "using default values", "TYPE::Getno #2" );

        list_of_undefined[nUndef] = no;
        nUndef++;
      }

      return &deflt;
    }

    // ----------------------------------------------------------------------------------
    // Setno()

    void Setno( int no )
    {
      _no = no;
    }

    // ----------------------------------------------------------------------------------
    // return the class members <_id> und <_no>
    int id()
    {
      return _id;
    }

    int no( int type )
    {
      if( _no <= 0 )  return abs( type );
      else            return _no;
    }

    // ----------------------------------------------------------------------------------
    // Module Friction.cpp
    // ----------------------------------------------------------------------------------
    // The method bottom() computes the total friction from the superposition of bottom
    // roughness and non-submerged vegetation
    //
    // Parameter:  int    rtype = roughness type
    // ---------   double rc    = roughness parameter depending on rType
    //             double Ur    = resultant velocity
    //             double h     = flow depth
    //             double ka    = von Karman's constant ( = 0.41 )
    //             double vk    = kinematic viscosity
    //             double dw    = wall distance
    //             double g     = gravity acceleration ( = 9.81 )
    double bottom( double Us, double h, double ka, double vk, double g,
                   double rho, double rhob, double d50, double d90 );
    double friction( int rtype, double rc, double Us, double h,
                     double ka, double vk, double g );
    void Dune( double Us, double H,
               double ka, double vk, double g, double rho,
               double rhob, double d50, double d90,
               double* Hd, double* Ld, double* Hr, double* Lr );

    // ----------------------------------------------------------------------------------
    // Module Lindner.cpp
    // ----------------------------------------------------------------------------------
    // The method lindner() computes the flow resistance of
    // non-submerged vegetation
    //
    // Parameter:  double dp = diameter of vegetation
    //             double sp = spacing of vegetation
    //             double va = velocity
    //             double ha = flow depth
    //             double cf = friction coefficient for bottom roughness
    //             double vk = kinemtic viscosity
    //             double g  = gravity acceleration
    double lindner( double dp, double sp, double va, double ha,
                    double cf, double vk, double g );
    void   statisLindner();
    double dragCoeff( double, double, double );
    int    cubeEquation( double, double, double, double, double [3] );

    // ----------------------------------------------------------------------------------
    // Module Velzen.cpp
    // ----------------------------------------------------------------------------------
    // The method velzen() computes the flow resistance of
    // submerged vegetation in case of h > hp and of
    // non-submerged vegetation according to lindner() in case of h < hp
    //
    // Parameter:  double hp  = height of vegetation
    //             double dp = diameter of vegetation
    //             double sp = spacing of vegetation
    //             double va = velocity
    //             double ha = flow depth
    //             double cf = friction coefficient for bottom roughness
    //             double ka = von Karman's constant ( = 0.41 )
    //             double vk = kinemtic viscosity
    //             double g  = gravity acceleration
    double velzen( double hp, double dp, double sp, double va, double ha,
                   double cf, double ka, double vk, double g );
};

#endif
