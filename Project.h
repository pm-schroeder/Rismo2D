// ======================================================================================
//                                   P R O J E C T
// ======================================================================================
// This class implements the project class.
// ======================================================================================
//
// Copyright (C) 1992-2013  by  P.M. SCHROEDER
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
// 01.01.1998     sc     first implementation / first concept
// 13.10.2011     sc     cycles "..._dt" + "..._tr" removed, stationary flow computation
//                       will be established by the key $STATIONARY im the tmiestep-file
// 28.10.2011     sc     cleaning up keys structure (RISKEY)

// ======================================================================================

#include <string.h>

#ifndef PROJECT_INCL
#define PROJECT_INCL

#include "Datkey.h"
#include "Defs.h"
#include "Eqs.h"
#include "EqsSL2D.h"
#include "EqsBL2D.h"
#include "EqsDisp.h"
#include "EqsDz.h"
#include "EqsK2D.h"
#include "EqsD2D.h"
#include "EqsKD2D.h"
#include "EqsKL2D.h"
#include "EqsPPE2D.h"
#include "EqsUVS2D.h"
#include "EqsUVS2D_TM.h"
#include "EqsUVS2D_AI.h"
#include "EqsUVS2D_LV.h"
#include "EqsUVS2D_TMAI.h"
#include "Sed.h"
#include "Solver.h"
#include "Scale.h"
#include "Subdom.h"
#include "Timeint.h"

class TYPE;
class MODEL;
class SECTION;
class ASCIIFILE;
class STATIST;


// cycles ------------------------------------------------------------------------------------------

#define kUVSCyc           1  // 2D: UVS quadratic/linear, steady+unsteady flow
#define kUVS_TMCyc        2  // 2D: UVS quadratic/linear, unsteady flow, time prediction

#define kUVS_LVCyc        4  // 2D: UVS linear/constant, steady flow

#define kDispCurv2D      20  // compute dispersion from streamline curvature

#define kSurfaceCyc      30  // initialize water surface elevation
#define kSmoothS         31  // smoothing of free surface
#define kSurfaceToVol    35  // compute surface in volumes from nodes

#define kDivCyc          40  // zero-divergence cycle, quadratic/linear

#define kQbCyc           50  // compute sediment transport rates

#define kSLCyc           51  // suspended load cycle
#define kBLCyc           52  // bed load cycle
#define kBSLCyc          53  // bed and suspended load cycle

#define kDiffBLCyc       55  // difference in bed load used for bottom evolution

//#define kKDCyc_Pkv       60  // KD initialization, bottom friction
//#define kKDCyc_Pr_Pkv    61  // KD initialization, bottom friction + velocity gradients
//#define kKDCyc_Pr        62  // Prandtls mixing model with prescribed lm
//#define kKDCyc_LES       63  // Prandtls mixing model, grid depending lm (LES)

#define kKLCyc           71  // one equation cycle

#define kKDInitCyc       80  // KD initialization, bottom friction
#define kKDCyc           81  // k-epsilon cycle, steady flow
#define kKD_LCyc         84  // k-epsilon cycle, steady flow
#define kKD_QCyc         87  // k-epsilon cycle, steady flow (quartered elements)

#define kDryRewet        90  // dry/rewet cycle
#define kReOrderCyc      95  // reorder cycle
#define kScaledOutputCyc 98  // write scaled output files
#define kOutputCyc       99  // write output files

#define kKCyc           181  // k cycle, steady flow
#define kKCyc_dt        182  // k cycle
#define kKCyc_tr        183  // k cycle, steady flow, time relaxed

#define kDCyc           281  // epsilon cycle, steady flow
#define kDCyc_dt        282  // epsilon cycle
#define kDCyc_tr        283  // epsilon cycle, steady flow, time relaxed


// file names ----------------------------------------------------------------------------

struct NAME
{
  enum { kLength = 500};

  char inputFile[kLength];       // name of input file
  char materialFile[kLength];
  char reportFile[kLength];

  int  ascii_initial;
  int  ascii_prevTime;
  int  ascii_restart;
  int  ascii_statist;

  char timeStepFile[kLength];
  char regionFile[kLength];
  char controlFile[kLength];
  char initialFile[kLength];
  char subdomFile[kLength];
  char sectionFile[kLength];
  char geometryFile[kLength];
  char restartFile[kLength];
  char rgAvsFile[kLength];
  char ctAvsFile[kLength];
  char redAvsFile[kLength];
  char wetAvsFile[kLength];
  char inputStatistFile[kLength];
  char outputStatistFile[kLength];
  char rgAvsStatistFile[kLength];
  char optOutputFile[kLength];

  NAME()
  {
    ascii_initial  = 1;
    ascii_prevTime = 1;
    ascii_restart  = 1;
    ascii_statist  = 1;

    inputFile[0]         = '\0';
    timeStepFile[0]      = '\0';
    materialFile[0]      = '\0';
    reportFile[0]        = '\0';
    sectionFile[0]       = '\0';
    regionFile[0]        = '\0';
    controlFile[0]       = '\0';
    initialFile[0]       = '\0';
    subdomFile[0]        = '\0';
    geometryFile[0]      = '\0';
    restartFile[0]       = '\0';
    rgAvsFile[0]         = '\0';
    ctAvsFile[0]         = '\0';
    redAvsFile[0]        = '\0';
    wetAvsFile[0]        = '\0';
    inputStatistFile[0]  = '\0';
    outputStatistFile[0] = '\0';
    rgAvsStatistFile[0]  = '\0';
    optOutputFile[0]     = '\0';
  };

  int Length()
  {
    return kLength;
  };
};


// k-epsilon constants -------------------------------------------------------------------

struct KDCONST
{
  double cm;                 // standard k-epsilon constants
  double cd;                 //   cm  = 0.09;  cd  = 1.0
  double sK;                 //   c1D = 1.43;  c2D = 1.92
  double sD;                 // depending on von Karman's constant:
  double c1D;                // kappa = 0.43: sk = 1.0; sD = 1.3
  double c2D;                //         0.41: sk = 1.2; sD = 1.2
};


// class PROJECT declaration -------------------------------------------------------------

#define kIndInflow           0     // subscripts for KDBcon array
#define kIndOutflow          1
#define kIndSide             2

#define kErr_no_error        0     // error levels
#define kErr_some_errors     1
#define kErr_interrupt       2

#define kErr_no_conv_nr      4     // kind of errors
#define kErr_no_conv_cg      8


struct VALIST                      // list of variables
{
  int         id;
  int         dim;                 // internal dimension        e.g. UV : dim = 2
  int         vec;                 // vector length in ucd-file           vec = 3
  const char* name;
  const char* unit;
};

struct TMSER
{
  int    ntm;
  int*   tmlist;                        // array of time steps for output
  int    first;
  int    last;
  int    step;

  int    vdata;
  int    vcomp;
  int   *voutlist;                      // list of node variables for output

  char   filename[NAME::kLength];       // file name for output

  TMSER()
  {
    ntm  = 0;
    vdata = 0;
    vcomp = 0;

    filename[0] = '\0';
  }

//  TMSER( const TMSER &t )
//  {
//
//  }

  ~TMSER()
  {
    if( ntm ) delete[] tmlist;
  }

  TMSER& operator =( const TMSER &t )
  {
    if( this == &t ) return *this;
    if( this->ntm ) delete[] this->tmlist;

    this->ntm   = t.ntm;

    this->first = t.first;
    this->last  = t.last;
    this->step  = t.step;

    this->vdata = t.vdata;
    this->vcomp = t.vcomp;

    this->tmlist = new int[this->ntm];
    for( int i=0; i<this->ntm; i++ ) this->tmlist[i] = t.tmlist[i];

    this->voutlist = new int[this->vcomp];
    for( int i=0; i<this->vcomp; i++ ) this->voutlist[i] = t.voutlist[i];

    strcpy( this->filename, t.filename );

    return *this;
  }
};


class PROJECT
{
  private:
    int     nkey;
    DATKEY* datkey;

    enum RISKEY
    {
      kVERSION=0,        kMACRO,            kTITLE,            kCOMMENT,
      kINPUTFILE,        kREPORTLEVEL,      kREPORTFILE,       kREPORTTIME,

      kTIMESTEPFILE,     kMATERIALFILE,     kREGIONFILE,       kCONTROLFILE,
      kASC_INITFILE,     kBIN_INITFILE,     kSTA_INITFILE,     kSUBDOMFILE,
      kOUTPUTPATH,

      kASC_RESTFILE,     kBIN_RESTFILE,     kSTA_RESTFILE,
      kRG_UCDFILE,       kCT_UCDFILE,       kCN_UCDFILE,       kWN_UCDFILE,
      kST_UCDFILE,       kTIMESERIES,       kOUTPUTVARS,

      kSOLVER,           kCONVERGENCE,      kMIN_VT,           kMIN_VT_XX,
      kMIN_VT_YY,        kMIN_VT_KD,        kMIN_VT_C,         kMIN_K,
      kMIN_D,            kMIN_C,            kMAX_Us,           kRELAX,
      kDRYREW,           kSMOOTH,           kKDBOUNDARY,       kGPDEGREE,
      kDELTAZLIMIT,      kTEMPERATURE,      kVISCOSITY,        kDENSITY,
      kGRAVITY,          kVON_KARMAN,       kEARTH_ROTATION,   kLATITUDE,
      kKDCONST,

      kMUE_SF,           kMINU_SF,          kMAXTAN_SF,

      kSCALE,

      kSED_RHOB,         kSED_M,            kSED_TAUC,         kSED_TAUS,
      kSED_US,           kSED_D50,          kSED_D90,          kSED_POR,
      kSED_PHIR,         kSED_LOADEQ,       kSED_LS,           kSED_SLOPE,
      kSED_MINQB,        kSED_MAXDZ,        kSED_EXNEREQ,      kSED_ZB_INIT,

      // deprecated keys
      kMINMAX,

      kSZ_RISKEY
    };

    int     iTM;              // number of current time step (loop variable)

  public:
    enum VARS
    {
      kUV=0,    kS,      kK,      kD,     kC,    kCB,    kQB,   kQBE,   kLS,      kSEDUV,
      kDUVDT,   kDSDT,   kZ,      kH,     kUS,   kCF,    kRC,   kD90,   kD50,     kKD,
      kHR,      kHD,     kHP,     kDP,    kSP,   kTAU,   kMAN,  kVT,    kEst,     kEXX,
      kEYY,     kDUU,    kDUV,    kDVV,   kVSEC, kUVBOT, kDZ,   kRE,    kFR,      kPE,
      kCU,      kPHI,    kROT,    kCURV,  kDZDS, kDZDN,  kDZMX, kDHDS,  kMEANUV,  kMEANS,
      kMEANUS,  kMEANH,  kMEANVT, kVARU,  kVARV, kVARUV, kKINE, kSDEVH, kVARVT,   kKINER,
      kFLDRATE, kKINRATIO, kMAXUS, kMAXTAU, kSZ_VARS
    };

    int     nval;
    VALIST* valist;

    int     vpdata;
    int     vpcomp;
    int*    vpoutlist;        // list of node variables for output in AVS-UCS-Files

    int     vedata;
    int     vecomp;
    int*    veoutlist;        // list of element variables for output in AVS-UCS-Files

    int     ntmser;
    TMSER  *tmser;

public:
    static int release;

    enum     { kMaxMacro = 20 };
    int        nMacro;
    char       macro_key[kMaxMacro][100];
    char       macro_value[kMaxMacro][100];

    SUBDOM     subdom;                  // class for domain decomposition
    char       outputPath[500];

    NAME       name;

    MODEL*     M2D;

    TIMEINT    timeint;

    SCALE      scale;

    int        nSection;
    SECTION*   section;

    int        elemKind;                // flag to include boundary elements

    // ----------------------------------- differential equation systems -----------------
    EQS*           eqs;

    EQS_SL2D       eqs_sl2d;            // equation system for suspended transport
    EQS_BL2D       eqs_bl2d;            // equation system for bed load transport
    EQS_DZ         eqs_dz;              // equation system for bottom evolution
    EQS_DISP       eqs_disp;            // equation system for dispersion coefficients
    EQS_K2D        eqs_k2d;             // equation system for k
    EQS_D2D        eqs_d2d;             // equation system for epsilon
    EQS_KD2D       eqs_kd2d;            // equation system for k-epsilon
    EQS_KL2D       eqs_kl2d;            // equation system for k
    EQS_PPE2D      eqs_ppe2d;           // equation system for divergence free flow field
    EQS_UVS2D      eqs_uvs2d;           // equation system for shallow water flow
    EQS_UVS2D_TM   eqs_uvs2d_tm;        // equation system for shallow water flow
    EQS_UVS2D_AI   eqs_uvs2d_ai;        // equation system for shallow water flow
    EQS_UVS2D_LV   eqs_uvs2d_lv;
    EQS_UVS2D_TMAI eqs_uvs2d_tmai;      // equation system for shallow water flow

    int            errLevel;

    // ----------------------------------- arrays ----------------------------------------
    double*  lmm;                       // lumped mass matrix

    // ----------------------------------- constants -------------------------------------
    double   g;                         // gravity
    double   earthVel;                  // earth rotation velocity
    double   latitude;                  // latitude of coordinate y = 0
    double   celsius;
    double   vk;                        // kinematic viscosity
    double   rho;                       // density of water
    double   kappa;                     // von Karman's constant

    KDCONST  KD;                        // k-epsilon constants

    double   hmin;                      // minimum value for flow depth

    double   mueSf;                     // dispersion coefficient for adaption length
    double   minUSf;                    // min for mean flow velocity
    double   maxTanSf;                  // max for ratio secondary flow / mean flow


    // ----------------------------------- iteration/convergence control -----------------
    double   convUV;                    // convergence parameters
    double   convS;
    double   convKD;

    int      relaxMethod;               // relaxation method (1 or 2)
    double   relaxMin,                  // minimum relaxation parameter  (< 1.0)
             relaxMax,                  // maximum relaxation parameter  (= 1.0)
             maxDeltaUV,                // maximum allowed change of velocity
             maxDeltaS,                 // maximum allowed change of water elevation
             maxDeltaKD;                // maximum allowed changes of K and D

    int      smoothPassesBC;            // number of smoothing passes for bc.
    int      smoothPassesKD;            // number of smoothing passes for KD
    int      smoothPassesVT;            // number of smoothing passes for vt

    // -------------------------------- control parameter --------------------------------
    int      KDBcon[3];                 // boundary model for K and D
    int      GPdeg;                     // degree of GAUSS point integration

    unsigned
    int      fix[kSimDF];               // flags to fix equations

    double   dep_minVt;                 // minimum vt
    double   dep_minVtxx;               // minimum vt in flow direction (Elder)
    double   dep_minVtyy;               // minimum vt in lateral direction (Elder)
    double   minVtKD;                   // minimum vt in KD-equations
    double   minVtC;                    // minimum vt in sediment transport

    double   minK;                      // minimum K
    double   minD;                      // minimum D
    double   minC;                      // minimum concentration

    double   maxUs;                     // maximum flow velocity

    // -------------------------------- cycle parameter ----------------------------------
    int      coefsFunc;
    int      theCycle;                  // actual cycle
    int      actualStat;                // stationary computation
    int      actualTurb;                // turbulence model
    int      actualDisp;                // dispersion model
    int      actualCycit;               // actual number of iterations per cycle
    SOLVER*  actualSolver;              // actual solver type

    int      flowFlag;                  // flag to recompute flow field

    // -------------------------------- sediment parameter -------------------------------
    SED      sed;

    int      statistics;
    STATIST* statist;

    // -------------------------------- gauges controlling outlet water level -----------------------
    int             ngauge;             // number of gauges (over all subdomains)
    int            *gauge;              // name of nodes which are a gauge
    double         *gaugeS;             // current water elevation at gauge
    double         *gaugeSo;            // targeted water elevation at gauge

    int             ngct;               // number of gauge controlled nodes on outlet
    BCVAL::GAUGECT *gct;                // list of gauge controlled nodes on outlet

  public:
    PROJECT();
    ~PROJECT();

    int     GetTimeStep() { return iTM; }

    void    ReplaceMacro( char* name );
    void    ReplaceMacro(char* src, char *dst, int tmstep=0, int pid=0 );
    bool    ReplaceKey( char* src, char *dst, char* key, char* val );
    void    ReplaceAllKeys( char* src, char *dst, char* key, char* val );
    void    SetPartKey( char *filename );

    void    Input();
    void    Input_00000( ASCIIFILE* file );
    void    Input_30900( ASCIIFILE* file );

    int Cube( double a, double b, double c, double d, double x[3] )
    {
      int    real;

      double ba = b / a / 3.0;
      double ca = c / a;

      double P  = ca/3.0 - ba*ba;
      double Q  = ba*ba*ba - ba*ca/2.0 + d/a/2.0;

      double Q2P3 = Q*Q + P*P*P;

      if( Q2P3 > 0.0 )
      {
        real = 1;

        double tmp, U, V;
        double expo = 1.0/3.0;

        tmp  = -Q + sqrt(Q2P3);
        if( tmp < 0.0 )  U = -pow(-tmp,expo );
        else             U =  pow( tmp,expo );

        tmp  = -Q - sqrt(Q2P3);
        if( tmp < 0.0 )  V = -pow(-tmp,expo );
        else             V =  pow( tmp,expo );

        x[0] = (U + V) - ba;
      }
      else
      {
        real = 3;

        double phi;
        double tmp = -Q / pow( -P, 1.5 );

        if( tmp >= 1.0 )        phi = 0.0;
        else if( tmp <= -1.0 )  phi = PI;
        else                    phi = acos( tmp );

        x[0] = 2.0 * sqrt(-P) * cos(phi/3.0)           -  ba;
        x[1] = 2.0 * sqrt(-P) * cos((phi+2.0*PI)/3.0)  -  ba;
        x[2] = 2.0 * sqrt(-P) * cos((phi+4.0*PI)/3.0)  -  ba;
      }

      return real;
    };

    // Compute.cpp -----------------------------------------------------------------------
    void    Compute();

    // Cycle.cpp -------------------------------------------------------------------------
    int     NextCycle( BCONSET*, int );
    void    PrintTheCycle( int );
};

#endif

