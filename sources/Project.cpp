// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class PROJECT
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
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
// /////////////////////////////////////////////////////////////////////////////////////////////////

#include "Defs.h"
#include "Asciifile.h"
#include "Report.h"
#include "Scale.h"
#include "Type.h"
#include "Shape.h"
#include "Node.h"
#include "Elem.h"
#include "Grid.h"
#include "Model.h"
#include "Section.h"
#include "Statist.h"
#include "Solver.h"
#include "Front.h"
#include "Frontm.h"
#include "Bicgstab.h"
#include "P_bcgstabd.h"
#include "P_fgmresd.h"

#include "Project.h"


double PROJECT::release = 40521.001;


PROJECT::PROJECT()
{
  lmm = NULL;

  errLevel = kErr_no_error;

  nSection = 0;

  smoothPassesBC = 2;
  smoothPassesKD = 0;
  smoothPassesVT = 0;

  convUV = 1.0e-3;
  convS  = 1.0e-4;
  convKD = 1.0e-4;

  dep_minVt   = 0.0;
  dep_minVtxx = 0.0;
  dep_minVtyy = 0.0;

  minVtKD = 0.0;
  minVtC  = 0.0;
  minK    = 0.0;
  minD    = 0.0;
  minC    = 0.0;
  maxUs   = 100.0;

  relaxMethod = 3;
  relaxMin    = 0.001;
  relaxMax    = 0.1;

  maxDeltaUV  = 0.1;
  maxDeltaS   = 0.1;
  maxDeltaKD  = 0.1;

  mueSf       = 1.0;
  maxTanSf    = 0.5;
  minUSf      = 0.01;

  KDBcon[kIndInflow]  = 2;
  KDBcon[kIndOutflow] = 1;
  KDBcon[kIndSide]    = 1;

  GPdeg = 5;

  // the kinematic viscosity is a function of temperature and pressure;
  // in case of water (10 degree Celsius) it is approximately
  //
  // vk = 1.78e-6 / (1.0 + 0.0337*T + 2.2e-4*T*T)

  celsius   = 10.0;
  vk        = 1.78e-6 / (1.0 + 0.0337*celsius + 0.00022*celsius*celsius);
  rho       = 1000.0;
  g         = 9.81;
  kappa     = 0.41;
  earthVel  = 0.0001458;
  latitude  = 0.0;

  KD.cm     = 0.09;
  KD.cd     = 1.00;
  KD.sK     = 1.20;
  KD.sD     = 1.20;
  KD.c1D    = 1.43;
  KD.c2D    = 1.92;

  // -------------------------------------------------------------------------------------

  static DATKEY dk[] =
  {
    kVERSION,         "VERSION",            // 00
    kMACRO,           "MACRO",              // 01
    kTITLE,           "TITLE",              // 02
    kCOMMENT,         "COMMENT",            // 03

    kINPUTFILE,       "INPUTFILE",          // 04

    kREPORTLEVEL,     "REPORTLEVEL",        // 05
    kREPORTFILE,      "REPORTFILE",         // 06
    kREPORTTIME,      "REPORTTIME",         // 07

    kTIMESTEPFILE,    "TIMESTEPFILE",       // 08
    kMATERIALFILE,    "MATERIALFILE",       // 09
    kREGIONFILE,      "REGIONFILE",         // 10
    kCONTROLFILE,     "CONTROLFILE",        // 11
    kASC_INITFILE,    "ASC_INITFILE",       // 12
    kBIN_INITFILE,    "BIN_INITFILE",       // 13
    kSTA_INITFILE,    "STA_INITFILE",       // 14
    kSUBDOMFILE,      "SUBDOMFILE",         // 15
    kOUTPUTPATH,      "OUTPUTPATH",         // 16

    kASC_RESTFILE,    "ASC_RESTFILE",       // 17
    kBIN_RESTFILE,    "BIN_RESTFILE",       // 18
    kSTA_RESTFILE,    "STA_RESTFILE",       // 19
    kRG_UCDFILE,      "RG_UCDFILE",         // 20
    kCT_UCDFILE,      "CT_UCDFILE",         // 21
    kCN_UCDFILE,      "CN_UCDFILE",         // 22
    kWN_UCDFILE,      "WN_UCDFILE",         // 23
    kST_UCDFILE,      "ST_UCDFILE",         // 24
    kTIMESERIES,      "TIMESERIES",         // 25

    kOUTPUTVARS,      "OUTPUTVARS",         // 26

    kSOLVER,          "SOLVER",             // 27
    kCONVERGENCE,     "CONVERGENCE",        // 28
    kMIN_VT,          "MIN_VT",             // 29
    kMIN_VT_XX,       "MIN_VT_XX",          // 30
    kMIN_VT_YY,       "MIN_VT_YY",          // 31
    kMIN_VT_KD,       "MIN_VT_KD",          // 32
    kMIN_VT_C,        "MIN_VT_C",           // 33
    kMIN_K,           "MIN_K",              // 34
    kMIN_D,           "MIN_D",              // 35
    kMIN_C,           "MIN_C",              // 36
    kMAX_Us,          "MAX_Us",             // 37
    kRELAX,           "RELAX",              // 38
    kDRYREW,          "DRYREW",             // 39
    kSMOOTH,          "SMOOTH",             // 40
    kKDBOUNDARY,      "KDBOUNDARY",         // 41
    kGPDEGREE,        "GPDEGREE",           // 42
    kDELTAZLIMIT,     "DELTAZLIMIT",        // 43
    kTEMPERATURE,     "TEMPERATURE",        // 44
    kVISCOSITY,       "VISCOSITY",          // 45
    kDENSITY,         "DENSITY",            // 46
    kGRAVITY,         "GRAVITY",            // 47
    kVON_KARMAN,      "VON_KARMAN",         // 48
    kEARTH_ROTATION,  "EARTH_ROTATION",     // 49
    kLATITUDE,        "LATITUDE",           // 50
    kKDCONST,         "KDCONST",            // 51

    kMUE_SF,          "MUE_SF",             // 52
    kMINU_SF,         "MINU_SF",            // 53
    kMAXTAN_SF,       "MAXTAN_SF",          // 54

    kSCALE,           "SCALE",              // 55

    kSED_RHOB,        "SED_RHOB",           // 56
    kSED_M,           "SED_M",              // 57
    kSED_TAUC,        "SED_TAUC",           // 58
    kSED_TAUS,        "SED_TAUS",           // 59
    kSED_US,          "SED_US",             // 60
    kSED_D50,         "SED_D50",            // 61
    kSED_D90,         "SED_D90",            // 62
    kSED_POR,         "SED_POR",            // 63
    kSED_PHIR,        "SED_PHIR",           // 64
    kSED_LOADEQ,      "SED_LOADEQ",         // 65
    kSED_LS,          "SED_LS",             // 66
    kSED_SLOPE,       "SED_SLOPE",          // 67
    kSED_MINQB,       "SED_MINQB",          // 68
    kSED_MAXDZ,       "SED_MAXDZ",          // 69
    kSED_EXNEREQ,     "SED_EXNEREQ",        // 70
    kSED_ZB_INIT,     "SED_ZB_INIT",        // 71

    // depreciated keys (recognized for compatibility reasons)
    kMINMAX,          "MINMAX",             // 72

    // key with changed names (recognized for compatibility reasons)
    kASC_INITFILE,    "ASC_INIFILE",        // 73
    kBIN_INITFILE,    "BIN_INIFILE",        // 74
    kSTA_INITFILE,    "STA_INIFILE",        // 75
    kASC_RESTFILE,    "ASC_RESTARTFILE",    // 76
    kBIN_RESTFILE,    "BIN_RESTARTFILE",    // 77
    kSTA_RESTFILE,    "STA_OUTFILE",        // 78
    kCN_UCDFILE,      "RED_UCDFILE",        // 79
    kWN_UCDFILE,      "WET_UCDFILE",        // 80
    kST_UCDFILE,      "STA_UCDFILE",        // 81

    kRG_UCDFILE,      "GEO_UCDFILE",        // 82

    kOUTPUTPATH,      "SUBDOMPATH",         // 83

    kREPORTLEVEL,     "REPPORTLEVEL",       // 84
    kREPORTFILE,      "REPPORTFILE"         // 85
 };

  nkey   = kSZ_RISKEY + 13;
  datkey = dk;

  //---------------------------------------------------------------------------------------

  nval = kSZ_VARS;                          // kSZ_VARS = ...

  valist    = new VALIST[nval];
  vpoutlist = new int[nval];
  veoutlist = new int[nval];

  if( !valist || !vpoutlist || !veoutlist )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (MODEL::MODEL #1)" );

  VALIST* vl = valist;
  vl->id = kUV;      vl->name = "UV";      vl->unit = "m/s";    vl->dim = 2;  vl->vec = 3;  vl++;  //  0
  vl->id = kS;       vl->name = "S";       vl->unit = "mNN";    vl->dim = 1;  vl->vec = 1;  vl++;  //  1
  vl->id = kK;       vl->name = "K";       vl->unit = "m2/s2";  vl->dim = 1;  vl->vec = 1;  vl++;  //  2
  vl->id = kD;       vl->name = "D";       vl->unit = "m2/s3";  vl->dim = 1;  vl->vec = 1;  vl++;  //  3
  vl->id = kC;       vl->name = "C";       vl->unit = "kg/m3";  vl->dim = 1;  vl->vec = 1;  vl++;  //  4
  vl->id = kCB;      vl->name = "Cb";      vl->unit = "kg/m2";  vl->dim = 1;  vl->vec = 1;  vl++;  //  5
  vl->id = kQB;      vl->name = "qb";      vl->unit = "m3/m/s"; vl->dim = 1;  vl->vec = 1;  vl++;  //  6
  vl->id = kQBE;     vl->name = "qbe";     vl->unit = "m3/m/s"; vl->dim = 1;  vl->vec = 1;  vl++;  //  7
  vl->id = kLS;      vl->name = "Ls";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  //  8
  vl->id = kSEDUV;   vl->name = "SedUV";   vl->unit = "m/s";    vl->dim = 2;  vl->vec = 3;  vl++;  //  9
  vl->id = kDUVDT;   vl->name = "dUVdt";   vl->unit = "m/s2";   vl->dim = 2;  vl->vec = 3;  vl++;  // 10
  vl->id = kDSDT;    vl->name = "dSdt";    vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 11
  vl->id = kZ;       vl->name = "Zb";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 12
  vl->id = kH;       vl->name = "H";       vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 13
  vl->id = kUS;      vl->name = "Us";      vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 14
  vl->id = kCF;      vl->name = "cf";      vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 15
  vl->id = kRC;      vl->name = "rc";      vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 16
  vl->id = kD90;     vl->name = "d90";     vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 17
  vl->id = kD50;     vl->name = "d50";     vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 18
  vl->id = kKD;      vl->name = "kd";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 19
  vl->id = kHR;      vl->name = "Hr";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 20
  vl->id = kHD;      vl->name = "Hd";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 21
  vl->id = kHP;      vl->name = "hp";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 22
  vl->id = kDP;      vl->name = "dp";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 23
  vl->id = kSP;      vl->name = "sp";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 24
  vl->id = kTAU;     vl->name = "tau";     vl->unit = "N/m2";   vl->dim = 1;  vl->vec = 1;  vl++;  // 25
  vl->id = kMAN;     vl->name = "man";     vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 26
  vl->id = kVT;      vl->name = "vt";      vl->unit = "m2/s";   vl->dim = 1;  vl->vec = 1;  vl++;  // 27
  vl->id = kEst;     vl->name = "Est";     vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 28
  vl->id = kEXX;     vl->name = "Exx";     vl->unit = "m2/s";   vl->dim = 1;  vl->vec = 1;  vl++;  // 29
  vl->id = kEYY;     vl->name = "Eyy";     vl->unit = "m2/s";   vl->dim = 1;  vl->vec = 1;  vl++;  // 30
  vl->id = kDUU;     vl->name = "Duu";     vl->unit = "m2/s";   vl->dim = 1;  vl->vec = 1;  vl++;  // 31
  vl->id = kDUV;     vl->name = "Duv";     vl->unit = "m2/s";   vl->dim = 1;  vl->vec = 1;  vl++;  // 32
  vl->id = kDVV;     vl->name = "Dvv";     vl->unit = "m2/s";   vl->dim = 1;  vl->vec = 1;  vl++;  // 33
  vl->id = kVSEC;    vl->name = "Vsec";    vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 34
  vl->id = kUVBOT;   vl->name = "UVbot";   vl->unit = "m/s";    vl->dim = 2;  vl->vec = 3;  vl++;  // 35
  vl->id = kDZ;      vl->name = "dz";      vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 36
  vl->id = kRE;      vl->name = "Re";      vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 37
  vl->id = kFR;      vl->name = "Fr";      vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 38
  vl->id = kPE;      vl->name = "Pe";      vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 39
  vl->id = kCU;      vl->name = "Cu";      vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 40
  vl->id = kPHI;     vl->name = "phi";     vl->unit = "1/s2";   vl->dim = 1;  vl->vec = 1;  vl++;  // 41
  vl->id = kROT;     vl->name = "rot";     vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 42
  vl->id = kCURV;    vl->name = "curv";    vl->unit = "1/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 43
  vl->id = kDZDS;    vl->name = "dzds";    vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 44
  vl->id = kDZDN;    vl->name = "dzdn";    vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 45
  vl->id = kDZMX;    vl->name = "dzmx";    vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 46
  vl->id = kDHDS;    vl->name = "dhds";    vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 47
  vl->id = kMEANUV;  vl->name = "meanUV";  vl->unit = "m/s";    vl->dim = 2;  vl->vec = 3;  vl++;  // 48
  vl->id = kMEANS;   vl->name = "meanS";   vl->unit = "mNN";    vl->dim = 1;  vl->vec = 1;  vl++;  // 49
  vl->id = kMEANUS;  vl->name = "meanUs";  vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 50
  vl->id = kMEANH;   vl->name = "meanH";   vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 51
  vl->id = kMEANVT;  vl->name = "meanVt";  vl->unit = "m2/s";   vl->dim = 1;  vl->vec = 1;  vl++;  // 52
  vl->id = kVARU;    vl->name = "varU";    vl->unit = "m2/s2";  vl->dim = 1;  vl->vec = 1;  vl++;  // 53
  vl->id = kVARV;    vl->name = "varV";    vl->unit = "m2/s2";  vl->dim = 1;  vl->vec = 1;  vl++;  // 54
  vl->id = kVARUV;   vl->name = "varUV";   vl->unit = "m2/s2";  vl->dim = 1;  vl->vec = 1;  vl++;  // 55
  vl->id = kKINE;    vl->name = "kinE";    vl->unit = "m2/s2";  vl->dim = 1;  vl->vec = 1;  vl++;  // 56
  vl->id = kSDEVH;   vl->name = "sdevH";   vl->unit = "m";      vl->dim = 1;  vl->vec = 1;  vl++;  // 57
  vl->id = kVARVT;   vl->name = "varVt";   vl->unit = "m4/s2";  vl->dim = 1;  vl->vec = 1;  vl++;  // 58
  vl->id = kKINER;   vl->name = "kinEr";   vl->unit = "m2/s2";  vl->dim = 1;  vl->vec = 1;  vl++;  // 59
  vl->id = kFLDRATE; vl->name = "fldRate"; vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 60
  vl->id = kKINRATIO;vl->name = "kinRatio";vl->unit = "-";      vl->dim = 1;  vl->vec = 1;  vl++;  // 61
  vl->id = kMAXUV;   vl->name = "maxUV";   vl->unit = "m/s";    vl->dim = 2;  vl->vec = 3;  vl++;  // 62
  vl->id = kMINUV;   vl->name = "minUV";   vl->unit = "m/s";    vl->dim = 2;  vl->vec = 3;  vl++;  // 63
  vl->id = kMAXUS;   vl->name = "maxUs";   vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 64
  vl->id = kMINUS;   vl->name = "minUs";   vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 65
  vl->id = kMAXTAU;  vl->name = "maxTau";  vl->unit = "N/m2";   vl->dim = 1;  vl->vec = 1;  vl++;  // 66
  vl->id = kMAXU;    vl->name = "maxU";    vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 67
  vl->id = kMINU;    vl->name = "minU";    vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 68
  vl->id = kMAXV;    vl->name = "maxV";    vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 69
  vl->id = kMINV;    vl->name = "minV";    vl->unit = "m/s";    vl->dim = 1;  vl->vec = 1;  vl++;  // 70


  vpdata = 0;
  vpcomp = 0;

  vedata = 0;
  vecomp = 0;

  ntmser = 0;
  tmser  = NULL;

  statistics = false;
  statist    = NULL;

  nMacro = 0;

  ngct = 0;
}


PROJECT::~PROJECT()
{
}


void PROJECT::ReplaceMacro( char *name )
{
  ReplaceMacro( name, name );
}


void PROJECT::ReplaceMacro( char *src, char *dst, int tmstep, int pid )
{
  char tmp[500];
  char frm[500];
  strncpy( tmp, src, 500 );

  for( int i=0; i<nMacro; )
  {
    bool found = false;
    int  tm    = strcmp( macro_key[i], "<tm>" );
    int  part  = strcmp( macro_key[i], "<part>" );
    int  hhmm  = strcmp( macro_key[i], "<time>" );

    // handle special macros for timestep-no and partition-id
    if( tmstep > 0  &&  tm == 0 )
    {
      found = ReplaceKey( tmp, frm, macro_key[i], macro_value[i] );
      sprintf( tmp, frm, tmstep );
    }
    else if( pid > 0  &&  part == 0 )
    {
      found = ReplaceKey( tmp, frm, macro_key[i], macro_value[i] );
      sprintf( tmp, frm, pid );
    }
    else if( hhmm == 0 )
    {
      time_t now  = time( NULL );
      struct tm *tminfo = localtime( &now );

      // some selected format specifiers used in strftime()
      //
      // --- date
      // %a  Abbreviated weekday name                                 * Thu
      // %A  Full weekday name                                        * Thursday
      // %b  Abbreviated month name                                   * Aug
      // %B  Full month name                                          * August
      // %c  Date and time representation                             * Thu Aug 23 14:55:02 2001
      // %d  Day of the month, zero-padded (01-31)                      23
      // %F  Short YYYY-MM-DD date, equivalent to %Y-%m-%d              2001-08-23
      // %g  Week-based year, last two digits (00-99)                   01
      // %G  Week-based year                                            2001
      // %j  Day of the year (001-366)                                  235
      // %m  Month as a decimal number (01-12)                          08
      // %u  ISO 8601 weekday as number with Monday as 1 (1-7)          4
      // %V  ISO 8601 week number (00-53)                               34
      // %w  Weekday as a decimal number with Sunday as 0 (0-6)         4
      // %x  Date representation                                      * 08/23/01
      // %y  Year, last two digits (00-99)                              01
      // %Y  Year                                                       2001
      //
      // --- time
      // %H  Hour in 24h format (00-23)                                 14
      // %I  Hour in 12h format (01-12)                                 02
      // %M  Minute (00-59)                                             55
      // %p  AM or PM designation                                       PM
      // %r  12-hour clock time                                       * 02:55:02 pm
      // %R  24-hour HH:MM time, equivalent to %H:%M                    14:55
      // %S  Second (00-61)                                             02
      // %T  ISO 8601 time format (HH:MM:SS), equivalent to %H:%M:%S    14:55:02
      // %X  Time representation                                      * 14:55:02
      //
      // --- others
      // %n  New-line character                                         '\n'
      // %t  Horizontal-tab character                                   '\t'
      // %%  A % sign                                                   %

      char   tmtext[50];
      strftime( tmtext, 50, macro_value[i], tminfo );

      found = ReplaceKey( tmp, tmp, macro_key[i], tmtext );
    }
    else if( tm != 0  &&  part != 0 )
    {
      found = ReplaceKey( tmp, tmp, macro_key[i], macro_value[i] );
    }

    if( !found ) i++;
  }

  strncpy( dst, tmp, 500 );
}


void PROJECT::ReplaceAllKeys( char *src, char *dst, char *key, char *val )
{
  // NOTE: On windows Rismo crashes when calling ReplaceKey from this point.
  // return;
  while( ReplaceKey(src,dst,key,val) );
}


bool PROJECT::ReplaceKey( char *src, char *dst, char *key, char *val )
{
  char tmp[500];

  char *c;
  char *t = tmp;
  bool  found = false;

  while( *src )
  {
    char *k = key;
    while( *src  &&  *src != *key ) *t++ = *src++;

    c = src;
    while( *c  &&  *k )
    {
      if( *c++ != *k++ ) break;
    }

    if( *k ) *t++ = *src++;             // macro key not found
    else      break;
  }

  if( *src )                            // replace key by val
  {
    found = true;

    src = c;
    while( *val ) *t++ = *val++;
  }

  while( *src ) *t++ = *src++;
  *t = '\0';

  strncpy( dst, tmp, 500 );

  return found;
}


void PROJECT::SetPartKey( char *filename )
{
  if( subdom.npr <= 1 )
  {
    char macro_part[10];
    char macro_zero[10];
    strcpy( macro_part, "<part>" );
    strcpy( macro_zero, "" );

    ReplaceAllKeys( filename, filename, macro_part, macro_zero );
    return;
  }

  if( strstr(filename,"<part>") )
  {
    ReplaceMacro( filename, filename, 0, subdom.pid+1 );
  }
  else
  {
    // copy the filename to tmp[]
    char tmp[500];
    strncpy( tmp, filename, 500 );

    char part[10];
    sprintf( part, "%03d_", subdom.pid+1 );

    // search for the last occurence of "/" in tmp
    char *s = tmp + strlen(tmp)-1;
    while( s != tmp )
    {
      if( *s == '/' ) break;
      s--;
    }

    // copy the first part of the filename from tmp
    char *c = filename;
    char *t = tmp;
    while( t <= s ) *c++ = *t++;

    // copy the partition id
    char *p = part;
    while( *p ) *c++ = *p++;

    // copy the rest of the filename from tmp
    while( *t ) *c++ = *t++;
  }
}


void PROJECT::Input()
{
  M2D->region->dryRew.method      = 2;
  M2D->region->dryRew.dryRewFreq  = 10;
  M2D->region->dryRew.dryLimit    = 0.001;
  M2D->region->dryRew.rewetLimit  = 0.005;
  M2D->region->dryRew.rewetPasses = 100;
  M2D->region->dryRew.countDown   = 5;


  // open the input file -----------------------------------------------------------------

  ASCIIFILE* file = new ASCIIFILE( name.inputFile, "r" );

  if( !file || !file->getid() )
    REPORT::rpt.Error( kOpenFileFault, "%s %s (PROJECT::Input #1)",
                       "can not open input file", name.inputFile );


  // check version of input file ---------------------------------------------------------

  char* textLine = file->next();

  int release = 0;
  sscanf( textLine, " $RISMO2D %d", &release );

  if( subdom.pid == 0 )
    REPORT::rpt.Screen( 2, "\n ### release of input file:              %d ###\n\n", release );

  if( release >= 30900 )   Input_30900( file );
  else                     Input_00000( file );
}

//////////////////////////////////////////////////////////////////////////////////////////

void PROJECT::Input_00000( ASCIIFILE* file )
{
  int   i;
  char* textLine;
  char  text[500];

  file->rewind();


  // read file names ---------------------------------------------------------------------

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    sscanf( textLine, "%s", name.timeStepFile );
  }
  else
  {
    REPORT::rpt.Error( kUserFault, "no time step file specified (PROJECT::Input #2)" );
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    sscanf( textLine, "%s", name.materialFile );
  }
  else
  {
    REPORT::rpt.Error( kUserFault, "no roghness table specified (PROJECT::Input #3)" );
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    sscanf( textLine, "%s", name.regionFile );
  }
  else
  {
    REPORT::rpt.Error( kUserFault, "no region file specified (PROJECT::Input #4)" );
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    sscanf( textLine, "%s", name.controlFile );
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    if( textLine[0] == 'b' )  name.ascii_initial = 0;
    else                      name.ascii_initial = 1;

    sscanf( textLine+2, "%s", name.initialFile );
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    sscanf( textLine, "%s", name.inputStatistFile );
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    sscanf( textLine, "%s", name.subdomFile );
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    int reportLevel = 3;

    if( subdom.pid > 0 )
    {
      sscanf( textLine, "%s %d", text, &reportLevel );
      sprintf( name.reportFile, "%03d_%s", subdom.pid+1, text );
    }
    else
    {
      sscanf( textLine, "%s %d", name.reportFile, &reportLevel );
    }

    if( reportLevel < 1 ) reportLevel = 1;

    REPORT::rpt.Open( name.reportFile, reportLevel, subdom.pid );     // open report file
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    if( subdom.npr > 1 )
    {
      sscanf( textLine, "%s", text );
      sprintf( name.geometryFile, "%03d_%s", subdom.pid+1, text );
    }
    else
    {
      sscanf( textLine, "%s", name.geometryFile );
    }
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    if( textLine[0] == 'b' )  name.ascii_restart = 0;
    else                      name.ascii_restart = 1;

    if( subdom.npr > 1 )
    {
      sscanf( textLine+2, "%s", text );
      sprintf( name.restartFile, "%03d_%s", subdom.pid+1, text );
    }
    else
    {
      sscanf( textLine+2, "%s", name.restartFile );
    }
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    if( subdom.npr > 1 )
    {
      sscanf( textLine, "%s", text );
      sprintf( name.outputStatistFile, "%03d_%s", subdom.pid+1, text );
    }
    else
    {
      sscanf( textLine, "%s", name.outputStatistFile );
    }
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    if( subdom.npr > 1 )
    {
      sscanf( textLine, "%s", text );
      sprintf( name.wetAvsFile, "%03d_%s", subdom.pid+1, text );
    }
    else
    {
      sscanf( textLine, "%s", name.wetAvsFile );
    }
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    if( subdom.npr > 1 )
    {
      sscanf( textLine, "%s", text );
      sprintf( name.rgAvsFile, "%03d_%s", subdom.pid+1, text );
    }
    else
    {
      sscanf( textLine, "%s", name.rgAvsFile );
    }
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    if( subdom.npr > 1 )
    {
      sscanf( textLine, "%s", text );
      sprintf( name.ctAvsFile, "%03d_%s", subdom.pid+1, text );
    }
    else
    {
      sscanf( textLine, "%s", name.ctAvsFile );
    }
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    if( subdom.npr > 1 )
    {
      sscanf( textLine, "%s", text );
      sprintf( name.rgAvsStatistFile, "%03d_%s", subdom.pid+1, text );
    }
    else
    {
      sscanf( textLine, "%s", name.rgAvsStatistFile );
    }
  }

  textLine = file->nextLine();
  if( textLine[0] != ' ' )
  {
    sscanf( textLine, "%s", name.optOutputFile );
  }


  // set variables for output in AVS-UCD-Files -------------------------------------------

  vpdata = 11;
  vpcomp =  9;

  vpoutlist[0] = kUV;
  vpoutlist[1] = kS;
  vpoutlist[2] = kK;
  vpoutlist[3] = kD;
  vpoutlist[4] = kC;
  vpoutlist[5] = kH;
  vpoutlist[6] = kUS;
  vpoutlist[7] = kTAU;
  vpoutlist[8] = kVT;

  vedata = 1;
  vecomp = 1;

  veoutlist[0] = kS;

  // -------------------------------------------------------------------------------------

  REPORT::rpt.Copyright( release, false );

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n  %30s   %s\n  %30s   %s\n",
                          " input...... time step file:",   name.timeStepFile,
                          "              material file:",   name.materialFile,
                          "                region file:",   name.regionFile,
                          "               control file:",   name.controlFile );

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n  %30s   %s\n\n",
                          "               initial file:",   name.initialFile,
                          "             statistic file:",   name.inputStatistFile,
                          "             subdomain file:",   name.subdomFile );

  REPORT::rpt.Message( 0, "  %30s   %s / %d\n  %30s   %s\n",
                          "output......... report file:",   name.reportFile, REPORT::rpt.level,
                          "         region output file:",   name.geometryFile );

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n",
                          "               restart file:",   name.restartFile,
                          "             statistic file:",   name.outputStatistFile );

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n  %30s   %s\n",
                          "             dry/rewet file:",   name.wetAvsFile,
                          "                region file:",   name.rgAvsFile,
                          "               control file:",   name.ctAvsFile);

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n\n",
                          "      statistic region file:",   name.rgAvsStatistFile,
                          "       optional output file:",   name.optOutputFile );

  REPORT::rpt.Line1( 0 );


  if( name.inputStatistFile[0]  ||  name.outputStatistFile[0]
                                ||  name.rgAvsStatistFile[0]  )
  {
    statistics = true;
  }


  // -------------------------------------------------------------------------------------
  // read number of solver specifications

  textLine = file->nextLine();
  sscanf( textLine, "%d", &SOLVER::m_neqs );

  SOLVER::m_solver = new SOLVER* [SOLVER::m_neqs];
  if( !SOLVER::m_solver )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (PROJECT::Input #5)" );


  // -------------------------------------------------------------------------------------
  // read switches for equation solver
  //                --- direct solvers
  //     types       1: kFront            frontal solver
  //
  //                --- conjugate gradient solvers
  //                 5: kBicgstab
  //                 6: kParmsBcgstabd    PARMS version of BiCGStab
  //                 7: kParmsFgmresd     PARMS version of Gmres

  sprintf( text, "   reading definition for %d solvers\n", SOLVER::m_neqs );
  REPORT::rpt.Output( text, 3 );

  for( i=0; i<SOLVER::m_neqs; i++ )
  {
    int no, type;

    textLine = file->nextLine();
    sscanf( textLine, "%d %d", &no, &type );

    SOLVER::m_solver[i] = NULL;

    switch( type )
    {
      default:
        REPORT::rpt.Error( kParameterFault, "solver type not supported (PROJECT::Input #6)" );

      case kFront:           SOLVER::m_solver[i] = new FRONT();      break;
      case kFrontm:          SOLVER::m_solver[i] = new FRONTM();     break;
      case kBicgstab:        SOLVER::m_solver[i] = new BICGSTAB();   break;
      case kParmsBcgstabd:   SOLVER::m_solver[i] = new P_BCGSTABD(); break;
      case kParmsFgmresd:    SOLVER::m_solver[i] = new P_FGMRESD();  break;
    }

    if( !SOLVER::m_solver[i] )
      REPORT::rpt.Error( kMemoryFault, "can not allocate memory (PROJECT::Input #7)" );

    SOLVER::m_solver[i]->no = no;

    switch( type )
    {
      case kFront:
        sscanf( textLine, "%d %d %d %d %s",
                &no, &type, &SOLVER::m_solver[i]->mfw,
                            &SOLVER::m_solver[i]->size,
                             SOLVER::m_solver[i]->path );

        sprintf( text, "\n %d. %s\n",
                 i+1, "solver specification: frontal solver" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %s\n",
                 "mfw:",   SOLVER::m_solver[i]->mfw,
                 "size:",  SOLVER::m_solver[i]->size,
                 "path:",  SOLVER::m_solver[i]->path );
        REPORT::rpt.Output( text, 4 );
        break;

      case kFrontm:
        sscanf( textLine, "%d %d %d %d %d %s",
                &no, &type, &SOLVER::m_solver[i]->mceq,
                            &SOLVER::m_solver[i]->mfw,
                            &SOLVER::m_solver[i]->size,
                             SOLVER::m_solver[i]->path );

        sprintf( text, "\n %d. %s\n",
                 i+1, "solver specification: frontal matrix solver" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %s\n",
                 "mceq:",  SOLVER::m_solver[i]->mceq,
                 "mfw:",   SOLVER::m_solver[i]->mfw,
                 "size:",  SOLVER::m_solver[i]->size,
                 "path:",  SOLVER::m_solver[i]->path );
        REPORT::rpt.Output( text, 4 );
        break;

      case kBicgstab:
        sscanf( textLine, "%d %d %d %d %d %d %lf",
                &no, &type, &SOLVER::m_solver[i]->preconType,
                            &SOLVER::m_solver[i]->proceed,
                            &SOLVER::m_solver[i]->mceq,
                            &SOLVER::m_solver[i]->maxIter,
                            &SOLVER::m_solver[i]->maxDiff );

        sprintf( text, "\n %d. %s\n",
                 i+1, "solver specification: BiCGStab" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %lf\n",
                 "precon type:", SOLVER::m_solver[i]->preconType,
                 "proceed:",     SOLVER::m_solver[i]->proceed,
                 "mceq:",        SOLVER::m_solver[i]->mceq,
                 "maxIter:",     SOLVER::m_solver[i]->maxIter,
                 "maxDiff:",     SOLVER::m_solver[i]->maxDiff );
        REPORT::rpt.Output( text, 4 );
        break;

      case kParmsBcgstabd:
        sscanf( textLine, "%d %d %d %d %d %d %lf",
                &no, &type, &SOLVER::m_solver[i]->preconType,
                            &SOLVER::m_solver[i]->proceed,
                            &SOLVER::m_solver[i]->mceq,
                            &SOLVER::m_solver[i]->maxIter,
                            &SOLVER::m_solver[i]->maxDiff );

        sprintf( text, "\n %d. %s\n",
                 i+1, "solver specification: PARMS - BiCGStab" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %lf\n",
                 "precon type:", SOLVER::m_solver[i]->preconType,
                 "proceed:",     SOLVER::m_solver[i]->proceed,
                 "mceq:",        SOLVER::m_solver[i]->mceq,
                 "maxIter:",     SOLVER::m_solver[i]->maxIter,
                 "maxDiff:",     SOLVER::m_solver[i]->maxDiff );
        REPORT::rpt.Output( text, 4 );
        break;

      case kParmsFgmresd:
        sscanf( textLine, "%d %d %d %d %d %d %d %lf",
                &no, &type, &SOLVER::m_solver[i]->preconType,
                            &SOLVER::m_solver[i]->proceed,
                            &SOLVER::m_solver[i]->mceq,
                            &SOLVER::m_solver[i]->mkyrl,
                            &SOLVER::m_solver[i]->maxIter,
                            &SOLVER::m_solver[i]->maxDiff );

        sprintf( text, "\n %d. %s\n",
                 i+1, "solver specification: PARMS - flexible Gmres" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %lf\n",
                 "precon type:", SOLVER::m_solver[i]->preconType,
                 "proceed:",     SOLVER::m_solver[i]->proceed,
                 "mceq:",        SOLVER::m_solver[i]->mceq,
                 "mkyrl:",       SOLVER::m_solver[i]->mkyrl,
                 "maxIter:",     SOLVER::m_solver[i]->maxIter,
                 "maxDiff:",     SOLVER::m_solver[i]->maxDiff );
        REPORT::rpt.Output( text, 4 );
        break;
    }
  }

  REPORT::rpt.OutputLine1( 3 );


  // read element type specifications ----------------------------------------------------

  sprintf( text, "   reading material file: %s\n", name.materialFile );
  REPORT::rpt.Output( text, 3 );

  TYPE::Import( name.materialFile );
  REPORT::rpt.OutputLine1( 3 );


  // -------------------------------------------------------------------------------------
  // read section coordinates for 1) reordering of elements and
  //                              2) initialization of water surface elevation

  textLine = file->nextLine();
  sscanf( textLine, " %d",  &nSection );

  section = new SECTION [nSection];
  if( !section )
    REPORT::rpt.Error( "can not allocate memory (PROJECT::input #8)" );


  sprintf( text, "   reading %d section coordinates\n", nSection );
  REPORT::rpt.Output( text, 3 );

  for( i=0; i<nSection; i++ )
  {
    double xs, ys, xe, ye, z;

    textLine = file->nextLine();
    sscanf( textLine, " %lf %lf %lf %lf %lf", &xs, &ys, &xe, &ye, &z );

    section[i].init( xs, ys, xe, ye, z );

    sprintf( text, "  %3d:  %12.3lf  %12.3lf  %12.3lf  %12.3lf    %8.3lf\n",
                   i+1, section[i].getxs(), section[i].getys(),
                        section[i].getxe(), section[i].getye(),
                        section[i].z );
    REPORT::rpt.Output( text, 5 );
  }

  REPORT::rpt.OutputLine1( 3 );


  // read convergence limits -------------------------------------------------------------

  textLine = file->nextLine();
  sscanf( textLine, " %lf %lf %lf", &(convUV), &(convS), &(convKD) );

  sprintf( text, "  %30s  %9.2le\n  %30s  %9.2le\n  %30s  %9.2le\n",
                 "convergence limits - UV:",  convUV,
                 "                      S:",  convS,
                 "                     KD:",  convKD );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // read min and max value for node parameters ------------------------------------------

  dep_minVt = 0.0;            // deprecated
  double maxVt, maxK, maxD;   // deprecated

  minK = minD = 0.0;
  maxUs = 100.0;
  minC  = 0.0;

  textLine = file->nextLine();
  sscanf( textLine, " %lf %lf %lf %lf %lf %lf %lf %lf",
                    &dep_minVt, &maxVt, &minK, &maxK, &minD, &maxD, &maxUs, &minC );

  sprintf( text, "  %30s  %le\n  %30s  %le\n  %30s  %le\n",
                 "vt - minimum value (depr.):",  dep_minVt,
                 " K - minimum value        :",  minK,
                 " D - minimum value        :",  minD );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %le\n  %30s  %le\n\n",
                 "maximum flow velocity:",  maxUs,
                 "minimum concentration:",  minC );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // read relaxation parameter -----------------------------------------------------------

  relaxMethod = 1;
  relaxMin    = 0.001;
  relaxMax    = 0.1;
  maxDeltaUV  = 0.01;
  maxDeltaS   = 0.1;
  maxDeltaKD  = 0.01;

  textLine = file->nextLine();
  sscanf( textLine, " %d %lf %lf %lf %lf %lf", &(relaxMethod),
                                               &(relaxMin),
                                               &(relaxMax),
                                               &(maxDeltaUV),
                                               &(maxDeltaS),
                                               &(maxDeltaKD) );

  sprintf( text, "  %30s  %d\n  %30s  %9.6lf\n  %30s  %9.6lf\n\n",
                 "relaxation...  method:",  relaxMethod,
                 "              minimum:",  relaxMin,
                 "              maximum:",  relaxMax );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.6lf\n  %30s  %9.6lf\n  %30s  %9.6lf\n\n",
                 "     max change of UV:",  maxDeltaUV,
                 "     max change of S :",  maxDeltaS,
                 "     max change of KD:",  maxDeltaKD );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // read parameters for dry/rewet algorithm ---------------------------------------------

  DRYREW* dryRew = &M2D->region->dryRew;

  textLine = file->nextLine();
  sscanf( textLine, " %d %d %lf %lf %d %d", &(dryRew->method),
          &(dryRew->dryRewFreq),
          &(dryRew->dryLimit),
          &(dryRew->rewetLimit),
          &(dryRew->rewetPasses),
          &(dryRew->countDown) );

  if( dryRew->rewetLimit < dryRew->dryLimit ) dryRew->rewetLimit = dryRew->dryLimit;

  sprintf( text, "  %30s  %4d\n  %30s  %4d\n  %30s  %4d\n  %30s  %4d\n\n",
                 "dry/rewet method:",         dryRew->method,
                 "frequency of dry/rewet:",   dryRew->dryRewFreq,
                 "number of rewet passes:",   dryRew->rewetPasses,
                 "count down for rewetting:", dryRew->countDown );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.6lf\n  %30s  %9.6lf\n\n",
                 "limit for drying:",    dryRew->dryLimit,
                 "limit for rewetting:", dryRew->rewetLimit );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );

  hmin = dryRew->dryLimit / 2.0;


  // read number of smoothing passes -----------------------------------------------------

  smoothPassesBC = 0;
  smoothPassesKD = 0;
  smoothPassesVT = 0;

  textLine = file->nextLine();
  sscanf( textLine, " %d %d %d", &(smoothPassesBC), &(smoothPassesKD), &(smoothPassesVT) );

  sprintf( text, "  %30s  %-d\n  %30s  %-d\n  %30s  %-d\n\n",
                 "BC smoothing passes:",  smoothPassesBC,
                 "KD smoothing passes:",  smoothPassesKD,
                 "vt smoothing passes:",  smoothPassesVT );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // read type of boundary model for different boundaries --------------------------------

  textLine = file->nextLine();
  sscanf( textLine, " %d %d %d", &(KDBcon[kIndInflow]),
                                 &(KDBcon[kIndOutflow]),
                                 &(KDBcon[kIndSide]) );

  sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n\n",
                 "KD boundary model for inflow:", KDBcon[kIndInflow],
                 "                     outflow:", KDBcon[kIndOutflow],
                 "                        side:", KDBcon[kIndSide] );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // read degree of GAUSS point integration ----------------------------------------------

  textLine = file->nextLine();
  sscanf( textLine, " %d", &(GPdeg) );

  sprintf( text, "  %30s  %4d\n\n",
           "degree of integration:", GPdeg );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  SHAPE::init( GPdeg );


  // read parameters for bookkeeping algorithm -------------------------------------------

  textLine = file->nextLine();
  sscanf( textLine, "%lf", &sed.maxDz );
  sprintf( text, "   %30s %6.3lf\n\n",
                 "max. change in bottom elev.: ", sed.maxDz );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // -------------------------------------------------------------------------------------
  // read constants
  // kinematic viscosity, density, gravity, von Karman's constant and  wall distance
  //
  // the kinematic viscosity is a function of material, temperature (T �C)
  // and pressure; in case of water (10 degree Celsius) it is approximately
  //
  //             vk = 1.78e-6 / (1.0 + 0.0337*T + 2.2e-4*T*T)

  vk       = 1.78e-6 / (1.0 + 0.337 + 0.022);
  rho      = 1000.0;
  g        = 9.81;
  kappa    = 0.41;
  earthVel = 0.0001458;

  textLine = file->nextLine();
  sscanf(textLine, " %lf %lf %lf %lf %lf", &(vk), &(rho), &(g),&(kappa), &(earthVel) );

  sprintf( text, "  %30s  %9.2le\n  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "kinematic viscosity:",   vk,
          "density:",               rho,
          "gravity:",               g );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.4lf\n",
          "von Karman's constant:", kappa );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.7lf\n\n",
          "earth rotation velocity:", earthVel );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // read k-epsilon constants ------------------------------------------------------------

  textLine = file->nextLine();
  sscanf(textLine, " %lf %lf %lf %lf %lf %lf", &(KD.cm), &(KD.cd),
                   &(KD.sK), &(KD.sD), &(KD.c1D), &(KD.c2D) );

  sprintf(text, "  %30s  %9.4lf\n  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "k-epsilon constants:  cm:",  KD.cm,
          "                      cd:",  KD.cd,
          "                      sK:",  KD.sK );
  REPORT::rpt.Output( text, 3 );

  sprintf(text, "  %30s  %9.4lf\n  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "                      sD:",  KD.sD,
          "                     c1D:",  KD.c1D,
          "                     c2D:",  KD.c2D );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // read scaling factors ------------------------------------------------------------------

  int    typeOfScale = 0;
  double lScale = 1.0;
  double hScale = 1.0;
  double kScale = 1.0;
  double vScale = 1.0;

  textLine = file->nextLine();
  sscanf( textLine, " %d", &(typeOfScale) );

  switch( typeOfScale )
  {
    case 1:
      sscanf( textLine, " %d %lf %lf", &(typeOfScale),
                                       &(lScale),
                                       &(vScale) );
      scale.init( typeOfScale, lScale, vScale );
      break;

    case 2:
      sscanf( textLine, " %d %lf %lf %lf", &(typeOfScale),
                                           &(lScale),
                                           &(hScale),
                                           &(kScale) );
      scale.init( typeOfScale, lScale, hScale, kScale );
      break;
  }

  sprintf( text, "  %30s  %9d\n  %30s  %9.4lf\n",
                 "scaling         :",  scale.gettype(),
                 "          lScale:",  scale.getlScale() );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "                 hScale:",  scale.gethScale(),
          "                ksScale:",  scale.getkScale() );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "                 tScale:",  scale.gettScale(),
          "                 vScale:",  scale.getvScale() );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );

  delete file;
}


//////////////////////////////////////////////////////////////////////////////////////////

void PROJECT::Input_30900( ASCIIFILE* file )
{
  char* textLine;
  char  text[500];
  char  cdummy[100];

  int   nsolver = 0;
  int   reportLevel = 3;

  int   _minVt   = false;
  int   _minVtxx = false;
  int   _minVtyy = false;
  int   _minVtKD = false;
  int   _minVtC  = false;

  outputPath[0] = '\0';

  while( !feof(file->getid()) )
  {
    if( !(textLine = file->nextLine()) )  break;

#   ifdef kDebug
    if( subdom.pid == 0 )  printf( textLine );
#   endif

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

    char macro_tm[10];
    char macro_part[10];
    char macro_zero[10];

    strcpy( macro_tm, "<tm>" );
    strcpy( macro_part, "<part>" );
    strcpy( macro_zero, "" );

    switch( ikey )
    {
      // ---------------------------------------------------------------------------------
      // Macros
      case kCOMMENT:
        {
          char *com = strstr( textLine, "$COMMENT" );
          REPORT::rpt.Output( &com[8], 0 );
        }
        break;
        
      case kMACRO:
        if( nMacro < kMaxMacro )
        {
          sscanf( textLine, "$MACRO %s %s", macro_key[nMacro], macro_value[nMacro] );
          nMacro++;
        }
        break;

      // ---------------------------------------------------------------------------------
      case kREPORTLEVEL:
        {
          sscanf( textLine, "%s %d", cdummy, &reportLevel );
          REPORT::rpt.Setlevel( reportLevel );
          if( reportLevel < 1 ) reportLevel = 1;
        }
        break;

      case kREPORTTIME:
        {
          int  first, last, step;
          int n = sscanf( textLine, "$REPORTTIME %d to %d step %d", &first, &last, &step );
          if( n == 3 )
          {
            REPORT::rpt.first = first;
            REPORT::rpt.last  = last;
            REPORT::rpt.step  = step;

            REPORT::rpt.ntm    = last + 1;
            REPORT::rpt.tmlist = new int[REPORT::rpt.ntm];
            if( !REPORT::rpt.tmlist )
              REPORT::rpt.Error( kMemoryFault, "%s (PROJECT::Input #9)",
                                               "can not allocate memory" );

            for( int i=0;     i<=last; i++    )  REPORT::rpt.tmlist[i] = false;
            for( int i=first; i<=last; i+=step ) REPORT::rpt.tmlist[i] = true;
          }
        }
        break;

      // ---------------------------------------------------------------------------------
      // Filenames for INPUT-Data
      case kTIMESTEPFILE:
        sscanf( textLine, "%s %s", cdummy, name.timeStepFile );
        ReplaceMacro( name.timeStepFile );
        ReplaceAllKeys( name.timeStepFile, name.timeStepFile, macro_tm, macro_zero );
        ReplaceAllKeys( name.timeStepFile, name.timeStepFile, macro_part, macro_zero );
        break;

      case kMATERIALFILE:
        sscanf( textLine, "%s %s", cdummy, name.materialFile );
        ReplaceMacro( name.materialFile );
        ReplaceAllKeys( name.materialFile, name.materialFile, macro_tm, macro_zero );
        ReplaceAllKeys( name.materialFile, name.materialFile, macro_part, macro_zero );
        break;

      case kREGIONFILE:
        sscanf( textLine, "%s %s", cdummy, name.regionFile );
        ReplaceMacro( name.regionFile );
        ReplaceAllKeys( name.regionFile, name.regionFile, macro_tm, macro_zero );
        ReplaceAllKeys( name.regionFile, name.regionFile, macro_part, macro_zero );
        break;

      case kCONTROLFILE:
        sscanf( textLine, "%s %s", cdummy, name.controlFile );
        ReplaceMacro( name.controlFile );
        ReplaceAllKeys( name.controlFile, name.controlFile, macro_tm, macro_zero );
        ReplaceAllKeys( name.controlFile, name.controlFile, macro_part, macro_zero );
        break;

      case kASC_INITFILE:
        sscanf( textLine, "%s %s", cdummy, name.initialFile );
        name.ascii_initial = true;
        ReplaceMacro( name.initialFile );
        ReplaceAllKeys( name.initialFile, name.initialFile, macro_tm, macro_zero );
        ReplaceAllKeys( name.initialFile, name.initialFile, macro_part, macro_zero );
        break;

      case kBIN_INITFILE:
        sscanf( textLine, "%s %s", cdummy, name.initialFile );
        name.ascii_initial = false;
        ReplaceMacro( name.initialFile );
        ReplaceAllKeys( name.initialFile, name.initialFile, macro_tm, macro_zero );
        ReplaceAllKeys( name.initialFile, name.initialFile, macro_part, macro_zero );
        break;

      case kSTA_INITFILE:
        sscanf( textLine, "%s %s", cdummy, name.inputStatistFile );
        ReplaceMacro( name.inputStatistFile );
        ReplaceAllKeys( name.inputStatistFile, name.inputStatistFile, macro_tm, macro_zero );
        ReplaceAllKeys( name.inputStatistFile, name.inputStatistFile, macro_part, macro_zero );
        statistics = true;
        break;

      case kSUBDOMFILE:
        sscanf( textLine, "%s %s", cdummy, name.subdomFile );
        ReplaceMacro( name.subdomFile );
        ReplaceAllKeys( name.subdomFile, name.subdomFile, macro_tm, macro_zero );
        ReplaceAllKeys( name.subdomFile, name.subdomFile, macro_part, macro_zero );
        break;

      // ---------------------------------------------------------------------------------
      // Filenames for OUTPUT-Data
      case kOUTPUTPATH:
        sscanf( textLine, "%s %s", cdummy, outputPath );
        ReplaceMacro( outputPath );
        break;

      case kREPORTFILE:
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.reportFile, "%s%s", outputPath, text );

        ReplaceAllKeys( name.reportFile, name.reportFile, macro_tm, macro_zero );
        SetPartKey( name.reportFile );

        REPORT::rpt.Open( name.reportFile, reportLevel, subdom.pid );
        REPORT::rpt.Copyright( release, false );
        break;

      case kASC_RESTFILE:
        name.ascii_restart = true;
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.restartFile, "%s%s", outputPath, text );
        SetPartKey( name.restartFile );
        break;

      case kBIN_RESTFILE:
        name.ascii_restart = false;
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.restartFile, "%s%s", outputPath, text );
        SetPartKey( name.restartFile );
        break;

      case kSTA_RESTFILE:
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.outputStatistFile, "%s%s", outputPath, text );
        SetPartKey( name.outputStatistFile );
        statistics = true;
        break;

      case kRG_UCDFILE:
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.rgAvsFile, "%s%s", outputPath, text );
        SetPartKey( name.rgAvsFile );
        break;

      case kCT_UCDFILE:
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.ctAvsFile, "%s%s", outputPath, text );
        SetPartKey( name.ctAvsFile );
        break;

      case kCN_UCDFILE:
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.redAvsFile, "%s%s", outputPath, text );
        SetPartKey( name.redAvsFile );
        break;

      case kWN_UCDFILE:
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.wetAvsFile, "%s%s", outputPath, text );
        SetPartKey( name.wetAvsFile );
        break;

      case kST_UCDFILE:
        sscanf( textLine, "%s %s", cdummy, text );
        ReplaceMacro( text );
        sprintf( name.rgAvsStatistFile, "%s%s", outputPath, text );
        SetPartKey( name.rgAvsStatistFile );
        statistics = true;
        break;

      case kTIMESERIES:
        {
          TMSER *newtmser = new TMSER[ntmser+1];
          if( !newtmser )
            REPORT::rpt.Error( kMemoryFault, "%s (PROJECT::Input #10)",
                                             "can not allocate memory" );
          for( int i=0; i<ntmser; i++ )
          {
            newtmser[i] = tmser[i];
          }

          if( ntmser ) delete[] tmser;

          ntmser++;
          tmser = newtmser;

          int  first, last, step;
          char filename[500], vars[200];
          int  n = sscanf( textLine, "$TIMESERIES %s %s %d to %d step %d",
                                     filename, vars, &first, &last, &step );
          ReplaceMacro( filename );
          sprintf( tmser[ntmser-1].filename, "%s%s", outputPath, filename );
          SetPartKey( tmser[ntmser-1].filename );

          char  list[100];
          char  seps[] = " ,\t\n\r";
          char* token;

          strcpy( list, vars );
          token = strtok( list, seps );

          while( token != NULL )
          {
            for( int i=0; i<nval; i++ )
            {
              if( strcmp(token, valist[i].name) == 0 )
              {
                tmser[ntmser-1].vcomp++;
                tmser[ntmser-1].vdata += valist[i].dim;
                break;
              }
            }

            token = strtok( NULL, seps );
          }

          tmser[ntmser-1].voutlist = new int[tmser[ntmser-1].vcomp];
          if( !tmser[ntmser-1].voutlist )
            REPORT::rpt.Error( kMemoryFault, "%s (PROJECT::Input #11)",
                                             "can not allocate memory" );

          strcpy( list, vars );
          token = strtok( list, seps );

          int j = 0;
          while( token != NULL )
          {
            for( int i=0; i<nval; i++ )
            {
              if( strcmp(token, valist[i].name) == 0 )
              {
                tmser[ntmser-1].voutlist[j] = i;
                break;
              }
            }

            token = strtok( NULL, seps );
            j++;
          }

          if( n == 5 )
          {
            tmser[ntmser-1].first = first;
            tmser[ntmser-1].last  = last;
            tmser[ntmser-1].step  = step;

            tmser[ntmser-1].ntm    = last + 1;
            tmser[ntmser-1].tmlist = new int[tmser[ntmser-1].ntm];
            if( !tmser[ntmser-1].tmlist )
              REPORT::rpt.Error( kMemoryFault, "%s (PROJECT::Input #12)",
                                               "can not allocate memory" );

            for( int i=0;     i<=last; i++    )  tmser[ntmser-1].tmlist[i] = false;
            for( int i=first; i<=last; i+=step ) tmser[ntmser-1].tmlist[i] = true;
          }
        }
        break;

      // ---------------------------------------------------------------------------------
      case kOUTPUTVARS:  // read variables for output in AVS-UCD-Files
        {
          char  list[200];
          char  seps[] = " ,\t\n\r";
          char* token;

          strcpy( list, textLine+11 );
          token = strtok( list, seps );

          while( token != NULL )
          {
            for( int i=0; i<nval; i++ )
            {
              if( strcmp(token, valist[i].name) == 0 )
              {
                vpoutlist[vpcomp] = i;
                vpcomp++;
                vpdata += valist[i].vec;

                if(    valist[i].id == kS  || valist[i].id == kDZ
                    || valist[i].id == kRC || valist[i].id == kD90 || valist[i].id == kD50
                    || valist[i].id == kDP || valist[i].id == kSP )
                {
                  veoutlist[vecomp] = i;
                  vecomp++;
                  vedata += valist[i].vec;
                }

                if(    valist[i].id == kMEANUV || valist[i].id == kMEANH
                    || valist[i].id == kMEANVT || valist[i].id == kVARU
                    || valist[i].id == kVARV   || valist[i].id == kVARUV
                    || valist[i].id == kKINE   || valist[i].id == kSDEVH
                    || valist[i].id == kVARVT  || valist[i].id == kMEANUS
                    || valist[i].id == kMEANS  || valist[i].id == kMAXTAU
                    || valist[i].id == kMAXUV  || valist[i].id == kMINUV
                    || valist[i].id == kMAXUS  || valist[i].id == kMINUS )
                {
                  statistics = true;
                }
                break;
              }
            }

            token = strtok( NULL, seps );
          }
        }
        break;

      // ---------------------------------------------------------------------------------
      case kSOLVER:  // read solver specifications
        {
          if( SOLVER::m_neqs >= nsolver )
          {
            SOLVER** old_solver = NULL;
            if( nsolver > 0 )  old_solver = SOLVER::m_solver;

            SOLVER::m_solver = new SOLVER* [nsolver+10];
            if( !SOLVER::m_solver )
              REPORT::rpt.Error( kMemoryFault, "%s (PROJECT::Input #13)",
                                               "can not allocate memory" );

            for( int i=0; i<nsolver; i++ )  SOLVER::m_solver[i] = old_solver[i];
            nsolver += 10;
          }


          // -----------------------------------------------------------------------------
          // read switches for equation solver
          //                --- direct solvers
          //     types       1: kFront            frontal solver
          //
          //                --- conjugate gradient solvers
          //                 5: kBicgstab
          //                 6: kParmsBcgstabd    PARMS version of BiCGStab
          //                 7: kParmsFgmresd     PARMS version of Gmres

          int no, type;

          sscanf( textLine, "$SOLVER %d %d", &no, &type );

          SOLVER::m_solver[SOLVER::m_neqs] = NULL;

          switch( type )
          {
            default:
              REPORT::rpt.Error( kParameterFault, "%s (PROJECT::Input #15)",
                                 "solver type not supported" );

            case kFront:           SOLVER::m_solver[SOLVER::m_neqs] = new FRONT();      break;
            case kFrontm:          SOLVER::m_solver[SOLVER::m_neqs] = new FRONTM();     break;
            case kBicgstab:        SOLVER::m_solver[SOLVER::m_neqs] = new BICGSTAB();   break;
            case kParmsBcgstabd:   SOLVER::m_solver[SOLVER::m_neqs] = new P_BCGSTABD(); break;
            case kParmsFgmresd:    SOLVER::m_solver[SOLVER::m_neqs] = new P_FGMRESD();  break;
          }

          SOLVER::m_solver[SOLVER::m_neqs]->no = no;

          if( !SOLVER::m_solver[SOLVER::m_neqs] )
            REPORT::rpt.Error( kMemoryFault, "%s (PROJECT::Input #16)",
                               "can not allocate memory" );

          switch( type )
          {
            case kFront:
              sscanf( textLine, "$SOLVER %d %d %d %d %s",
                      &no, &type, &SOLVER::m_solver[SOLVER::m_neqs]->mfw,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->size,
                                   SOLVER::m_solver[SOLVER::m_neqs]->path );
              break;

            case kFrontm:
              sscanf( textLine, "$SOLVER %d %d %d %d %d %s",
                      &no, &type, &SOLVER::m_solver[SOLVER::m_neqs]->mceq,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->mfw,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->size,
                                  SOLVER::m_solver[SOLVER::m_neqs]->path );
              break;

            case kBicgstab:
              sscanf( textLine, "$SOLVER %d %d %d %d %d %d %lf",
                      &no, &type, &SOLVER::m_solver[SOLVER::m_neqs]->preconType,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->proceed,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->mceq,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->maxIter,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->maxDiff );
              break;

            case kParmsBcgstabd:
              sscanf( textLine, "$SOLVER %d %d %d %d %d %d %lf",
                      &no, &type, &SOLVER::m_solver[SOLVER::m_neqs]->preconType,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->proceed,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->mceq,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->maxIter,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->maxDiff );
              break;

            case kParmsFgmresd:
              sscanf( textLine, "$SOLVER %d %d %d %d %d %d %d %lf",
                      &no, &type, &SOLVER::m_solver[SOLVER::m_neqs]->preconType,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->proceed,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->mceq,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->mkyrl,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->maxIter,
                                  &SOLVER::m_solver[SOLVER::m_neqs]->maxDiff );
              break;
          }

          SOLVER::m_neqs++;
        }
        break;

      // ---------------------------------------------------------------------------------
      case kCONVERGENCE:
        sscanf( textLine, "$CONVERGENCE %lf %lf %lf %lf",
                          &convUV, &convS, &convKD, &sed.convQb );
        break;

      // ---------------------------------------------------------------------------------
      case kMINMAX:
        sscanf( textLine, "$MINMAX %lf %lf %lf %lf %lf",
                          &dep_minVt, &minK, &minD, &maxUs, &minC );
        _minVt = _minVtKD = _minVtC = true;
        break;

      case kMIN_VT:
        sscanf( textLine, "$MIN_VT %lf", &dep_minVt );
        _minVt = true;
        break;

      case kMIN_VT_XX:
        sscanf( textLine, "$MIN_VT_XX %lf", &dep_minVtxx );
        _minVtxx = true;
        break;

      case kMIN_VT_YY:
        sscanf( textLine, "$MIN_VT_YY %lf", &dep_minVtyy );
        _minVtyy = true;
        break;

      case kMIN_VT_KD:
        sscanf( textLine, "$MIN_VT_KD %lf", &minVtKD );
        _minVtKD = true;
        break;

      case kMIN_VT_C:
        sscanf( textLine, "$MIN_VT_C %lf", &minVtC );
        _minVtC = true;
        break;

      case kMIN_K:
        sscanf( textLine, "$MIN_K %lf", &minK );
        break;

      case kMIN_D:
        sscanf( textLine, "$MIN_D %lf", &minD );
        break;

      case kMIN_C:
        sscanf( textLine, "$MIN_C %lf", &minC );
        break;

      case kMAX_Us:
        sscanf( textLine, "$MAX_Us %lf", &maxUs );
        break;

      // ---------------------------------------------------------------------------------
      case kRELAX:
        sscanf( textLine, "$RELAX %d %lf %lf %lf %lf %lf", &relaxMethod,
                                                           &relaxMin,
                                                           &relaxMax,
                                                           &maxDeltaUV,
                                                           &maxDeltaS,
                                                           &maxDeltaKD );
        break;

      // ---------------------------------------------------------------------------------
      case kDRYREW:
        {
          DRYREW* dryRew = &M2D->region->dryRew;
          sscanf( textLine, "$DRYREW %d %d %lf %lf %d %d", &dryRew->method,
                                                           &dryRew->dryRewFreq,
                                                           &dryRew->dryLimit,
                                                           &dryRew->rewetLimit,
                                                           &dryRew->rewetPasses,
                                                           &dryRew->countDown );
          if( dryRew->rewetLimit < dryRew->dryLimit ) dryRew->rewetLimit = dryRew->dryLimit;
        }
        break;

      // ---------------------------------------------------------------------------------
      case kSMOOTH:
          sscanf( textLine, "$SMOOTH %d %d %d", &(smoothPassesBC),
                                                &(smoothPassesKD),
                                                &(smoothPassesVT) );
        break;

      // ---------------------------------------------------------------------------------
      case kKDBOUNDARY:
        sscanf( textLine, "$KDBOUNDARY %d %d %d", &(KDBcon[kIndInflow]),
                                                  &(KDBcon[kIndOutflow]),
                                                  &(KDBcon[kIndSide]) );
        break;

      // ---------------------------------------------------------------------------------
      case kGPDEGREE:
        sscanf( textLine, "$GPDEGREE %d", &(GPdeg) );
        break;

      // ---------------------------------------------------------------------------------
      case kTEMPERATURE:
        sscanf( textLine, "$TEMPERATURE %lf", &celsius );
        vk = 1.78e-6 / (1.0 + 0.0337*celsius + 0.00022*celsius*celsius);
        break;

      case kVISCOSITY:
        sscanf( textLine, "$VISCOSITY %lf", &vk );
        break;

      case kDENSITY:
        sscanf( textLine, "$DENSITY %lf", &rho );
        break;

      case kGRAVITY:
        sscanf( textLine, "$GRAVITY %lf", &g );
        break;

      case kVON_KARMAN:
        sscanf( textLine, "$VON_KARMAN %lf", &kappa );
        break;

      case kEARTH_ROTATION:
        sscanf( textLine, "$EARTH_ROTATION %lf", &earthVel );
        break;

      case kLATITUDE:
        sscanf( textLine, "$LATITUDE %lf", &latitude );
        break;

      // ---------------------------------------------------------------------------------
      case kKDCONST:
        sscanf( textLine, "$KDCONST %lf %lf %lf %lf %lf %lf", &(KD.cm),  &(KD.cd),
                                                              &(KD.sK),  &(KD.sD),
                                                              &(KD.c1D), &(KD.c2D) );
        break;

      // ---------------------------------------------------------------------------------
      case kMUE_SF:
        sscanf( textLine, "$MUE_SF %lf", &mueSf );
        break;

      case kMINU_SF:
        sscanf( textLine, "$MINU_SF %lf", &minUSf );
        break;

      case kMAXTAN_SF:
        sscanf( textLine, "$MAXTAN_SF %lf", &maxTanSf );
        break;

      // ---------------------------------------------------------------------------------
      case kSED_RHOB:
        sscanf( textLine, "$SED_RHOB %d %lf", &sed.nfrac, &sed.rhob );
        break;

      case kSED_M:
        sscanf( textLine, "$SED_M %d %lf", &sed.nfrac, &sed.M );
        break;

      case kSED_TAUC:
        sscanf( textLine, "$SED_TAUC %d %lf", &sed.nfrac, &sed.tauc );
        break;

      case kSED_TAUS:
        sscanf( textLine, "$SED_TAUS %d %lf", &sed.nfrac, &sed.taus );
        break;

      case kSED_US:
        sscanf( textLine, "$SED_US %d %lf", &sed.nfrac, &sed.us );
        break;

      case kSED_D50:
        sscanf( textLine, "$SED_D50 %d %lf", &sed.nfrac, &sed.d50 );
        break;

      case kSED_D90:
        sscanf( textLine, "$SED_D90 %d %lf", &sed.nfrac, &sed.d90 );
        break;

      case kSED_POR:
        sscanf( textLine, "$SED_POR %d %lf", &sed.nfrac, &sed.por );
        break;

      case kSED_PHIR:
        sscanf( textLine, "$SED_PHIR %d %lf", &sed.nfrac, &sed.phir );
        sed.phir = tan( sed.phir * PI/180.0 );
        break;

      case kSED_LOADEQ:
        sscanf( textLine, "$SED_LOADEQ %d", &sed.loadeq );
        break;

      case kSED_LS:
        sscanf( textLine, "$SED_LS %d %lf %lf %lf", &sed.lsType,
                          &sed.minLs, &sed.factLs, &sed.alfaLs );
        break;

      case kSED_SLOPE:
        sscanf( textLine, "$SED_SLOPE %d %lf %lf %lf %lf %lf", &sed.slope,
                &sed.alfaSlope, &sed.betaSlope, &sed.gammaSlope, &sed.deltaSlope,
                &sed.maxSlope );
        sed.maxSlope = tan( sed.maxSlope * PI/180.0 );
        break;

      case kSED_MINQB:
        sscanf( textLine, "$SED_MINQB %lf", &sed.minQb );
        break;

      case kSED_MAXDZ:
        sscanf( textLine, "$SED_MAXDZ %lf", &sed.maxDz );
        break;

      case kSED_EXNEREQ:
        sscanf( textLine, "$SED_EXNEREQ %d", &sed.exnerEq );
        break;

      case kSED_ZB_INIT:
        sscanf( textLine, "$SED_ZB_INIT %d", &sed.zb_init );
        break;

      // ---------------------------------------------------------------------------------
      case kSCALE:
        {
          int    typeOfScale = 0;
          double lScale = 1.0;
          double hScale = 1.0;
          double kScale = 1.0;
          double vScale = 1.0;
          sscanf( textLine, "$SCALE %d", &(typeOfScale) );

          switch( typeOfScale )
          {
            case 1:
              sscanf( textLine, "$SCALE %d %lf %lf", &(typeOfScale),
                                                     &(lScale),
                                                     &(vScale) );
              scale.init( typeOfScale, lScale, vScale );
              break;

            case 2:
              sscanf( textLine, "$SCALE %d %lf %lf %lf", &(typeOfScale),
                                                         &(lScale),
                                                         &(hScale),
                                                         &(kScale) );
              scale.init( typeOfScale, lScale, hScale, kScale );
              break;
          }
        }
        break;
    }
  }

  delete file;


  // if minVtKD or minVtC not present, use minimum value minVt ---------------------------

  if( _minVt )
  {
    if( !_minVtKD )  minVtKD = dep_minVt;
    if( !_minVtC )   minVtC  = dep_minVt;
  }

  if( !_minVtxx )  dep_minVtxx = dep_minVt;
  if( !_minVtyy )  dep_minVtyy = dep_minVt;


  // -------------------------------------------------------------------------------------
  // check and output file names

  if( name.regionFile[0] == ' ' )
  {
    REPORT::rpt.Error( kUserFault, "no region file specified (PROJECT::Input #17)" );
  }

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n  %30s   %s\n  %30s   %s\n",
                          " input...... time step file:",   name.timeStepFile,
                          "              material file:",   name.materialFile,
                          "                region file:",   name.regionFile,
                          "               control file:",   name.controlFile );

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n  %30s   %s\n  %30s   %s\n\n",
                          "               initial file:",   name.initialFile,
                          "             statistic file:",   name.inputStatistFile,
                          "             subdomain file:",   name.subdomFile,
                          "               section file:",   name.sectionFile );

  REPORT::rpt.Message( 0, "  %30s   %s / %d\n  %30s   %s\n",
                          "output......... report file:",   name.reportFile, REPORT::rpt.level,
                          "              geometry file:",   name.geometryFile );

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n",
                          "               restart file:",   name.restartFile,
                          "             statistic file:",   name.outputStatistFile );

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n",
                          "                region file:",   name.rgAvsFile,
                          "               control file:",   name.ctAvsFile);

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n",
                          "               reduced file:",   name.redAvsFile,
                          "                   wet file:",   name.wetAvsFile);

  REPORT::rpt.Message( 0, "  %30s   %s\n  %30s   %s\n\n",
                          "      statistic region file:",   name.rgAvsStatistFile,
                          "       optional output file:",   name.optOutputFile );

  REPORT::rpt.Message( 0, "  %30s   %s\n",
                          "     name of subdomain path:",  outputPath );

  REPORT::rpt.Line1( 0 );


  // -------------------------------------------------------------------------------------

  if( vpcomp == 0 )           // use default values for UCD-output
  {
    vpdata = 11;
    vpcomp =  9;

    vpoutlist[0] = kUV;
    vpoutlist[1] = kS;
    vpoutlist[2] = kK;
    vpoutlist[3] = kD;
    vpoutlist[4] = kC;
    vpoutlist[5] = kH;
    vpoutlist[6] = kUS;
    vpoutlist[7] = kTAU;
    vpoutlist[8] = kVT;

    vedata = 1;
    vecomp = 1;

    veoutlist[0] = kS;
  }

  REPORT::rpt.Message( 3, "\n   variables for output in AVS-UCD-Files\n" );

  for( int i=0; i<vpcomp; i++ )
  {
    int id = vpoutlist[i];
    REPORT::rpt.Message( 3, "%20s [%s]\n", valist[id].name, valist[id].unit );
  }

  REPORT::rpt.Line1( 3 );


  // -------------------------------------------------------------------------------------
  // report solver specifications

  if( SOLVER::m_neqs == 0 )   // define one default solver
  {
    SOLVER::m_solver = new SOLVER* [1];
    if( !SOLVER::m_solver )
      REPORT::rpt.Error( kMemoryFault, "%s (PROJECT::Input #14)",
                         "can not allocate memory" );

    SOLVER::m_solver[SOLVER::m_neqs] = new P_BCGSTABD();
    if( !SOLVER::m_solver[SOLVER::m_neqs] )
      REPORT::rpt.Error( kMemoryFault, "%s (PROJECT::Input #18)",
                         "can not allocate memory" );

    SOLVER::m_solver[SOLVER::m_neqs]->preconType =     1;
    SOLVER::m_solver[SOLVER::m_neqs]->proceed    =    -1;
    SOLVER::m_solver[SOLVER::m_neqs]->mceq       =    90;
    SOLVER::m_solver[SOLVER::m_neqs]->maxIter    =  1000;
    SOLVER::m_solver[SOLVER::m_neqs]->maxDiff    = 0.010;

    SOLVER::m_neqs++;
  }

  for( int i=0; i<SOLVER::m_neqs; i++ )
  {
    switch( SOLVER::m_solver[i]->solverType )
    {
      case kFront:
        sprintf( text, "\n   %d. %s\n",
                 i+1, "solver specification: frontal solver" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %s\n",
                 "mfw:",   SOLVER::m_solver[i]->mfw,
                 "size:",  SOLVER::m_solver[i]->size,
                 "path:",  SOLVER::m_solver[i]->path );
        REPORT::rpt.Output( text, 4 );
        break;

      case kFrontm:
        sprintf( text, "\n   %d. %s\n",
                 i+1, "solver specification: frontal matrix solver" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %s\n",
                 "mceq:",  SOLVER::m_solver[i]->mceq,
                 "mfw:",   SOLVER::m_solver[i]->mfw,
                 "size:",  SOLVER::m_solver[i]->size,
                 "path:",  SOLVER::m_solver[i]->path );
        REPORT::rpt.Output( text, 4 );
        break;

      case kBicgstab:
        sprintf( text, "\n   %d. %s\n",
                 i+1, "solver specification: BiCGStab" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %lf\n",
                 "precon type:", SOLVER::m_solver[i]->preconType,
                 "proceed:",     SOLVER::m_solver[i]->proceed,
                 "mceq:",        SOLVER::m_solver[i]->mceq,
                 "maxIter:",     SOLVER::m_solver[i]->maxIter,
                 "maxDiff:",     SOLVER::m_solver[i]->maxDiff );
        REPORT::rpt.Output( text, 4 );
        break;

      case kParmsBcgstabd:
        sprintf( text, "\n   %d. %s\n",
                 i+1, "solver specification: PARMS - BiCGStab" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %lf\n",
                 "precon type:", SOLVER::m_solver[i]->preconType,
                 "proceed:",     SOLVER::m_solver[i]->proceed,
                 "mceq:",        SOLVER::m_solver[i]->mceq,
                 "maxIter:",     SOLVER::m_solver[i]->maxIter,
                 "maxDiff:",     SOLVER::m_solver[i]->maxDiff );
        REPORT::rpt.Output( text, 4 );
        break;

      case kParmsFgmresd:
        sprintf( text, "\n   %d. %s\n",
                 i+1, "solver specification: PARMS - flexible Gmres" );
        REPORT::rpt.Output( text, 4 );

        sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %d\n  %30s  %lf\n",
                 "precon type:", SOLVER::m_solver[i]->preconType,
                 "proceed:",     SOLVER::m_solver[i]->proceed,
                 "mceq:",        SOLVER::m_solver[i]->mceq,
                 "mkyrl:",       SOLVER::m_solver[i]->mkyrl,
                 "maxIter:",     SOLVER::m_solver[i]->maxIter,
                 "maxDiff:",     SOLVER::m_solver[i]->maxDiff );
        REPORT::rpt.Output( text, 4 );
        break;
    }
  }

  REPORT::rpt.OutputLine1( 4 );


  // -------------------------------------------------------------------------------------
  // read element type specifications

  if( name.materialFile[0] )
  {
    sprintf( text, "   reading material file: %s\n", name.materialFile );
    REPORT::rpt.Output( text, 4 );

    TYPE::Import( name.materialFile );
    REPORT::rpt.OutputLine1( 4 );
  }
  else
  {
    TYPE::Import( name.inputFile );
    REPORT::rpt.OutputLine1( 4 );
  }


  // -------------------------------------------------------------------------------------
  // read section coordinates for 1) reordering of elements and
  //                              2) initialization of water surface elevation

  ASCIIFILE* sfl = NULL;
  int release = 30900;

  if( name.sectionFile[0] )
  {
    sfl = new ASCIIFILE( name.sectionFile, "r" );

    if( !sfl || !sfl->getid() )
      REPORT::rpt.Error( kOpenFileFault, "%s %s (PROJECT::Input #19)",
                         "can not open section file", name.sectionFile );

    textLine = sfl->next();

    sscanf( textLine, " $RISMO2D %d", &release );

    REPORT::rpt.Screen( 2, "\n ### release of section file:            %d ###\n\n", release );

    if( release < 30900 )
    {
      sfl->rewind();

      textLine = sfl->nextLine();
      sscanf( textLine, " %d",  &nSection );

      section = new SECTION [nSection];
      if( !section )
        REPORT::rpt.Error( "can not allocate memory (PROJECT::Input #20)" );

      sprintf( text, "   reading %d section coordinates\n", nSection );
      REPORT::rpt.Output( text, 3 );

      for( int i=0; i<nSection; i++ )
      {
        double xs = 0.0;
        double ys = 0.0;
        double xe = 0.0;
        double ye = 0.0;
        double z  = 0.0;

        textLine = sfl->nextLine();
        sscanf( textLine, " %lf %lf %lf %lf %lf", &xs, &ys, &xe, &ye, &z );

        section[i].init( xs, ys, xe, ye, z );

        sprintf( text, "  %3d:  %12.3lf  %12.3lf  %12.3lf  %12.3lf    %8.3lf\n",
                       i+1, section[i].getxs(), section[i].getys(),
                            section[i].getxe(), section[i].getye(),
                            section[i].z );
        REPORT::rpt.Output( text, 5 );
      }
    }
  }
  else
  {
    sfl = new ASCIIFILE( name.inputFile, "r" );

    if( !sfl || !sfl->getid() )
      REPORT::rpt.Error( kOpenFileFault, "%s %s (PROJECT::Input #21)",
                         "can not open input file", name.inputFile );
  }


  // check version of section file

  if( release >= 30900 )
  {
    for( int rd=0; rd<2; rd++ )  // read the file two times
    {
      while( !feof(sfl->getid()) )
      {
        if( !(textLine = sfl->nextLine()) )  break;

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
      }

      sfl->rewind();
    }
  }

  REPORT::rpt.OutputLine1( 3 );

  delete sfl;


  // report convergence limits -----------------------------------------------------------

  sprintf( text, "  %30s  %9.2le\n  %30s  %9.2le\n  %30s  %9.2le\n",
                 "convergence limits - UV:",  convUV,
                 "                      S:",  convS,
                 "                     KD:",  convKD );
  REPORT::rpt.Output( text, 3 );
  REPORT::rpt.OutputLine1( 3 );


  // report min and max value for node parameters ----------------------------------------

  sprintf( text, "  %30s  %le\n  %30s  %le\n  %30s  %le\n",
                 "vt    - minimum value (depr.):",  dep_minVt,
                 "vt_KD - minimum value        :",  minVtKD,
                 "vt_C  - minimum value        :",  minVtC );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %le\n  %30s  %le\n  %30s  %le\n",
                 "    K - minimum value:",  minK,
                 "    D - minimum value:",  minD,
                 "    C - minimum value:",  minC );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %le\n\n",
                 "   Us - maximum value:",  maxUs );
  REPORT::rpt.Output( text, 3 );
  REPORT::rpt.OutputLine1( 3 );


  // report relaxation parameter ---------------------------------------------------------

  sprintf( text, "  %30s  %d\n  %30s  %9.6lf\n  %30s  %9.6lf\n\n",
                 "relaxation...  method:",  relaxMethod,
                 "              minimum:",  relaxMin,
                 "              maximum:",  relaxMax );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.6lf\n  %30s  %9.6lf\n  %30s  %9.6lf\n\n",
                 "     max change of UV:",  maxDeltaUV,
                 "     max change of S :",  maxDeltaS,
                 "     max change of KD:",  maxDeltaKD );
  REPORT::rpt.Output( text, 3 );
  REPORT::rpt.OutputLine1( 3 );


  // report parameters for dry/rewet algorithm -------------------------------------------

  DRYREW* dryRew = &M2D->region->dryRew;

  sprintf( text, "  %30s  %4d\n  %30s  %4d\n  %30s  %4d\n  %30s  %4d\n\n",
                 "dry/rewet method:",         dryRew->method,
                 "frequency of dry/rewet:",   dryRew->dryRewFreq,
                 "number of rewet passes:",   dryRew->rewetPasses,
                 "count down for rewetting:", dryRew->countDown );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.6lf\n  %30s  %9.6lf\n\n",
                 "limit for drying:",    dryRew->dryLimit,
                 "limit for rewetting:", dryRew->rewetLimit );
  REPORT::rpt.Output( text, 3 );
  REPORT::rpt.OutputLine1( 3 );

  hmin = dryRew->dryLimit / 2.0;


  // report number of smoothing passes ---------------------------------------------------

  sprintf( text, "  %30s  %-d\n  %30s  %-d\n  %30s  %-d\n\n",
                 "BC smoothing passes:",  smoothPassesBC,
                 "KD smoothing passes:",  smoothPassesKD,
                 "vt smoothing passes:",  smoothPassesVT );
  REPORT::rpt.Output( text, 3 );
  REPORT::rpt.OutputLine1( 3 );


  // report type of boundary model for different boundaries ------------------------------

  sprintf( text, "  %30s  %d\n  %30s  %d\n  %30s  %d\n\n",
                 "KD boundary model for inflow:", KDBcon[kIndInflow],
                 "                     outflow:", KDBcon[kIndOutflow],
                 "                        side:", KDBcon[kIndSide] );
  REPORT::rpt.Output( text, 3 );
  REPORT::rpt.OutputLine1( 3 );


  // report degree of GAUSS point integration --------------------------------------------

  sprintf( text, "  %30s  %4d\n\n", "degree of integration:", GPdeg );
  REPORT::rpt.Output( text, 3 );
  REPORT::rpt.OutputLine1( 3 );

  SHAPE::init( GPdeg );


  // -------------------------------------------------------------------------------------
  // report constants
  // kinematic viscosity, density, gravity, von Karman's constant and  wall distance

  sprintf( text, "  %30s  %9.2le\n  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "kinematic viscosity:",   vk,
          "density:",               rho,
          "gravity:",               g );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.4lf\n",
          "von Karman's constant:", kappa );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.7lf\n\n",
          "earth rotation velocity:", earthVel );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // report k-epsilon constants ----------------------------------------------------------

  sprintf(text, "  %30s  %9.4lf\n  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "k-epsilon constants:  cm:",  KD.cm,
          "                      cd:",  KD.cd,
          "                      sK:",  KD.sK );
  REPORT::rpt.Output( text, 3 );

  sprintf(text, "  %30s  %9.4lf\n  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "                      sD:",  KD.sD,
          "                     c1D:",  KD.c1D,
          "                     c2D:",  KD.c2D );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );


  // report scaling factors ----------------------------------------------------------------

  sprintf( text, "  %30s  %9d\n  %30s  %9.4lf\n",
                 "scaling         :",  scale.gettype(),
                 "          lScale:",  scale.getlScale() );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "                 hScale:",  scale.gethScale(),
          "                ksScale:",  scale.getkScale() );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "  %30s  %9.4lf\n  %30s  %9.4lf\n",
          "                 tScale:",  scale.gettScale(),
          "                 vScale:",  scale.getvScale() );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.OutputLine1( 3 );
}
