// ======================================================================================
//
// Copyright (C) 1992-2010  by  P.M. SCHROEDER
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

#include "Defs.h"
#include "Node.h"
#include "Elem.h"

#include "Bcon.h"

#include "Report.h"


DATKEY BCON::datkey[] =
{
    kNone,      "NONE",

    kInlet,     "INLET",          //   0: inlet boundary
    kOutlet,    "OUTLET",         //   1: outlet: water elevation specified
    kOpenBnd,   "OPEN",           //   2: open boundary; no boundary conditions

    kQInlet,    "QINLET",         //   3: discharge along control line (inlet)
    kQTInlet,   "QTINLET",        //   4: Q(t)-relation for inlet

    kSQOutlet,  "SQOUTLET",       //   5: S(Q)-relation for outlet
    kSTOutlet,  "STOUTLET",       //   6: S(t)-relation for outlet

    kSlip,      "SLIP",           //   7: slip velocity; flow direction
    kLoglaw,    "LOGLAW",         //   8: slip velocity; log law

    kSetUV,     "SET_UV",         //   9: Dirichlet velocity (U,V)
    kSetS,      "SET_S",          //  10: Dirichlet flow depth h

    kSource,    "SOURCE",         //  11: source/sink Q [m3/s]

    kSetKD,     "SET_KD",         //  12: Dirichlet K and D

    kSetC,      "SET_C",          //  13: Dirichlet sediment concentration C [kg/m3]
    kRateC,     "RATE_C",         //  14: Neumann; sediment transport rate [m3/s]
};

int BCON::nkey = 16;


BCON::BCON()
{
  kind = 0;
  mark = 0;
  rot  = 0;

  no   = 0;

  nx   = 0.0;
  ny   = 0.0;
  niox = 0.0;
  nioy = 0.0;

  val  = NULL;
}


BCON::~BCON()
{
}


// ---------------------------------------------------------------------------------------
// check for consistence of boundary conditions

int BCON::Consistent( int ki )
{
  switch( ki )
  {
    case BCON::kInlet:
      if( isFS(kind,BCON::kOutlet) )    break;
      if( isFS(kind,BCON::kSlip) )      break;
      if( isFS(kind,BCON::kSetUV) )     break;
      return true;

    case BCON::kInlet|BCON::kQInlet:
      if( isFS(kind,BCON::kOutlet) )    break;
      if( isFS(kind,BCON::kSlip) )      break;
      if( isFS(kind,BCON::kSetUV) )     break;
      return true;

    case BCON::kInlet|BCON::kQInlet|BCON::kQTInlet:
      if( isFS(kind,BCON::kOutlet) )    break;
      if( isFS(kind,BCON::kSlip) )      break;
      if( isFS(kind,BCON::kSetUV) )     break;
      return true;

    case BCON::kOpenBnd:
      if( isFS(kind,BCON::kInlet) )     break;
      if( isFS(kind,BCON::kQInlet) )    break;
      if( isFS(kind,BCON::kQTInlet) )   break;
      if( isFS(kind,BCON::kSource) )    break;
      return true;

    case BCON::kOutlet:
      if( isFS(kind,BCON::kInlet) )     break;
      if( isFS(kind,BCON::kQInlet) )    break;
      if( isFS(kind,BCON::kQTInlet) )   break;
      if( isFS(kind,BCON::kSource) )    break;
      return true;

    case BCON::kOutlet|BCON::kSQOutlet:
      if( isFS(kind,BCON::kInlet) )     break;
      if( isFS(kind,BCON::kQInlet) )    break;
      if( isFS(kind,BCON::kQTInlet) )   break;
      if( isFS(kind,BCON::kSource) )    break;
      return true;

    case BCON::kOutlet|BCON::kSTOutlet:
      if( isFS(kind,BCON::kInlet) )     break;
      if( isFS(kind,BCON::kQInlet) )    break;
      if( isFS(kind,BCON::kQTInlet) )   break;
      if( isFS(kind,BCON::kSource) )    break;
      return true;

    case BCON::kSlip:
      if( isFS(kind,BCON::kInlet) )     break;
      if( isFS(kind,BCON::kQInlet) )    break;
      if( isFS(kind,BCON::kQTInlet) )   break;
      if( isFS(kind,BCON::kSetUV) )     break;
      return true;

    case BCON::kLoglaw:
      return true;

    case BCON::kSetUV:
      if( isFS(kind,BCON::kInlet) )     break;
      if( isFS(kind,BCON::kQInlet) )    break;
      if( isFS(kind,BCON::kQTInlet) )   break;
      if( isFS(kind,BCON::kSlip) )      break;
      return true;

    case BCON::kSource:
      if( isFS(kind,BCON::kOutlet) )    break;
      if( isFS(kind,BCON::kSQOutlet) )  break;
      if( isFS(kind,BCON::kSTOutlet) )  break;
      if( isFS(kind,BCON::kSetS) )      break;
      return true;

    case BCON::kSetS:
      if( isFS(kind,BCON::kSource) )    break;
      return true;

    case BCON::kSetKD:
      return true;

    case BCON::kSetC:
      return true;

    case BCON::kRateC:
      return true;
  }

  // -------------------------------------------------------------------------------------
  // print warning for inconsistent boundary condition

  int bc[2] = { 0, 0 };

  for( int i=0; i<nkey; i++ )
  {
    if( isFS(ki,   datkey[i].id) )  bc[0] = i;
    if( isFS(kind, datkey[i].id) )  bc[1] = i;
  }

  REPORT::rpt.Warning( kParameterFault,
                       "%s %s + %s %s %d\n",
                       "inconsistent boundary conditions",
                       datkey[bc[0]].name, datkey[bc[1]].name, "at node", no+1 );
  return false;
}


// ---------------------------------------------------------------------------------------
// determine element (r,c) of rotation matrix

double BCON::Getrot( int r, int c )
{
  switch( rot )
  {
    case kTangentFlowRot:   // rotation matrix on closed boundaries
           if( r == 0 && c == 0 )  return  ny;
      else if( r == 0 && c == 1 )  return  nx;
      else if( r == 1 && c == 0 )  return -nx;
      else if( r == 1 && c == 1 )  return  ny;
      else                         return 0.0;
      break;

    case kNormalFlowRot:    // rotation matrix on open boundaries (inlet & outlet)
           if( r + c == 0 )  return  niox;
      else if( r + c == 1 )  return  nioy;
      else if( r + c == 2 )  return -niox;
      else                   return   0.0;
      break;

    case kSlipFlowRot:      // specified slip flow boundary
           if( r + c == 0 )  return  val->U;
      else if( r + c == 1 )  return  val->V;
      else if( r + c == 2 )  return -val->U;
      else                   return     0.0;
      break;

    default:
      return 0.0;
  }
}
