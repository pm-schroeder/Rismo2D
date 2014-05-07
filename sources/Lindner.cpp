// ======================================================================================
//
// Copyright (C) 1992-2012  by  P.M. SCHROEDER
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
// Please feel free to contact me:   Dr. P.M. Schroeder
//                                   Sperberweg 22-1
//                                   75045 Walzbachtal
//                                   Germany
//
//                                   e-mail: michael.schroeder@hnware.de
//
// ======================================================================================
//
// --------------------------------------------------------------------------------------
// Compute friction coefficient for non-submerged vegetation from parameters
//           dp = diameter of vegetation
//           sp = spacing of vegetation
//           va = velocity
//           ha = flow depth
//           cf = friction coefficient for bottom roughness
//           vk = kinemtic viscosity
//           g  = gravity acceleration
//
// The algorithm was developed by LINDNER (1982) and PASCHE (1986).
// I have made some changes to the computation of depth ratio (1990)
// and to the wake length equation (1992): the slope of energy line
// is no longer requested (this essential for use in 2D methods).
//
// Michael Schroeder in November, 1992
// --------------------------------------------------------------------------------------

#include "Defs.h"
#include "Type.h"

#define kMaxIter    200
#define kPrecision  1.0e-3

int TYPE::cWRav    = 0;
int TYPE::cWRmax   = 0;
int TYPE::cWRcount = 0;
int TYPE::aNLav    = 0;
int TYPE::aNLmax   = 0;
int TYPE::aNLcount = 0;
int TYPE::itErr    = 0;

double TYPE::_cWR  = 0.0;      // drag coefficient of vegetation
double TYPE::_aNL  = 0.0;      // wake length
double TYPE::_Ust  = 0.0;      // ratio of flow velocity
double TYPE::_Hst  = 0.0;      // ratio of flow depth


double TYPE::lindner( double dp,
                      double sp,
                      double va,
                      double ha,
                      double cf,
                      double vk,
                      double g )
{
  int i;
  int icWR,               // iteration counter: cWR
      iaNL;               // iteration counter: aNL
  int realRoots;

  double cW, cWR, RcWR, dcWR, aNL, aNB, RaNL;
  double cWR1, cWR2, dcWR1, dcWR2;
  double lambda, Fr;
  double x[3], vRatio, hRatio;
  double alfa, aCof, bCof, cCof, dCof;

  if(     dp < 1.0e-4  ||  sp < 1.0e-4
      ||  va < 1.0e-3  ||  ha < 1.0e-3  )  return cf;

  // initializations
  cWR = 1.0;              // drag coefficient
  aNL = sp/2.0;           // wake length of a cylinder

  cWR1  = cWR2  = 1.0;
  dcWR1 = dcWR2 = 0.0;


  // ------------------------------------------------------------------------------------
  // start of iteration for cWR

  for( icWR=0; icWR<kMaxIter; icWR++ )
  {
    // superposed friction coefficient
    lambda = 8.0 * cf  +  4.0 * cWR * ha * dp / sp / sp;

    // drag coefficient cW for one cylinder
    cW = dragCoeff(va, dp, vk);

    // wake length of a cylinder (iterative computation)
    for( iaNL=0; iaNL<kMaxIter; iaNL++ )
    {
      double tmp1, tmp2;

      tmp1 = 1.0  +  aNL * lambda / 4.0 / ha;
      tmp2 = 30.0 / pow(fabs(tmp1), 1.5);
      RaNL = cW * dp * pow(fabs(tmp2), 1.429);

      // test for convergence
      if( fabs((RaNL-aNL)/RaNL) < kPrecision )
      {
        aNL   = RaNL;
        iaNL *= -1;
        break;
      }

      // revised wake length
      aNL = 0.5 * (RaNL + aNL);
    }

    // ----------------------------------------------------------------------------------
    // statistics of cWR iteration

    // if( iaNL > 0 )  return -1.0;  // number of iterations exceeded

    if( iaNL > 0 )
    {
      aNL = sp/2.0;
    }

    else
    {
      iaNL = abs(iaNL);
      aNLcount++;
      aNLav += iaNL;
      if ( iaNL > aNLmax ) aNLmax = iaNL;
    }

    // ----------------------------------------------------------------------------------

    // wake width
    aNB = 0.24 * pow(fabs(aNL), 0.59) * pow(fabs(cW*dp), 0.41);

    // ratio of velocity in front of and behind cylinder
    vRatio = 1.151 * pow (fabs(aNL/sp), -0.483)
            +  0.5 * pow (fabs(aNB/sp), 1.1);

    // ratio of flow depth
    Fr = va / sqrt ( g * ha );           // Froude number

    alfa = dp / sp;

    aCof =  Fr * Fr * (1.0 - alfa * cWR/2.0);
    bCof = -Fr * Fr - (1.0 - alfa) / 2.0;
    cCof =  0.0;
    dCof = (1.0 - alfa) / 2.0;

    hRatio = 1.0;

    if( fabs(aCof) < 1.0e-10 )
    {
      hRatio = sqrt( -dCof / bCof);
    }
    else
    {
      realRoots = cubeEquation(aCof, bCof, cCof, dCof, x);

      for( i=0; i<realRoots; i++ )
      {
        if( x[i] > 0.0  &&  x[i] < 1.0 )
        {
          hRatio = x[i];
          break;
        }
      }
    }


    // revise drag coefficient cWR
    RcWR = 1.3124 * cW * vRatio  +  2.0 * (1.0 - hRatio) / Fr / Fr;

    // test for convergence
    if( fabs((RcWR-cWR)/RcWR) < kPrecision )
    {
      icWR *= -1;
      break;
    }

    // ----------------------------------------------------------------------------------
    // use PEGASUS algorithm for cWR iteration

    dcWR = RcWR - cWR;

    if( icWR >= 3  &&  dcWR1*dcWR2 < 0.0 )
    {
      if( dcWR2*dcWR < 0.0 )
      {
        dcWR1 *= dcWR2 / (dcWR2+dcWR);
      }
      else
      {
        cWR1  = cWR2;
        dcWR1 = dcWR2;
      }

      cWR2  = cWR;
      dcWR2 = dcWR;
      cWR   = cWR2 - dcWR2 * (cWR2-cWR1) / (dcWR2-dcWR1);
    }
    else
    {
      cWR1 = cWR2;   dcWR1 = dcWR2;
      cWR2 = cWR;    dcWR2 = dcWR;

      if( icWR >= 2  &&  dcWR1*dcWR2 < 0.0 )
      {
        cWR = cWR2 - dcWR2 * (cWR2-cWR1) / (dcWR2-dcWR1);
      }
      else
      {
        cWR = RcWR;
      }
    }
  }

  if( icWR > 0 )
  {
    itErr++;
    return -1.0;
  }

  // ------------------------------------------------------------------------------------
  // statistics of cWR iteration
  icWR = -icWR;
  cWRcount++;
  cWRav += icWR;
  if( icWR > cWRmax ) cWRmax = icWR;

  // ------------------------------------------------------------------------------------
  // return superposed friction coefficient: cf_P + cf_So

  _cWR = cWR;
  _aNL = aNL;
  _Ust = vRatio;
  _Hst = hRatio;

  return lambda / 8.0;
}


double TYPE::dragCoeff( double v, double d, double vk )
{
  double Re;

  Re = v * d / vk;

  if( Re <= 800.0 )   return 3.07 / pow (Re, 0.168);

  if( Re <= 6000.0 )  return ( 1.0 );

  if( Re <= 11000.0 ) return 1.0  +  0.2 * (Re - 6000.0) / 5000.0;

  return 1.2;
}


int TYPE::cubeEquation( double aCof, double bCof, double cCof, double dCof, double x[3] )
{
  int    real;
  double ba, ca, P, Q, Q2P3, U, V;

  ba = bCof / aCof / 3.0;
  ca = cCof / aCof;

  P  = ca/3.0 - ba*ba;
  Q  = ba*ba*ba - ba*ca/2.0 + dCof/aCof/2.0;

  Q2P3 = Q*Q + P*P*P;


  if( Q2P3 > 0.0 )
  {
    double expo, sign, tmp;

    real = 1;

    expo = 1.0/3.0;

    tmp  = -Q + sqrt(Q2P3);
    sign = tmp / fabs(tmp);
    U    = sign * pow (fabs(tmp), expo);

    tmp  = -Q - sqrt(Q2P3);
    sign = tmp / fabs(tmp);
    V    = sign * pow (fabs(tmp), expo);

    x[0] = (U + V) - ba;
  }

  else
  {
    double tmp, phi;

    real = 3;

    tmp = -Q / pow (-P, 1.5);

    if( tmp >= 1.0 )
    {
      phi = 0.0;
    }

    else if( tmp <= -1.0 )
    {
      phi = PI;
    }

    else
    {
      phi = acos (tmp);
    }

    x[0] = 2.0 * sqrt (-P) * cos (phi/3.0)           -  ba;
    x[1] = 2.0 * sqrt (-P) * cos ((phi+2.0*PI)/3.0)  -  ba;
    x[2] = 2.0 * sqrt (-P) * cos ((phi+4.0*PI)/3.0)  -  ba;
  }

  return real;
}


void TYPE::statisLindner()
{
  char text[200];

  if( itErr > 0 )
  {
    sprintf( text, "\n ((TYPE::lindner)        %s %d-times exceeded\n",
                   "limit of iterations", itErr );
    REPORT::rpt.Output( text, 4 );
  }

  if( aNLcount != 0 )
  {
    sprintf( text, "\n (TYPE::lindner)         %s\n",
                   "statistics of iterations:         maximum/average" );
    REPORT::rpt.Output( text, 4 );

    sprintf( text, "%s %8d calls: %2d  /  %2d\n",
                   "                         wake length      aNL",
                   aNLcount, aNLmax, aNLav/aNLcount );
    REPORT::rpt.Output( text, 4 );
  }

  if( cWRcount != 0 )
  {
    sprintf( text, "%s %8d calls: %2d  /  %2d\n",
                   "                         drag coefficient cWR",
                   cWRcount, cWRmax, cWRav/cWRcount );
    REPORT::rpt.Output( text, 4 );
  }

  aNLmax = aNLav = aNLcount = 0;
  cWRmax = cWRav = cWRcount = 0;
  itErr  = 0;
}
