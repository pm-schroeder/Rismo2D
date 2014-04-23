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
#include "Defs.h"
#include "Type.h"

double TYPE::_cb = 0.0;        // fricition coefficient of bottom roughness
double TYPE::_cp = 0.0;        // fricition coefficient of non-submerged vegetation

double TYPE::_Hd = 0.0;        // height and
double TYPE::_Ld = 0.0;        // length of dunes
double TYPE::_Hr = 0.0;        // height and
double TYPE::_Lr = 0.0;        // length of ripples

// --------------------------------------------------------------------------------------
// Methode TYPE::bottom()
// --------------------------------------------------------------------------------------

double TYPE::bottom( double  Us,        // scalar velocity
                     double  H,         // flow depth
                     double  ka,        // von Karman's constant ( = 0.41 )
                     double  vk,        // kinematic viscosity
                     double  g,         // gravity acceleration ( = 9.81 )
                     double  rho,       // density of water
                     double  rhob,      // density of sediment
                     double  d50,       // 50% diameter of grain
                     double  d90 )      // 90% diameter of grain
{
  // ------------------------------------------------------------------------------------
  // compute form roughness, test for presence of vegetation
  double kd = 0.0;
  double kr = 0.0;

  if( form > 0  &&  dp[0] < 1.0e-4  &&  dp[1] < 1.0e-4 )
  {
    // compute height and length of ripples and dunes -----------------------------------
    Dune( Us, H, ka, vk, g, rho, rhob, d50, d90, &_Hd, &_Ld, &_Hr, &_Lr );

    if( _Ld > 0.01 )
    {
      kd = duneCoef * _Hd * ( 1.0 - exp(-25.0*_Hd/_Ld) );
//    if( _Lr > 0.01 )  kr = 20.0 * 0.7 * _Hr * _Hr / _Lr;
      if( _Lr > 0.01 )  kr = duneCoef * _Hr * ( 1.0 - exp(-25.0*_Hr/_Lr) );
    }
    else if( _Lr > 0.01 )
    {
      kr = duneCoef * _Hr * ( 1.0 - exp(-25.0*_Hr/_Lr) );
    }
  }

  // ------------------------------------------------------------------------------------
  // compute friction coefficient for bottom roughness
  double cb = 0.0;

  // Colebrook-Whites or Nikuradses law chosen for grain roughness ...
  if( rtype < kCHEZ  &&  kslaw > 0 )
  {
    cb = friction( rtype, ksfact*d90 + kd + kr, Us, H, ka, vk, g );
  }

  // ... Chezys or Mannings law chosen for grain roughness
  else
  {
    cb  = friction( rtype, rcoef, Us, H, ka, vk, g );
    cb += friction( kNIKU, kd + kr, Us, H, ka, vk, g );
  }

  _cb = cb;          // remind cb in static class variable _cb

  // ------------------------------------------------------------------------------------
  // compute friction coefficient for non-submerged vegetation
  double cp = lindner( dp[0], sp[0], Us, H, cb, vk, g );

  if( cp < -0.9 )    // means: cp == -1.0
  {
    // approximately set cwr = 1.5
    cp  =  cb  +  0.75 * H * dp[0] / sp[0] / sp[0];
  }

  _cp = cp;          // remind cp in static class variable _cp

  // ------------------------------------------------------------------------------------
  // compute friction coefficient for submerged vegetation
  double cv = velzen( hp, dp[1], sp[1], Us, H, cb, ka, vk, g );

  // ------------------------------------------------------------------------------------
  // return the assembled friction coefficient for non-submerged and submerged vegetation
  return cp + cv;
}


// --------------------------------------------------------------------------------------
// Methode TYPE::friction()
// --------------------------------------------------------------------------------------

double TYPE::friction( int    rtype,        // roughness type
                       double rcoef,        // roughness parameter depending on rType
                       double Us,           // scalar velocity
                       double h,            // flow depth
                       double ka,           // von Karman's constant ( = 0.41 )
                       double vk,           // kinematic viscosity
                       double g )           // gravity acceleration ( = 9.81 )
{
  int    iter;
  double term;
  double Re;
  double cf = 0.0;

  switch( rtype )
  {
    // ----------------------------------------------------------------------------------
    // friction losses computed with the Colebrook-White formula
    // ----------------------------------------------------------------------------------

    case kCOWI:
      Re = Us * 4.0 * h / vk;

      // the original condition for laminar/turbulent flow
      // could not be hold (problems during NR iteration):

      if( Re > 500.0 )
      {

        iter = 0;

        // -----------------------------------------------------------------------------------------
        // NO FURTHER NEED FOR THIS!
        // test condition: roughness less than flow depth?
        // if( rcoef < 0.5 * h )
        {
          double oldCf;
          double term1 = 4.4 / Re;
          // double term2 = rcoef / 14.84 / h; // changed on 23.11.2012 (sc)
          double term2 = rcoef / (14.84 * h + rcoef);

          cf = 2.5;                            // initialize cf=1/sqrt(cf) for iteration

          do
          {
            oldCf = cf;
            cf    = -2.03 * log10( oldCf*term1 + term2 );
            iter++;

          } while( fabs((cf-oldCf)/cf) > 1.0e-6  &&  iter < 50 );

          if( iter >= 50 )
          {
            cf = -2.03 * log10(term2);
            //REPORT::rpt.Warning( kValueFault, "%s (TYPE::friction #2)",
            //                     "Colebrook-White's roghness law did not converge" );
          }

          cf = 1.0 / cf / cf / 8.0;
        }
        // -----------------------------------------------------------------------------------------
        // NO FURTHER NEED FOR THIS!
        // else
        // {
        //   // changed on 11.02.2008 and 07.10.2009; sc
        //   // compute Manning's n to gain a smooth change of friction for h <= ks
        //   double rc = pow(2.0*rcoef,0.166667) / sqrt(8.0*g) / 2.03 / log10(29.68);
        //   cf = (rc * rc) * g / pow( h, 0.333333 );
        // }
      }

      else if( Re > 100.0 )
      {
        cf = 8.0 / Re;
      }

      else
      {
        cf = 0.08;
      }
      break;


    // ----------------------------------------------------------------------------------
    // friction losses computed with the Nikuradse formula
    // ----------------------------------------------------------------------------------

    case kNIKU:
      term = ka / log( 1.0  +  12.0 * h / rcoef );
      cf   = term * term;

      // -----------------------------------------------------------------------------------------
      // NO FURTHER NEED FOR THIS!
      // if( rcoef < 0.5 * h )
      // {
      //   term = ka / log( 12.0 * h / rcoef );
      //   cf   = term * term;
      // }
      // else
      // {
      //   // changed on 11.02.2008, sc
      //   // compute Manning's n to gain a smooth change of friction for h <= ks
      //   // rc = nDef;                        // default Manning's n
      //   double rc = ka * pow(2.0*rcoef,0.166667) / sqrt(g) / log(24.0);
      //   cf = (rc * rc) * g / pow( h, 0.333333 );
      // }
      break;


    // ----------------------------------------------------------------------------------
    // friction losses computed with the Manning's n
    // ----------------------------------------------------------------------------------

    case kMANN:
      cf = g * (rcoef * rcoef) / pow( h, 0.333333 );
      break;


    // ----------------------------------------------------------------------------------
    // friction losses computed with the Chezy's C
    // ----------------------------------------------------------------------------------

    case kCHEZ:
      cf = g / (rcoef * rcoef);
      break;
  }


  return cf;
}


// ---------------------------------------------------------------------------------------
// compute height and length of dunes and/or ripples
// ---------------------------------------------------------------------------------------

void TYPE::Dune( double  Us,        // flow velocity
                 double  H,         // flow depth
                 double  ka,        // von Karman's constant ( = 0.41 )
                 double  vk,        // kinematic viscosity
                 double  g,         // gravity acceleration ( = 9.81 )
                 double  rho,       // density of water
                 double  rhob,      // density of sediment
                 double  d50,       // 50% diameter of grain
                 double  d90,       // 90% diameter of grain
                 double* Hd,        // height and
                 double* Ld,        // length of dunes
                 double* Hr,        // height and
                 double* Lr )       // length of ripples
{
  double Utau = Us * ka / log( 12.0 * H / d90 );

  // shear stress according to the grain
  double tau = rho * Utau * Utau;

  // relative grain size under buoyancy conditions
  double rr  = rhob/rho - 1.0;

  // sedimentological particle diameter
  double Dst = d50 * pow( rr*g/vk/vk, 0.3333 );

  // critical shields parameter
  double thetacr = 0.055;

  if(                     Dst <=   6.0 ) thetacr = 0.109 * pow( Dst, -0.50 );
  else if(   6.0 < Dst && Dst <=  10.0 ) thetacr = 0.14  * pow( Dst, -0.64 );
  else if(  10.0 < Dst && Dst <=  20.0 ) thetacr = 0.04  * pow( Dst, -0.10 );
  else if(  20.0 < Dst && Dst <= 150.0 ) thetacr = 0.013 * pow( Dst,  0.29 );
  else if( 150.0 < Dst )                 thetacr = 0.055;

  // critical shear stress
  double taucr = thetacr * (rhob - rho) * g * d50;

  // Tranportaparameter
  double Tst = (tau - taucr) / taucr;


  // -------------------------------------------------------------------------------------
  // compute height and length of ripples and dunes

  *Hr = *Lr = 0.0;
  *Hd = *Ld = 0.0;

  switch( form )
  {
    // VAN RIJN (1993) -------------------------------------------------------------------
    case 1:
      if( Tst > 0.0  &&  Dst >= 1.0 )
      {
        // ripples and/or dunes
        if( Dst <= 10.0 )
        {
          if( Tst <= 3.0 )        // mini-ripples
          {
            *Hr = 100.0 * d50;
            *Lr = 700.0 * d50;
          }
          else if( Tst < 10.0 )  // mega-ripples and dunes
          {
            *Hr = H * 0.02 * (1.0 - exp(-0.1*Tst)) * (10.0 - Tst);
            *Lr = 0.5 * H;

            *Hd = H * 0.11 * pow(d50/H, 0.3) * (1.0 - exp(-0.5*Tst)) * (25.0 - Tst);
            *Ld = 7.3 * H;
          }
          else if( Tst < 25.0 )  // only dunes
          {
            *Hd = H * 0.11 * pow(d50/H, 0.3) * (1.0 - exp(-0.5*Tst)) * (25.0 - Tst);
            *Ld = 7.3 * H;
          }
        }

        // dunes
        else if( Tst < 25.0 )
        {
          *Hd = H * 0.11 * pow(d50/H, 0.3) * (1.0 - exp(-0.5*Tst)) * (25.0 - Tst);
          *Ld = 7.3 * H;
        }
      }
      break;

    // YALIN (1980) ----------------------------------------------------------------------
    case 2:
      if( tau > taucr )
      {
        *Hd = 0.023 * Tst * exp(1.0 - Tst/12.84);
        *Ld = 6.3 * H;
      }
      break;
  }
}
