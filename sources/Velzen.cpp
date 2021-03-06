// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class TYPE
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

// ======================================================================================
// Description
// ======================================================================================
//
// Compute friction coefficient for submerged or non-submerged vegetation
// parameters:  hp = height of vegetation
//              dp = diameter of vegetation
//              sp = spacing of vegetation
//              va = velocity
//              ha = flow depth
//              cf = friction coefficient for bottom roughness
//              ka = von Karman's constant ( = 0.41 )
//              vk = kinemtic viscosity
//              g  = gravity acceleration
//
// The method was described by VAN VELZEN et al. (2003). Additionally to the method
// of VAN VELZEN the drag coeefficient cWR is computed with the method TYPE::lindner().
//
// Michael Schroeder in May, 2008
//
// ======================================================================================
// History of changes
// ======================================================================================
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.05.2008     sc     first implementation
// 08.08.2013     sc     run-time optimization
// --------------------------------------------------------------------------------------

#include "Defs.h"
#include "Type.h"

double TYPE::velzen( double hp,
                     double dp,
                     double sp,
                     double va,
                     double ha,
                     double cf,
                     double ka,
                     double vk,
                     double g )
{
  if(     dp < 1.0e-4  ||  sp < 1.0e-4
      ||  hp < 1.0e-4  ||  va < 1.0e-3  ||  ha < 1.0e-3  )  return 0.0;

  // ------------------------------------------------------------------------------------
  // check the flow depth ha against the vegetation height hp
  double cv;

  if( ha > hp )
  {
    // ------------------------------------------------------------------------------------
    // compute friction coefficient with method TYPE::lindner()
    double cp = lindner( dp, sp, va, hp, cf, vk, g );
    if( cp < -0.9 ) cp = cf  +  0.75 * hp * dp / sp / sp;

    // ------------------------------------------------------------------------------------
    // case ha  > hp (submerged) : separate the flow profile into zones according
    //                             to the proposal of VAN VELZEN et al.
    double ho  = ha - hp;
    double hro = ho / ha;
    double hrp = hp / ha;
    double ko  = 1.6 * pow( hp, 0.7 );

    double co  = ka / log( 1.0 + 12.0*ho/ko );
    double coi = sqrt(hrp)/sqrt(cp) + hro*sqrt(hro)/co;

    cv = 1.0 / coi / coi;
  }
  else
  {
    // ------------------------------------------------------------------------------------
    // compute friction coefficient with method TYPE::lindner()
    cv = lindner( dp, sp, va, ha, cf, vk, g );
    if( cv < -0.9 ) cv = cf  +  0.75 * ha * dp / sp / sp;
  }

  // ------------------------------------------------------------------------------------
  // return superposed friction coefficient: cf + cp + cv
  return cv;
}
