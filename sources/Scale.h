// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// S C A L E
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Scale.h   : definition file of the class.
// Scale.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the scaling of model data.
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
//  01.01.1998    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SCALE_INCL
#define SCALE_INCL

class  BCONSET;
class  DRYREW;
class  NODE;
class  PROJECT;
class  SECTION;
class  TIMEINT;
class  TYPE;


class SCALE
{
  private:
    int    type;               // 0 = no scaling
                               // 1 = Froude model
    double lScale;             // scaling of length
    double hScale;             // scaling of height
    double kScale;             // scaling of roughness height ks

    double tScale;             // depending: scaling of time
    double vScale;             //            scaling of velocity
    double vtScale;            //            scaling of viscosity

  public:
    SCALE();
    ~SCALE();

    int    gettype();
    double getlScale();
    double gethScale();
    double getkScale();
    double gettScale();
    double getvScale();
    double getvtScale();

    void init( int, double, double, double );
    void init( int, double, double );

    void scale( PROJECT* );
};

#endif
