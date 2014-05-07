// ======================================================================================
//                                      S C A L E
// ======================================================================================
// This class implements the scaling of model data.
// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
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
//
// ======================================================================================

// read this header file only once ...

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
