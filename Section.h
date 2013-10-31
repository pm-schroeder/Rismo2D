// ======================================================================================
//                                    S E C T I O N
// ======================================================================================
// This class implements sections.
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

#ifndef SECTION_INCL
#define SECTION_INCL

class ELEM;

class SECTION
{
  private:
    double xs, ys;
    double xe, ye;

  public:
    double tx, ty;
    double nx, ny;
    double z;

    ELEM   *elem;

  public:
    SECTION();
    ~SECTION();

    void   init( double, double, double, double, double );

    double getxs();
    double getys();
    double getxe();
    double getye();

    double length();
    double normalDistance( double, double );
    double tangentialDistance( double, double );
    int    interpolate( SECTION*, double, double, double* );

    void   scaleFR( double, double );
};


#endif

