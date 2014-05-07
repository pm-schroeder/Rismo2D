// ======================================================================================
//                                       E L E M
// ======================================================================================
// This class implements a finite element.
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
// 01.01.1992     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once

#ifndef ELEM_INCL
#define ELEM_INCL

#include "Defs.h"
#include "Node.h"

#define kMaxNodes1D       3   // maximum number of nodes at 1D-element
#define kMaxNodes2D       8   // maximum number of nodes at 2D-element


class SHAPE;
class TYPE;

class ELEM
{
  private:
    int no;
    int name;

  public:
    enum
    {
      kInlet        =          1,       //  1: inflow boundary
      kOutlet       =          2,       //  2: outflow; water elevation specified
      kOpenBnd      =          4,       //  3: open boundary; no boundary conditions
      kSlip         =          8,       //  4: slip velocity; flow direction specified

      kBound        =         16,       //  5: boundary element (1D/2D)
      kRegion       =         32,       //  6: normal element
      kControl      =         64,       //  7: element for continuity check

      kDry          =        128,       //  8: dry flag
      kDryPrev      =        256,       //  9:
      kMarsh        =        512,       // 10:
      kMarshPrev    =       1024        // 11:
    };

    unsigned int flag    : 15;          // flag
    unsigned int mark    :  1;          // mark elements/nodes (temporary used)
    unsigned int rigid   :  1;

  public:
    ELEM*    link;

    int      shape;                               // shape of element (kLine, KTri, kQuad)
    int      type;                                // physical element type (MAT-ID)
    NODE*    nd[kMaxNodes2D];                     // nodes at element
    int      isLast[kSimDF][kMaxNodes2D];         // TRUE if last occurence of equation

    double   areaFact;                            // factor for inclined elements

    double   U;                                   // velocity U
    double   V;                                   //          V
    double   P;                                   // pressure
    double   dz;                                  // change of bed elevation

//    SHAPE*   lShape;                              // linear shape specification
//    SHAPE*   bShape;                              // hyperbolic shape specification
//    SHAPE*   qShape;                              // quadratic shape specification

  public:
    ELEM();
    ~ELEM();

    ELEM    operator =( const ELEM& );

    void    Setno( int n )   { this->no = n; };
    void    Setname( int n ) { this->name = n; };
    int     Getno()          { return this->no; };
    int     Getname()        { return this->name; };

    void    Setshape( int );
    int     GetshapeID();

    SHAPE*  GetLShape();
    SHAPE*  GetQShape();
    SHAPE*  GetBShape();

    int     Getncn();
    int     Getnnd();
    NODE*   Getnode( int );

    void    center( double*, double* );
    double  area();
    double  perimeter();
};

#endif
