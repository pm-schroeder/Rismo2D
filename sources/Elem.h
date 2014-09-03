// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// E L E M
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Elem.h   : definition file of the class.
// Elem.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements a finite element.
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
//  01.01.1992     sc     first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ELEM_INCL
#define ELEM_INCL

#include "Defs.h"
#include "Node.h"

#define kMaxNodes1D       3   // maximum number of nodes of 1D-element
#define kMaxNodes2D       8   // maximum number of nodes of 2D-element


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

    double   U, Uo;                               // velocity U
    double   V, Vo;                               //          V
    double   P;                                   // pressure
    double   dz;                                  // change of bed elevation

  public:
    ELEM();
    ~ELEM();

    ELEM    operator =( const ELEM& );

    void    Setno( int n )   { this->no = n; }
    void    Setname( int n ) { this->name = n; }
    int     Getno()          { return this->no; }
    int     Getname()        { return this->name; }

    void    Setshape( int );
    int     GetshapeID();

    SHAPE*  GetLShape();
    SHAPE*  GetQShape();
    SHAPE*  GetBShape();

    int     Getncn();
    int     Getnnd();
    NODE*   Getnode( int );

    void    Center( double*, double* );
    void    Center( double*, double*, double* );
    double  area();
    double  perimeter();
};

#endif
