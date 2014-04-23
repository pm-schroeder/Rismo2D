// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// R E O R D E R
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Reorder.h   : definition file of the class.
// Reorder.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the reordering of Finite Element Grids.
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
//  01.01.1994    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef REORDER_INCL
#define REORDER_INCL


class NODE;
class ELEM;
class GRID;


class REOFRONT
{
  public:
    REOFRONT* next;
    NODE*     node;
    int       hist;
    int       miss;
};

class REORDER
{
  private:
    int    nel;            // Anzahl der verbleibenden Elemente

    REOFRONT* frontPtr;    // FRONT array for each node

    REOFRONT* frnt;        // Liste der aktiven Knoten
    ELEM*  prev_el;        // letztes erfasstes Element
    int    hist_cnt;       // Zaehler fuer Fronthistory
    int    num;            // Anzahl der aktiven Knoten

  public:
    int    minDepth;       // minimale Rekursionstiefe
    int    maxDepth;       // maximale Rekursionstiefe

    int    max;            // maximale Frontweite, Rueckgabewert
    double square;         // Quadratsumme der Frontweiten

  private:
    REORDER();

  public:
    REORDER( GRID* );
    ~REORDER();

    void initFront( ELEM* );
    void start();
    int  frontListOfElems( REOFRONT*, ELEM** );
    int  updateFront( ELEM*, REOFRONT**, int , int );
    int  removeFront( ELEM*, REOFRONT**, int );
    int  frontWidth( REOFRONT**, int, int, double* );
    ELEM *chooseElem( ELEM**, int );
};

#endif
