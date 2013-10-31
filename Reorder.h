// ======================================================================================
//                                    R E O R D E R
// ======================================================================================
// This class implements the reordering of Finite Element Grids.
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
// 01.01.1994     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once ...

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
