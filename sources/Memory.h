// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// M E M O R Y
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Memory.h   : definition file of the class.
// Memory.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements functionality to allocate memory for 1D- and 2D-arrays.
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
//  18.09.2004    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MEMORY_INCL
#define MEMORY_INCL


// =======================================================================================
// Includes

#include "Defs.h"


class MEMORY
{
  // =====================================================================================
  //                               V A R I A B L E S
  // =====================================================================================

  public:
    static MEMORY memo;

  protected:

  private:
    union ITEM
    {
      int    i;
      double d;
      void*  v;
    };

    int           m_array;

    ITEM**        m_temp;     // teporary used arrays
    unsigned int* m_size;     // size of arrays
    unsigned int* m_flag;     // flag to mark arrays ...
    enum {
      kUsed    =  1,          // ... in use
      kDelete  =  2,          // ... can be deleted
      kNd      =  4,
      kEl      =  8,
      kEq      = 16,
      kImatrix = 32,          // ... 2-dimensional int array
      kDmatrix = 64           // ... 2-dimensional double array
    };
                              // maximum size of arrays ...
    unsigned int m_max_nnd;   // length of list: GRID::np
    unsigned int m_max_nel;   //                 GRID::ne
    unsigned int m_max_neq;   //                 EQS::neq


  // =====================================================================================
  //                                 M E T H O D S
  // =====================================================================================

  public:
    // -----------------------------------------------------------------------------------
    // constructor
    MEMORY( int array =200 );

    // -----------------------------------------------------------------------------------
    // destructor
    ~MEMORY();

    // -----------------------------------------------------------------------------------
    // return a pointer to an available array with size m_nd, m_el or m_eq ...

    void*    Array_nd( unsigned int nnd );
    void*    Array_el( unsigned int nel );
    void*    Array_eq( unsigned int neq );

    void*    Array( unsigned int n, unsigned int flag =0 );

    void     Detach( void* ptr );
    void     Delete( void* ptr );

    //char*    Carray( unsigned int n );
    //int*     Iarray( unsigned int n );
    //double*  Darray( unsigned int n );
    //void*    Parray( unsigned int n );

    int**    Imatrix( unsigned int rows, unsigned int cols );
    double** Dmatrix( unsigned int rows, unsigned int cols );

    void     PrintInfo();

  // =====================================================================================
  private:
    void     Delete( int i );

  // =====================================================================================
  protected:
};

#endif
