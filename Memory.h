// =======================================================================================
//                                     M E M O R Y
// ======================================================================================
// This class implements functionality to allocate memory for 1- and 2D arrays.
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
//    date               short description
// ----------   ------   -----------------------------------------------------------------
// 18.09.2004     sc     first implementation
//
// =======================================================================================

// read this header file only once ...

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
