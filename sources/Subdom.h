// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// S U B D O M
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Subdom.h   : definition file of the class.
// Subdom.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the domain decomposition for parallel computing.
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
//  27.05.2004    sc    first implementation (domains with overlapping element)
//  16.07.2004    sc    domains without overlapping elements
//  25.09.2004    sc    new inner classes SD_NODE and SD_ELEM
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SUBDOM_INCL
#define SUBDOM_INCL

// =======================================================================================
// include files

#include <stdio.h>
#include "Defs.h"

class GRID;
class NODE;
class ELEM;
class SUB;


class INFACE
{
  public:
    int     exist;

    int     np;                        // number of interface nodes
    NODE**  node;                      // list of interface node pointers

    double* recv;                      // arrays for receiving and sending data
    double* send;

    char*   ria1;                      // integer arrays to receive / send data
    char*   sia1;

    int*    ria2;
    int*    sia2;


  public:
    // constructor -----------------------------------------------------------------------
    INFACE()
    {
      np    = 0;
      exist = false;
    };

    // destructor ------------------------------------------------------------------------
    ~INFACE()
    {
      if( np )
      {
        delete[] node;
        delete[] recv;
        delete[] send;
        delete[] ria1;
        delete[] ria2;
        delete[] sia1;
        delete[] sia2;

        np   = 0;
        node = NULL;
      }
    };
};


class SUBDOM
{
  // =====================================================================================
  //                               V A R I A B L E S
  // =====================================================================================

  public:
    int pid;                            // subdomain-id
    int npr;                            // total number of subdomains

    int np;                             // total number of nodes (size of sdnd[])
    int ne;                             // total number of elements (size of sel[])

    int npdom;                          // number of nodes in the domain
    int nedom;                          // number of elements in the domain

    NODE** node;                        // pointer to nodes in this subdomain or NULL
    ELEM** elem;                        // pointer to elements in this subdomain or NULL

    SUB*   subbuf;                      // buffer to hold pointers from nodes to SUB

    class SD_NODE
    {
      public:
        enum {
          kInface       =       1,      // 1: interface
          kInface_UP    =       2,      // 2: upstream interface
          kInface_DN    =       4       // 3: downstream interface
        };

        unsigned int mark     :  1;     // mark nodes contained in this domain
        unsigned int flag     :  3;     // general flags

        int sub;

      public:
        SD_NODE()
        {
          mark = false;
          flag = 0;
          sub  = 0;
        };
    };
    SD_NODE* sdnd;                      // node list with length of total grid

    class SD_ELEM
    {
      public:
        unsigned int mark     :  1;     // mark elements contained in this domain

        int sub;

      public:
        SD_ELEM()
        {
          mark = false;
          sub  = 0;
        };
    };
    SD_ELEM* sdel;                      // element list with length of total grid

    INFACE*  inface;                    // list of interfaces

  protected:
  private:

  // =====================================================================================
  //                                M E T H O D S
  // =====================================================================================

  public:
    // -----------------------------------------------------------------------------------
    // constructor
    SUBDOM();

    // -----------------------------------------------------------------------------------
    // destructor
    ~SUBDOM();

    // -----------------------------------------------------------------------------------
    void   Input( char* subdomFileName, char* regionFileName, GRID* region );

    void   SetInface( GRID* region );

    int    Mpi_max( int num );
    double Mpi_max( double num );
    double Mpi_maxabs( double num );
    int    Mpi_min( int num );
    double Mpi_min( double num );
    int    Mpi_sum( int data );
    double Mpi_sum( double data );

    void   Mpi_assemble( double* vec );
    void   Mpi_assemble( int* cnt );
    void   Mpi_average( double* vec );
    void   Mpi_max( double* vec );

  // =====================================================================================
  private:

  // =====================================================================================
  protected:
};

#endif
