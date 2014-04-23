// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// F R O M A T
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Fromat.h   : definition file of the class.
// Fromat.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the frontal matrix.
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
//  29.06.2005    sc    first implementation / first concept
//  13.07.2005    sc    method MulVec implemented
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FROMAT_INCL
#define FROMAT_INCL

class EQS;
class PROJECT;


class FROW
{
  public:
    unsigned int mark :  1;
    unsigned int no   : 31;

    double* eq;

  public:
    FROW()
    {
      no = 0;
      eq = NULL;
    };

    ~FROW()
    {
    };
};

class EQL
{
  public:
    unsigned int cpl :  1;    // flag: equation is complete
    unsigned int no  : 31;    // equation number

    int  ind;                 // ind >=  0 : pointer to partial eq system m_frow[}
                              // ind == -1 : equation ist not in the partial eq system
    EQL* prev;                // pointered list to remember ordering of equations
    EQL* next;

  public:
    EQL()
    {
      cpl  = false;
      no   = 0;
      ind  = -1;
      next = NULL;
      prev = NULL;
    };

    //////////////////////////////////////////////////////////////////////////////////////

    void Apppend( EQL* eql )  // append "this" after element eql
    {
      this->prev = eql;
      this->next = eql->next;

      if( eql->next )  eql->next->prev = this;
      eql->next = this;
    }

    void Insert( EQL* eql )   // insert "this" before element eql
    {
      this->next = eql;
      this->prev = eql->prev;

      if( eql->prev )  eql->prev->next = this;
      eql->prev = this;
    }

    void Remove()             // remove "this" from the linked list
    {
      if( this->prev )  this->prev->next = this->next;
      if( this->next )  this->next->prev = this->prev;
    }

    EQL* First()              // find topmost element in the linked list
    {
      EQL* first = this;
      while( first->prev )  first = first->prev;
      return first;
    }

    EQL* Last()               // find last element in the linked list
    {
      EQL* last = this;
      while( last->next )  last = last->next;
      return last;
    }
};

class FROMAT
{
  private:
    double*  m_eqbuf;

  public:
    int      m_neq;           // total nuber of equations
    int      m_maxFW;         // max front width (max size of partial eq system)
    int      m_actFW;         // actual front width

    FROW*    m_frow;          // array of length "m_maxFW": holds the partial eq system
    EQL*     m_eql;           // array of length "m_neq":   m_eql[eq number].ind is
                              //                            a pointer to m_frow[]

    double*  m_U;
    double*  m_L;

  public:
    FROMAT( int maxFW, int neq );
    ~FROMAT();

    //void Resize( int maxFW );

    void Insert( int e, int* width, int** index, REALPR** A );
    void Add( int e, int cnt, int* eqno, double* eq );

    void Eliminate( int e );
    void Eliminate( int e, double* B );

    void MulVec( double* x, double* r, EQS* eqs, SUBDOM* subdom );

    void Output();
};

#endif
