// ======================================================================================
//                                      F R O M A T
// ======================================================================================
// This class implements the frontal matrix.
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
// 29.06.2005     sc     first implementation / first concept
// 13.07.2005     sc     method MulVec implemented
//
// ======================================================================================

// read this header file only once ...

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
