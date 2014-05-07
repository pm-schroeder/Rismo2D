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
//////////////////////////////////////////////////////////////////////////////////////////
//
// CRSMat.cpp: implementation of the CRSMAT class.
//
//////////////////////////////////////////////////////////////////////////////////////////
//
#include "Defs.h"
#include "Report.h"
#include "Eqs.h"
#include "Project.h"
#include "Fromat.h"

#include "CRSMat.h"


//////////////////////////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////////////////////////

CRSMAT::CRSMAT()
{
  m_user    = true;
  m_shadow  = false;

  m_neq     = 0;
  m_ceq     = 0;
  m_entries = 0;

  m_width   = NULL;
  m_index   = NULL;
  m_A       = NULL;

  m_nbuf    = 0;
  m_Ibuf    = NULL;
  m_Abuf    = NULL;
}


CRSMAT::CRSMAT( int n, int m, int entries )
{
  m_user    = false;
  m_shadow  = false;

  m_neq     = n;
  m_ceq     = m;

  m_entries = 0;

  m_nbuf    = 0;
  m_bufsz   = n;

  // -------------------------------------------------------------------------------------

  m_width   = new int     [m_neq];
  m_index   = new int*    [m_neq];
  m_A       = new REALPR* [m_neq];

  if( !m_width || !m_index || !m_A )
    REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - CRSMAT::CRSMAT(2)" );

  for( int i=0; i<m_neq; i++ )  m_width[i] = 0;

  // -------------------------------------------------------------------------------------

  if( m_ceq > 0 )
  {
    Alloc_Buf();
    Alloc_I();

    if( entries > 0 )
    {
      for( int i=0; i<m_neq; i++ )  m_width[i] = m_ceq;
      Alloc_A();
    }
  }
  else
  {
    Alloc_Buf( k_nbuf );
  }
}


CRSMAT::CRSMAT( CRSMAT* A )
{
  m_user    = false;
  m_shadow  = true;

  m_neq     = A->m_neq;
  m_ceq     = A->m_ceq;

  m_neq_up  = A->m_neq_up;
  m_neq_dn  = A->m_neq_dn;

  m_width   = A->m_width;
  m_index   = A->m_index;

  m_entries = 0;

  m_A = new REALPR* [m_neq];
  if( !m_A )
    REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - CRSMAT::CRSMAT(3)" );

  Alloc_Buf();
  Alloc_A();
}


CRSMAT::~CRSMAT()
{
  if( !m_user )
  {
    if( !m_shadow )
    {
      if( m_width )  delete[] m_width;
      if( m_index )  delete[] m_index;

      for( int i=0; i<m_nbuf; i++)  if( m_Ibuf[i] )  delete[] m_Ibuf[i];
      delete[] m_Ibuf;
    }

    if( m_A )  delete[] m_A;

    for( int i=0; i<m_nbuf; i++)  if( m_Abuf[i] )  delete[] m_Abuf[i];
    delete[] m_Abuf;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////
// methods
//////////////////////////////////////////////////////////////////////////////////////////

void CRSMAT::Alloc_Buf( int nb )
{
  m_nbuf = nb;
  if( m_nbuf <= 0 )  m_nbuf = 1;

  m_size = new long    [m_nbuf];
  m_Ibuf = new int*    [m_nbuf];
  m_Abuf = new REALPR* [m_nbuf];

  if( !m_size || !m_Ibuf || !m_Abuf )
    REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - CRSMAT::Alloc_Buf(1)" );

  for( int i=0; i<m_nbuf; i++ )
  {
    m_size[i] = 0;
    m_Ibuf[i] = NULL;
    m_Abuf[i] = NULL;
  }
}

void CRSMAT::Alloc_I()
{
  m_Ibuf[0] = new int [m_neq * m_ceq];

  if( !m_Ibuf[0] )
    REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - CRSMAT::Alloc_I(1)" );

  m_index[0] = m_Ibuf[0];
  for( int i=1; i<m_neq; i++ )  m_index[i] = m_index[i-1] + m_ceq;
}

void CRSMAT::Alloc_A()
{
  m_entries = 0;
  for( int i=0; i<m_neq; i++ )  m_entries += m_width[i];

  m_Abuf[0] = new REALPR [m_entries];

  if( !m_Abuf[0] )
    REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - CRSMAT::Alloc_A(2)" );

  m_A[0] = m_Abuf[0];
  for( int i=1; i<m_neq; i++ )  m_A[i] = m_A[i-1] + m_width[i-1];
}


void CRSMAT::Init()
{
  for( int i=0; i<m_nbuf; i++ )
  {
    for( int j=0; j<m_entries; j++ )  m_Abuf[i][j] = 0.0;
  }
}


void CRSMAT::Append( int eqno, double pivot, FROMAT* fromat )
{
  int width = fromat->m_actFW + 1;

  m_width[eqno] = 0;

  for( int i=0; i<m_nbuf; i++ )
  {
    if( m_size[i] >= width )
    {
      m_index[eqno] = m_Ibuf[i] + m_bufsz - m_size[i];
      m_A[eqno]     = m_Abuf[i] + m_bufsz - m_size[i];
      m_size[i]    -= width;

      m_width[eqno] = width;
      break;
    }
  }

  while( m_width[eqno] <= 0 )
  {
    for( int i=0; i<m_nbuf; i++ )
    {
      if( !m_Ibuf[i] )
      {
        m_size[i] = m_bufsz;

        m_Ibuf[i] = new int    [m_bufsz];
        m_Abuf[i] = new REALPR [m_bufsz];

        if( !m_Ibuf[i] || !m_Abuf[i] )
          REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - CRSMAT::Append(1)" );

        m_index[eqno] = m_Ibuf[i];
        m_A[eqno]     = m_Abuf[i];
        m_size[i]    -= width;

        m_width[eqno] = width;
        break;
      }
    }

    if( m_width[eqno] <= 0 )
    {
      long*    size = new long    [m_nbuf + k_nbuf];
      int**    Ibuf = new int*    [m_nbuf + k_nbuf];
      REALPR** Abuf = new REALPR* [m_nbuf + k_nbuf];

      if( !Ibuf || !Abuf )
        REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - CRSMAT::Append(2)" );

      for( int i=0; i<m_nbuf; i++ )
      {
        size[i] = m_size[i];
        Ibuf[i] = m_Ibuf[i];
        Abuf[i] = m_Abuf[i];
      }

      for( int i=m_nbuf; i<m_nbuf+k_nbuf; i++ )
      {
        size[i] = 0;
        Ibuf[i] = NULL;
        Abuf[i] = NULL;
      }

      if( m_nbuf )
      {
        delete[] m_size;
        delete[] m_Ibuf;
        delete[] m_Abuf;
      }

      m_size = size;
      m_Ibuf = Ibuf;
      m_Abuf = Abuf;

      m_nbuf += k_nbuf;
    }
  }

  // -------------------------------------------------------------------------------------

  m_index[eqno][0] = eqno;
  m_A[eqno][0]     = (REALPR)pivot;

  for( int i=0; i<fromat->m_actFW; i++ )
  {
    m_index[eqno][i+1] = fromat->m_frow[i].no;
    m_A[eqno][i+1]     = (REALPR)fromat->m_U[i];
  }
}


void CRSMAT::Addrow( int from, int to, REALPR factor )
{
  for( int i=0; i<m_width[from]; i++ )
  {
    for( int j=0; j<m_width[to]; j++ )
    {
      if( m_index[to][j] == m_index[from][j] )  m_A[to][j] += factor * m_A[from][i];
    }
  }
}

double CRSMAT::ScaleL2Norm( double* b, SUBDOM* subdom )
{
  double scale = 0.0;

//for( int i=0; i<m_neq_dn; i++ )  scale += m_A[i][0] * m_A[i][0];    // scaling type 1
  for( int i=0; i<m_neq_dn; i++ )  scale += b[i] * b[i];              // scaling type 2

# ifdef _MPI_
  scale = subdom->Mpi_sum( scale );
# endif

  scale = sqrt( scale );

  if( fabs(scale) > kZero )
  {
    for( int i=0; i<m_neq; i++ )
    {
      b[i] /= scale;

      for( int j=0; j<m_width[i]; j++)
      {
        m_A[i][j] /= (REALPR) scale;
      }
    }
  }

  return scale;
}

void CRSMAT::ScaleDiag( double* b )
{
  for( int i=0; i<m_neq; i++ )
  {
    double scale = m_A[i][0];

    if( fabs(scale) > kZero )
    {
      b[i] /= scale;

      for( int j=0; j<m_width[i]; j++)
      {
        m_A[i][j] /= (REALPR) scale;
      }
    }
  }
}

int CRSMAT::Getneq()
{
  return m_neq;
}


int CRSMAT::Getwidth( int r )
{
# ifdef kRangeCheck
  if( r < 0  ||  r >= m_neq )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Getwidth(1)" );
# endif

  return m_width[r];
}


int CRSMAT::Getindex( int r, int i )
{
# ifdef kRangeCheck
  if( r < 0  ||  r >= m_neq )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Getindex(1)" );
  if( i < 0  ||  i >= m_width[r] )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Getindex(2)" );
# endif

  return m_index[r][i];
}


void CRSMAT::Setwidth( int r, int w )
{
# ifdef kRangeCheck
  if( r < 0  ||  r >= m_neq )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Setwidth(1)" );
# endif

  m_width[r] = w;
}


void CRSMAT::Setindex( int r, int i, int c )
{
# ifdef kRangeCheck
  if( r < 0  ||  r >= m_neq )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Setindex(1)" );
  if( i < 0  ||  i >= m_width[r] )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Setindex(2)" );
# endif

  m_index[r][i] = c;
}


int CRSMAT::Find( int r, int c )
{
# ifdef kRangeCheck
  if( r < 0  ||  r >= m_neq )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Find(1)" );
# endif

  for( int i=0; i<m_width[r]; i++ )
  {
    if( m_index[r][i] == c )  return i;
  }

  return -1;
}


REALPR CRSMAT::Get( int r, int i )
{
# ifdef kRangeCheck
  if( r < 0  ||  r >= m_neq )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Get(1)" );
  if( i < 0  ||  i >= m_width[r] )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Get(2)" );
# endif

  return m_A[r][i];
}


void CRSMAT::Set( int r, int i, REALPR v )
{
# ifdef kRangeCheck
  if( r < 0  ||  r >= m_neq )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Set(1)" );
  if( i < 0  ||  i >= m_width[r] )
    REPORT::rpt.Error( kRangeFault, "range fault in CRSMAT::Set(2)" );
# endif

  m_A[r][i] = v;
}


//////////////////////////////////////////////////////////////////////////////////////////
// (matrix * vector) - multiplication:  r = A * x
//////////////////////////////////////////////////////////////////////////////////////////

double* CRSMAT::MulVec( double* x, double* r, PROJECT* project, EQS* eqs )
{
  for( int i=0; i<m_neq; i++ )
  {
    // multiplicate row "i" of "A" with "x"

    int*    IPtr = m_index[i];
    REALPR* APtr = m_A[i];
    int     w    = m_width[i];

    double  s    = APtr[0] * x[i];

    for( int j=1; j<w; j++ )
    {
      int k = IPtr[j];

      s += APtr[j] * x[k];
    }

    r[i] = s;
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // assemble local vector r from all adjacent subdomains
# ifdef _MPI_
  eqs->Mpi_assemble( r, project );
# endif
  ////////////////////////////////////////////////////////////////////////////////////////

  return r;
}
