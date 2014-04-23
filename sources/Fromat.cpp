// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class FROMAT
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
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
// /////////////////////////////////////////////////////////////////////////////////////////////////

#include "Defs.h"
#include "Report.h"
#include "Memory.h"
#include "Eqs.h"
#include "Node.h"
#include "Subdom.h"

#include "Fromat.h"

#define kDebug


FROMAT::FROMAT( int maxFW, int neq )
{
  m_neq   = neq;
  m_maxFW = maxFW;
  m_actFW = 0;

  m_frow = new FROW[m_maxFW];
  if( !m_frow )
    REPORT::rpt.Error( "can not allocate memory - FROMAT::FROMAT(1)" );

  int totsz = m_maxFW * m_maxFW;

  m_eqbuf = new double[totsz];
  if( !m_eqbuf )
    REPORT::rpt.Error( "can not allocate memory - FROMAT::FROMAT(2)" );

  for( int i=0; i<m_maxFW; i++ )
  {
    m_frow[i].eq = m_eqbuf + i*m_maxFW;
  }

  for( int i=0; i<totsz; i++ )  m_eqbuf[i] = 0.0;

  // -------------------------------------------------------------------------------------

  m_U = new double[m_maxFW];
  m_L = new double[m_maxFW];
  if( !m_U  ||  !m_L )
    REPORT::rpt.Error( "can not allocate memory - FROMAT::FROMAT(3)" );

  // -------------------------------------------------------------------------------------

  m_eql = new EQL[m_neq];
  if( !m_eql )
    REPORT::rpt.Error( "can not allocate memory - FROMAT::FROMAT(4)" );

  if( m_neq > 0 )
  {
    for( int i=0; i<m_neq; i++ )  m_eql[i].no = i;

    m_eql[0].prev       = NULL;
    m_eql[m_neq-1].next = NULL;

    for( int i=0; i<m_neq-1; i++ )  m_eql[i].next = &m_eql[i+1];
    for( int i=1; i<m_neq;   i++ )  m_eql[i].prev = &m_eql[i-1];
  }
}


FROMAT::~FROMAT()
{
  delete[] m_eqbuf;
  delete[] m_frow;
  delete[] m_U;
  delete[] m_L;
  delete[] m_eql;
}


//////////////////////////////////////////////////////////////////////////////////////////
// Insert equation "e" into the frontal equation matrix.
// 1. Increase the frontal matrix by the equation numbers connected to "e".
// 2. Fill in all new columns and rows.
//////////////////////////////////////////////////////////////////////////////////////////

void FROMAT::Insert( int e, int* width, int** index, REALPR** A )
{
  if( m_eql[e].cpl )  return;           // equation is already inserted


  // increase number of existing equations to hold equation "e" --------------------------

  int oldFW = m_actFW;

  for( int i=0; i<width[e]; i++ )
  {
    int r = index[e][i];

    if( m_eql[r].cpl )  continue;       // insert equations only once

    if( m_eql[r].ind < 0 )              // append new equation to front
    {
      if( m_actFW >= m_maxFW )
      {
        REPORT::rpt.Error( "maximum front width exceeded - FROMAT::Insert(1)" );
      //Resize( m_maxFW + 50 );
      }

      m_frow[m_actFW].no = r;
      m_eql[r].ind = m_actFW;
      m_actFW++;
    }
  }


  // insert new columns in existing equations --------------------------------------------

  for( int i=0; i<oldFW; i++ )
  {
    int     r  = m_frow[i].no;
    double* ep = m_frow[i].eq;

    for( int j=0; j<width[r]; j++ )
    {
      int k = m_eql[ index[r][j] ].ind;
      if( k >= oldFW )  ep[k] = A[r][j];
    }
  }


  // insert new equations ----------------------------------------------------------------

  for( int i=oldFW; i<m_actFW; i++ )
  {
    int     r  = m_frow[i].no;
    double* ep = m_frow[i].eq;

    for( int j=0; j<width[r]; j++ )
    {
      int k = m_eql[ index[r][j] ].ind;
      if( k >= 0 )  ep[k] = A[r][j];
    }
  }


  m_eql[e].cpl = true;
}


//////////////////////////////////////////////////////////////////////////////////////////
// Add equation "e" into the frontal equation matrix.
// Increase the frontal matrix by the equation numbers connected to "e".
//////////////////////////////////////////////////////////////////////////////////////////

void FROMAT::Add( int e, int n, int* eqno, double* eq )
{
  // increase number of existing equations to hold equation "e" --------------------------

  int oldFW = m_actFW;

  for( int i=0; i<n; i++ )
  {
    int r = eqno[i];

    if( m_eql[r].ind < 0 )              // append new equation to front
    {
      if( m_actFW >= m_maxFW )
        REPORT::rpt.Error( "maximum front width exceeded - FROMAT::Insert(2)" );
      //Resize( m_maxFW + 50 );

      m_frow[m_actFW].no = r;
      m_eql[r].ind = m_actFW;
      m_actFW++;
    }
  }


  // fill in equation "e" ----------------------------------------------------------------

  for( int i=0; i<m_actFW; i++ )
  {
    if( e == m_frow[i].no )
    {
      double* ep = m_frow[i].eq;

      for( int j=0; j<n; j++ )
      {
        int k = m_eql[ eqno[j] ].ind;
        if( k >= 0 )  ep[k] += eq[j];
      }
      break;
    }
  }


  m_eql[e].cpl = true;
}


//////////////////////////////////////////////////////////////////////////////////////////
// Gauss elimination of equation e.
//////////////////////////////////////////////////////////////////////////////////////////

void FROMAT::Eliminate( int e )
{
  int ind = m_eql[e].ind;               // front index to elimination equation

# ifdef kDebug
  if( ind == -1  ||  m_eql[e].cpl == false )
    REPORT::rpt.Error( "%s - FROMAT::Eliminate(1)",
                       "unexpected fault: trying to eliminate incomplete equation" );
# endif


  m_actFW--;                            // decrease the actual front width

  double* eeq   = m_frow[ind].eq;       // pointer to elimination equation
  double  pivot = eeq[ind];             // pivot


  // exchange elimination equation "e" with last equation "m_actFW"  ---------------------

  m_frow[ind].eq     = m_frow[m_actFW].eq;
  m_frow[m_actFW].eq = eeq;

  m_frow[ind].no     = m_frow[m_actFW].no;

  m_eql[m_frow[m_actFW].no].ind = ind;


  // eliminate the pivot -----------------------------------------------------------------

  eeq[ind] = eeq[m_actFW];
  m_eql[e].ind = -1;                    // remove equation "e"


  // Gauss elimination -------------------------------------------------------------------

  for( int i=0; i<m_actFW; i++ )        // loop on remaining equations
  {
    double* ieq  = m_frow[i].eq;        // pointer to equation i
    double  fac  = -ieq[ind] / pivot;   // elimination factor

    ieq[ind]     = ieq[m_actFW];        // eliminate column
    ieq[m_actFW] = 0.0;                 // and initialize

    for( int j=0; j<m_actFW; j++ )
    {
      ieq[j] += fac * eeq[j];
    }

    m_L[i] = fac;                       // save elimination factor
  }

  m_L[m_actFW] = 0.0;


  // save elimination equation "e" to m_U and initialize eeq[] ---------------------------

  for( int i=0; i<m_actFW; i++ )
  {
    m_U[i] = eeq[i];
    eeq[i] = 0.0;
  }

  m_U[m_actFW] = pivot;
  eeq[m_actFW] = 0.0;
}

//////////////////////////////////////////////////////////////////////////////////////////

void FROMAT::Eliminate( int e, double* B )
{
  int ind = m_eql[e].ind;               // front index to elimination equation

# ifdef kDebug
  if( ind == -1  ||  m_eql[e].cpl == false )
    REPORT::rpt.Error( "%s - FROMAT::Eliminate(1)",
                       "unexpected fault: trying to eliminate incomplete equation" );
# endif


  m_actFW--;                            // decrease the actual front width

  double* eeq   = m_frow[ind].eq;       // pointer to elimination equation
  double  pivot = eeq[ind];             // pivot


  // exchange elimination equation "e" with last equation "m_actFW"  ---------------------

  m_frow[ind].eq     = m_frow[m_actFW].eq;
  m_frow[m_actFW].eq = eeq;

  m_frow[ind].no     = m_frow[m_actFW].no;

  m_eql[m_frow[m_actFW].no].ind = ind;


  // eliminate the pivot -----------------------------------------------------------------

  eeq[ind] = eeq[m_actFW];
  m_eql[e].ind = -1;                    // remove equation "e"


  // Gauss elimination -------------------------------------------------------------------

  for( int i=0; i<m_actFW; i++ )        // loop on remaining equations
  {
    double* ieq  = m_frow[i].eq;        // pointer to equation i
    double  fac  = -ieq[ind] / pivot;   // elimination factor

    ieq[ind]     = ieq[m_actFW];        // eliminate column
    ieq[m_actFW] = 0.0;                 // and initialize

    for( int j=0; j<m_actFW; j++ )
    {
      ieq[j] += fac * eeq[j];
    }

    B[m_frow[i].no] += fac * B[e];      // RHS
  }


  // save elimination equation "e" to m_U and initialize eeq[] ---------------------------

  for( int i=0; i<m_actFW; i++ )
  {
    m_U[i] = eeq[i];
    eeq[i] = 0.0;
  }

  m_U[m_actFW] = pivot;
  eeq[m_actFW] = 0.0;
}


void FROMAT::MulVec( double* x, double* b, EQS* eqs, SUBDOM* subdom )
{
  for( int r=0; r<m_actFW; r++ )
  {
    b[r] = 0.0;

    for( int c=0; c<m_actFW; c++ )
    {
      b[r] += m_frow[r].eq[c] * x[c];
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // assemble local vector b from all adjacent subdomains

# ifdef _MPI_
  if( subdom->npr > 1 )
  {
    INFACE* inface = subdom->inface;

    int df = eqs->dfcn;


    // copy vector "b[]" to a temporary array --------------------------------------------

    double* tmp = (double*) MEMORY::memo.Array_eq( m_actFW );

    memcpy( tmp, b, m_actFW*sizeof(double) );


    // loop on all interfaces: exchange vector data --------------------------------------

    for( int s=0; s<subdom->npr; s++ )
    {
      MPI_Status status;

      int np = inface[s].np;

      if( np > 0 )
      {
        int cnt = 0;

        for( int n=0; n<np; n++ )
        {
          NODE* nd = inface[s].node[n];

          if( !isFS(nd->flag, NODE::kDry) )       // nothing to do if the node is dry...
          {
            SUB* sub = nd->sub;
            while( sub )
            {
              if( sub->no == s )  break;
              sub = sub->next;
            }

            if( !sub->dry )                       // ...or if the node is dry in
            {                                     // the adjacent subdomain s
              for( int e=0; e<df; e++ )
              {
                int eqno = eqs->GetEqno( nd, e );

                if( eqno >= 0 )
                {
#                 ifdef kDebug
                  if( m_eql[eqno].ind < 0 )
                    REPORT::rpt.Error( "%s - FROMAT::MulVec(1)",
                                       "unexpected fault: equation is not at front" );
#                 endif

                  inface[s].send[cnt] = tmp[m_eql[eqno].ind];
                  cnt++;
                }
              }
            }
          }
        }

        MPI_Sendrecv( inface[s].send, cnt, MPI_DOUBLE, s, 1,
                      inface[s].recv, cnt, MPI_DOUBLE, s, 1,
                      MPI_COMM_WORLD, &status );

        cnt = 0;

        for( int n=0; n<np; n++ )
        {
          NODE* nd = inface[s].node[n];

          if( !isFS(nd->flag, NODE::kDry) )
          {
            SUB* sub = nd->sub;
            while( sub )
            {
              if( sub->no == s )  break;
              sub = sub->next;
            }

            if( !sub->dry )
            {
              for( int e=0; e<df; e++ )
              {
                int eqno = eqs->GetEqno( nd, e );

                if( eqno >= 0 )
                {
                  b[m_eql[eqno].ind] += inface[s].recv[cnt];
                  cnt++;
                }
              }
            }
          }
        }
      }
    }

    MEMORY::memo.Detach( tmp );
  }
# endif
  ////////////////////////////////////////////////////////////////////////////////////////
}


void FROMAT::Output()
{
  char text[50];

  REPORT::rpt.Line2( 0 );
  REPORT::rpt.Output( "         " );

  for( int r=0; r<m_actFW; r++ )
  {
    sprintf( text, "   %8d  ", m_frow[r].no );
    REPORT::rpt.Output( text );
  }

  REPORT::rpt.Output( "\n" );

  for( int r=0; r<m_actFW; r++ )
  {
    sprintf( text, "%8d ", m_frow[r].no );
    REPORT::rpt.Output( text );

    for( int c=0; c<m_actFW; c++ )
    {
      sprintf( text, " %12.4le", m_frow[r].eq[c] );
      REPORT::rpt.Output( text );
    }

    REPORT::rpt.Output( "\n" );
  }
}
