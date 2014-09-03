// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class MEMORY
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


MEMORY MEMORY::memo;


MEMORY::MEMORY( int array )
{
  m_max_nnd = 0;
  m_max_nel = 0;
  m_max_neq = 0;

  m_array = array;

  m_temp   = new ITEM* [m_array];
  m_size   = new unsigned int [m_array];
  m_flag   = new unsigned int [m_array];
  m_method = new char* [m_array];

  if( !m_temp || !m_size || !m_flag || !m_method )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory - MEMORY::MEMORY(1)" );

  for( int i=0; i<m_array; i++ )
  {
    m_temp[i]   = 0;
    m_size[i]   = 0;
    m_flag[i]   = 0;
    m_method[i] = new char[50];
  }
}


MEMORY::~MEMORY()
{
  for( int i=0; i<m_array; i++ ) Delete( i );
  for( int i=0; i<m_array; i++ ) delete[] m_method[i];

  delete[] m_temp;
  delete[] m_size;
  delete[] m_flag;
  delete[] m_method;
}


void* MEMORY::Array_nd( unsigned int nnd, char *method )
{
  if( nnd > m_max_nnd )
  {
    for( int i=0; i<m_array; i++ )
    {
      if( isFS(m_flag[i], kNd) )
      {
        if( !isFS(m_flag[i], kUsed) )
          Delete(i);
        else
          SF( m_flag[i], kDelete );
      }
    }

    m_max_nnd = nnd;
  }

  return Array( m_max_nnd, kNd, method );
}


void* MEMORY::Array_el( unsigned int nel, char *method )
{
  if( nel > m_max_nel )
  {
    for( int i=0; i<m_array; i++ )
    {
      if( isFS(m_flag[i], kEl) )
      {
        if( !isFS(m_flag[i], kUsed) )
          Delete(i);
        else
          SF( m_flag[i], kDelete );
      }
    }

    m_max_nel = nel;
  }

  return Array( m_max_nel, kEl, method );
}


void* MEMORY::Array_eq( unsigned int neq, char *method )
{
  if( neq > m_max_neq )
  {
    for( int i=0; i<m_array; i++ )
    {
      if( isFS(m_flag[i], kEq) )
      {
        if( !isFS(m_flag[i], kUsed) )
          Delete(i);
        else
          SF( m_flag[i], kDelete );
      }
    }

    m_max_neq = neq;
  }

  return Array( m_max_neq, kEq, method );
}


void* MEMORY::Array( unsigned int n, unsigned int flag, char *method )
{
  // -------------------------------------------------------------------------------------
  // look for an array that fits "n"

  for( int i=0; i<m_array; i++ )
  {
    if( !isFS(m_flag[i], kUsed) && isFS(m_flag[i], flag) )
    {
      //if( m_size[i] == n )
      if( m_size[i] >= n )
      {
        SF( m_flag[i], kUsed );
        if( method ) strcpy( m_method[i], method );
        return m_temp[i];
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // if no fitting array was found:
  // ... search for arrays that should be deleted

  for( int i=0; i<m_array; i++ )
  {
    if( !isFS(m_flag[i], kUsed)  &&  isFS(m_flag[i], kDelete) )
    {
      if( m_temp[i] ) delete[] m_temp[i];

      if( !(m_temp[i] = new ITEM[n]) )
        REPORT::rpt.Error( kMemoryFault, "can not allocate memory - MEMORY::Array(1)" );

      m_size[i] = n;

      SF( m_flag[i], kUsed );
      SF( m_flag[i], flag );
      if( method ) strcpy( m_method[i], method );
      return m_temp[i];
    }
  }

  // ... or search for unused arrays

  for( int i=0; i<m_array; i++ )
  {
    if( !m_temp[i] )
    {
      if( !(m_temp[i] = new ITEM[n]) )
        REPORT::rpt.Error( kMemoryFault, "can not allocate memory - MEMORY::Array(2)" );

      m_size[i] = n;

      SF( m_flag[i], kUsed );
      SF( m_flag[i], flag );
      if( method ) strcpy( m_method[i], method );
      return m_temp[i];
    }
  }

  // -------------------------------------------------------------------------------------

  REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Array(3)" );

  return NULL;
}


void MEMORY::Detach( void* temp )
{
  for( int i=0; i<m_array; i++ )
  {
    if( m_temp[i] == temp )
    {
      if( isFS(m_flag[i], kImatrix) )
      {
        Delete( i );
      }
      else if( isFS(m_flag[i], kDmatrix) )
      {
        Delete( i );
      }
      else
      {
        if( isFS(m_flag[i], kDelete) )  Delete( i );
        else                            CF( m_flag[i], kUsed );
      }
      return;
    }
  }

  REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Detach(1)" );
}


void MEMORY::Delete( void* temp )
{
  for( int i=0; i<m_array; i++ )
  {
    if( m_temp[i] == temp )  Delete( i );
    return;
  }

  REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Delete(1)" );
}


void MEMORY::Delete( int i )
{
  if( m_temp[i] )
  {
    if( isFS(m_flag[i], kImatrix) )
    {
      int** M = (int**) m_temp[i];
      delete[] M[0];
      delete[] M;
    }
    else if( isFS(m_flag[i], kDmatrix) )
    {
      double** M = (double**) m_temp[i];
      delete[] M[0];
      delete[] M;
    }
    else
    {
      delete[] m_temp[i];
    }

    m_temp[i] = NULL;
    m_size[i] = 0;
    m_flag[i] = 0;

    return;
  }

  REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Delete(1)" );
}


int** MEMORY::Imatrix( unsigned int rows, unsigned int cols )
{
  int** M = new int* [ rows ];
  if( !M )
    REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Imatrix(1)" );

  M[0] = new int[ rows * cols ];
  if( !M[0] )
    REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Imatrix(2)" );

  for( unsigned int i=1; i<rows; i++ )
  {
    M[i] = M[i-1] + cols;
  }

  for( int i=0; i<m_array; i++ )
  {
    if( !m_temp[i] )
    {
      m_temp[i] = (ITEM*) M;
      m_size[i] = rows * cols;

      SF( m_flag[i], kUsed );
      SF( m_flag[i], kImatrix );

      return M;
    }
  }


  REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Imatrix(3)" );

  return NULL;
}


double** MEMORY::Dmatrix( unsigned int rows, unsigned int cols )
{
  double** M = new double* [ rows ];
  if( !M )
    REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Dmatrix(1)" );

  M[0] = new double[ rows * cols ];
  if( !M[0] )
    REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Dmatrix(2)" );

  for( unsigned int i=1; i<rows; i++ )
  {
    M[i] = M[i-1] + cols;
  }

  for( int i=0; i<m_array; i++ )
  {
    if( !m_temp[i] )
    {
      m_temp[i] = (ITEM*) M;
      m_size[i] = rows * cols;

      SF( m_flag[i], kUsed );
      SF( m_flag[i], kDmatrix );

      return M;
    }
  }


  REPORT::rpt.Error( kUnexpectedFault, "unexpected internal fault - MEMORY::Dmatrix(3)" );

  return NULL;
}


void MEMORY::PrintInfo()
{
  unsigned int nused  = 0;
  unsigned int nalloc = 0;
  unsigned int nbytes = 0;

  for( int i=0; i<m_array; i++ )
  {
    if( m_temp[i] )  nalloc++;
    if( isFS(m_flag[i], kUsed) )  nused++;

         if( isFS(m_flag[i], kImatrix) )  nbytes += m_size[i] * sizeof(int);
    else if( isFS(m_flag[i], kDmatrix) )  nbytes += m_size[i] * sizeof(double);
    else                                  nbytes += m_size[i] * sizeof(ITEM);
  }


  char text[200];

  sprintf( text, "\n%-25s%u %s\n%-25s%u arrays | %u in use\n",
                 " (MEMORY)", nbytes,
                 "bytes of memory allocated",
                 "                       ", nalloc, nused );
  REPORT::rpt.Output( text, 3 );

  REPORT::rpt.Output( "\n", 5 );

  for( int i=0; i<m_array; i++ )
  {
    if( isFS(m_flag[i], kUsed)  &&  m_temp[i] )
    {
      sprintf( text, "%-25sarray[%3d] attached by method %s\n", " ", i, m_method[i] );
      REPORT::rpt.Output( text, 5 );
    }
    else if( m_temp[i] )
    {
      sprintf( text, "%-25sarray[%3d] detached by method %s\n", " ", i, m_method[i] );
      REPORT::rpt.Output( text, 5 );
    }
  }
}
