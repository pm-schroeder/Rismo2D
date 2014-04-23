// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS
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

// ---------------------------------------------------------------------------------------
// Set up the index matrix for compact stored matrices.
//
// The index matrix "index" is stored in index as one dimensional array.
// The twodimensional array "eqInd" is set up to address each row in "index".
//
//
// Storage of the matrix:
//   eqInd[i][j]  :    column number of a non-zero entry of equation "i"
//                     where: 0 <= j <= totaltWidth[i]
//
//
// input parameters for function: int *setIndexMat()
//   fstelem      :    pointer to element list
//   np           :    total number of nodes
//   node         :    array of NODE structures
//   neq          :    total number of equations (rows of matrix)
//   df           :    degree of freedom at nodes
//   mceq         :    maximum number of connected equations expected
//
// output parameters
//   width        :    width array: normal length of equations
//
// return value
//   eqInd        :    two-dimensional array of index matrix
//                     the allocated memory can be freed with:
//                          delete[] totalWidth;    *** only if allocated ***
//                          delete[] eqInd[0];
//                          delete[] eqInd;
//                          delete[] width;
//
// Michael Schroeder in October 1994
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Elem.h"
#include "Node.h"
#include "Model.h"
#include "Project.h"
#include "CRSMat.h"

#include "Eqs.h"


void EQS::KillCrsm()
{
  if( crsm )
  {
    delete crsm;
    crsm = NULL;
  }

  initStructure = true;
}


void EQS::SetIndexMat( MODEL* model, int mceq )
{
  int   i, j, k;
  int   l, n, ceq;
  int   row, col;
  int*  indPtr;
  char  text[100];
  NODE* ndR;
  NODE* ndC;


  // allocate memory for index matrix ----------------------------------------------------

  crsm = new CRSMAT( neq, mceq );

  if( !crsm )
    REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - EQS::setIndexMat(1)" );

  crsm->m_neq_up = neq_up;
  crsm->m_neq_dn = neq_dn;

  int*  width = crsm->m_width;
  int** index = crsm->m_index;


  // initialization of index matrix and width arrays -------------------------------------

  j = mceq * neq;

  for( i=0; i<j;   i++ )  index[0][i] = -1;
  for( i=0; i<neq; i++ )  width[i]    =  0;


  // set up index matrix -----------------------------------------------------------------

  for( int e=0; e<model->ne; e++ )
  {
    ELEM* elem = model->elem[e];

    int nnd = elem->Getnnd();

    for( i=0; i<nnd; i++ )
    {
      ndR = elem->nd[i];

      for( j=0; j<dfcn; j++ )
      {
        row = GetEqno( ndR, j );

        if( row >= 0 )
        {
          indPtr = index[row];

          if( !width[row] )
          {
            width[row] = 1;
            indPtr[0]  = row;
          }

          for( k=0; k<nnd; k++ )
          {
            ndC = elem->nd[k];

            for( l=0; l<dfcn; l++ )
            {
              col = GetEqno( ndC, l );

              if( col >= 0 )
              {
                for( n=0; n<width[row]; n++ )
                {
                  if( indPtr[n] == col )
                  {
                    n = -1;
                    break;
                  }
                }

                if( n >= 0 )
                {
                  if( width[row] >= mceq )
                    REPORT::rpt.Error( "overflow in maximum connection: increase mceq!" );

                  indPtr[ width[row] ] = col;

                  width[row]++;
                }
              }
            }
          }


          for( l=0; l<dfel; l++ )
          {
            col = GetEqno( elem, l );

            if( col >= 0 )
            {
              for( n=0; n<width[row]; n++ )
              {
                if( indPtr[n] == col )
                {
                  n = -1;
                  break;
                }
              }

              if( n >= 0 )
              {
                if( width[row] >= mceq )
                  REPORT::rpt.Error( "overflow in maximum connection: increase mceq!" );

                indPtr[ width[row] ] = col;

                width[row]++;
              }
            }
          }
        }
      }
    }


    for( j=0; j<dfel; j++ )
    {
      row = GetEqno( elem, j );

      if( row >= 0 )
      {
        indPtr = index[row];

        if( !width[row] )
        {
          width[row] = 1;
          indPtr[0]  = row;
        }

        for( k=0; k<nnd; k++ )
        {
          ndC = elem->nd[k];

          for( l=0; l<dfcn; l++ )
          {
            col = GetEqno( ndC, l );

            if( col >= 0 )
            {
              for( n=0; n<width[row]; n++ )
              {
                if( indPtr[n] == col )
                {
                  n = -1;
                  break;
                }
              }

              if( n >= 0 )
              {
                if( width[row] >= mceq )
                  REPORT::rpt.Error( "overflow in maximum connection: increase mceq!" );

                indPtr[ width[row] ] = col;

                width[row]++;
              }
            }
          }
        }


        for( l=0; l<dfel; l++ )
        {
          col = GetEqno( elem, l );

          if( col >= 0 )
          {
            for( n=0; n<width[row]; n++ )
            {
              if( indPtr[n] == col )
              {
                n = -1;
                break;
              }
            }

            if( n >= 0 )
            {
              if( width[row] >= mceq )
                REPORT::rpt.Error( "overflow in maximum connection: increase mceq!" );

              indPtr[ width[row] ] = col;
              width[row]++;
            }
          }
        }
      }
    }
  }
/*
  if( solverType == kBicgstab_imp_3 )
  {
    int* owidth = new int[neq];
    if( !owidth )
      REPORT::rpt.Error( kMemoryFault, "cannot allocate memory - EQS::setIndexMat(2)" );

    for( i=0; i<neq; i++ )  owidth[i] = width[i];

    for( i=0; i<neq; i++ )
    {
      for( j=1; j<owidth[i]; j++ )
      {
        int eqj = index[i][j];

        if( eqj > i  &&  width[eqj] < mceq )
        {
          for( k=1; k<owidth[i]; k++ )
          {
            if( width[eqj] >= mceq ) break;

            int eqk = index[i][k];

            for( l=0; l<width[eqj]; l++ )
            {
              if( index[eqj][l] == eqk ) break;
            }

            if( l == width[eqj] )
            {
              index[eqj][ width[eqj] ]= eqk;
              width[eqj]++;
            }
          }
        }
      }
    }

    delete[] owidth;
  }
*/
  ceq = 0;

  for( i=0; i<neq; i++ )
  {
    if( width[i] > ceq ) ceq = width[i];
  }

  sprintf( text, "\n (EQS::SetIndexMat)      maximum connected equations %d\n", ceq );
  REPORT::rpt.Output( text, 2 );


//  SortIndex( neq, width, index );
}


void EQS::SortIndex( int neq, int* width, int** index )
{
  // maximum width of index matrix -------------------------------------------------------

  int maxWidth = width[0];

  for( int i=1; i<neq; i++ )  if( width[i] > maxWidth )  maxWidth = width[i];

  int* tmpIndex = new int[maxWidth];
  if( !tmpIndex ) REPORT::rpt.Error( "can not allocate memory - EQS::ILU_decomp (1)" );


  // sort elements of index matrix -------------------------------------------------------

  for( int i=0; i<neq; i++ )
  {
    tmpIndex[0] = index[i][1];

    int n = 1;

    for( int j=2; j<width[i]; j++ )
    {
      int k;
      for( k=0; k<n; k++ )
      {
        if( index[i][j] < tmpIndex[k] )  break;
      }

      for( int l=n-1; l>=k; l-- )  tmpIndex[l+1] = tmpIndex[l];

      tmpIndex[k] = index[i][j];
      n++;
    }

    for( int j=1; j<width[i]; j++ )  index[i][j] = tmpIndex[j-1];
  }

  delete[] tmpIndex;
}
