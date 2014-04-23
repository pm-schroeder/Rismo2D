// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class FRONT
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
#include "Elem.h"
#include "Node.h"
#include "Model.h"
#include "Project.h"

#include "Front.h"


//#define kDebug
//#define kElemCount
#define kMinPivot   1.0e-60


FRONT::FRONT()
{
  id = NULL;
  solverType = kFront;
}


FRONT::~FRONT()
{
}


//////////////////////////////////////////////////////////////////////////////////////////

void FRONT::Direct( PROJECT* project, MODEL* model, EQS* eqs, double* vec )
{
  int   i, j, k, l;

  char  text[240];

  int  record     = 0;
  int  actualFW   = 0;
  int  maximumFW  = 0;
  int  eqCounter  = 0;
  long bufCounter = 0;


# ifdef kDebug
  FILE *dbgId = fopen( "front.Report", "w" );
# endif


  int neq  = eqs->neq;
  int dfcn = eqs->dfcn;
  int dfel = eqs->dfel;

  double*  force  = eqs->force;
  double** estifm = eqs->estifm;

  long bufsz;

  if(  size <= 0  ||  !id )  bufsz = mfw * neq;
  else                       bufsz = mfw * size;


  // allocate temporary used memory ------------------------------------------------------

  int* frind = new int [eqs->maxEleq];

  if( !frind )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory - FRONT::Direct(1)" );


  int* frontWidth = new int [neq];

  if( !frontWidth )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory - FRONT::Direct(2)" );


  // if open, position scratch file to the beginning -------------------------------------

  if( id ) fseek( id, 0L, 0 );


  // allocate dynamic front memory -------------------------------------------------------

  REALPR* eqBuf   = new REALPR [bufsz];
  int*    eqnoBuf = new int    [bufsz];

  if( !eqBuf  ||  !eqnoBuf )
  {
    REPORT::rpt.Error( kMemoryFault, "%s - %lu bytes - FRONT::Direct(3)",
                           "can not allocate memory",
                           bufsz * (sizeof(int) + sizeof(REALPR)) );
  }


  int*     complete = new int     [mfw];
  int*     frntEqno = new int     [mfw];
  double** frntEq   = new double* [mfw];

  if( !complete  ||  !frntEqno  ||  !frntEq )
  {
    REPORT::rpt.Error( kMemoryFault, "%s - %lu bytes - FRONT::Direct(4)",
                           "can not allocate memory",
                           mfw * (2 * sizeof(int) + sizeof(double *)) );
  }


  double* frntEqBuf = new double [mfw * mfw];

  if( !frntEqBuf )
  {
    REPORT::rpt.Error( kMemoryFault, "%s - %lu bytes - FRONT::Direct(5)",
                           "can not allocate memory",
                           mfw * (mfw * sizeof(double)) );
  }


  for( i=0; i<mfw; i++ ) frntEq[i] = frntEqBuf + i*mfw;

  j = mfw * mfw;
  for( i=0; i<j; i++ ) frntEqBuf[i] = 0.0;


  // initialize the force vector ---------------------------------------------------------

  for( i=0; i<neq; i++ ) vec[i] = 0.0;


  // -------------------------------------------------------------------------------------
  // loop on all elements: create diagonal form of equation matrix

# ifdef kElemCount
  int  elemCounter = 0;
  REPORT::rpt.Screen(5, "\n" );
# endif


  for( int e=0; e<model->ne; e++ )
  {
    ELEM* elem = model->elem[e];

    long    cind, rind;
    long    indPiv;
    double* rowEq;
    double* estifmPtr;
    double  pivot;
    double  factor;
    double* rowI;
    double* rowPiv;


    int nnd = elem->Getnnd();


    // compute element coefficients ------------------------------------------------------

    if( !eqs->Coefs(elem, project, estifm, force) )  continue;


#   ifdef kElemCount
    if( !(elemCounter % 50) )
    {
      REPORT::rpt.Screen( 5, "\r                                                          " );
      REPORT::rpt.Screen(5,  "\r %6d ", elemCounter );
    }

    elemCounter++;
    REPORT::rpt.Screen( 5, "." );
#   endif


    // determine index to front for element equations ------------------------------------

    for( i=0; i<dfcn; i++ )
    {
      for( j=0; j<nnd; j++ )
      {
        rind = i*nnd + j;

        int eqno = eqs->GetEqno( elem->nd[j], i );

        if( eqno < 0 )                // no equation? ...
        {
          frind[rind] = -1;
        }

        else                          // ... search at front for equation number <eqno>
        {
          int found = -1;

          for( k=0; k<actualFW; k++ )
          {
            if( frntEqno[k] == eqno )
            {
              found = k;
              break;
            }
          }

          if( found >= 0 )            // equation already at front? ...
          {
            frind[rind] = found;
          }

          else                        // ... append new one
          {
            frind[rind] = actualFW;

            complete[actualFW] = false;
            frntEqno[actualFW] = eqno;

            actualFW++;

            if( actualFW >= mfw )
              REPORT::rpt.Error( kParameterFault, "maximum front width exceeded - FRONT::Direct(6)" );
          }


          // check for last occurence of equations

          complete[frind[rind]] = elem->isLast[i][j];
        }
      }
    }


    for( i=0; i<dfel; i++ )
    {
      rind = dfcn*nnd + i;

      int eqno = eqs->GetEqno( elem, i );

      if( eqno < 0 )                   // no equation?
      {
        frind[rind] = -1;
      }

      else
      {
        frind[rind] = actualFW;

        complete[actualFW] = true;
        frntEqno[actualFW] = eqno;

        actualFW++;

        if( actualFW >= mfw )
          REPORT::rpt.Error( kParameterFault, "maximum front width exceeded - FRONT::Direct(7)" );
      }
    }


    if( actualFW > maximumFW ) maximumFW = actualFW;


    // insert element stiffness matrix <estifm> and force vector <force> -----------------

    for( i=0; i<nnd; i++ )
    {
      for( j=0; j<dfcn; j++ )
      {
        rind = i + j*nnd;

        if( frind[rind] >= 0 )
        {
          int eqno = eqs->GetEqno( elem->nd[i], j );

          rowEq     = frntEq[frind[rind]];
          estifmPtr = estifm[rind];

          vec[eqno] += force[rind];

#         ifdef kDebug
          fprintf( dbgId, "\ninserting NODE equation %d ...\n", rowEqno[j] + 1 );
#         endif

          for( k=0; k<nnd; k++ )
          {
            for( l=0; l<dfcn; l++ )
            {
              cind = k + l*nnd;

              if( frind[cind] >= 0 )
              {
                rowEq[frind[cind]] += estifmPtr[cind];

#               ifdef kDebug
                fprintf( dbgId, "          col %d: %lf\n",
                                elem->nd[k]->eqno[l] + 1,
                                estifmPtr[cind] );
#               endif
              }
            }
          }

          for( l=0; l<dfel; l++ )
          {
            cind = dfcn*nnd + l;

            if( frind[cind] >= 0 )
            {
              rowEq[frind[cind]] += estifmPtr[cind];

#             ifdef kDebug
              fprintf( dbgId, "          col %d: %lf\n",
                              elem->eqno[l] + 1,
                              estifmPtr[cind] );
#             endif
            }
          }
        }
      }
    }


    for( j=0; j<dfel; j++ )
    {
      rind = dfcn*nnd + j;

      if( frind[rind] >= 0 )
      {
        int eqno = eqs->GetEqno( elem, j );

        rowEq     = frntEq[frind[rind]];
        estifmPtr = estifm[rind];

        vec[eqno] += force[rind];

#       ifdef kDebug
        fprintf( dbgId, "\ninserting ELEM equation %d ...\n", rowEqno[j] + 1 );
#       endif

        for( k=0; k<nnd; k++ )
        {
          for( l=0; l<dfcn; l++ )
          {
            cind = k + l*nnd;

            if( frind[cind] >= 0 )
            {
              rowEq[frind[cind]] += estifmPtr[cind];

#             ifdef kDebug
              fprintf( dbgId, "          col %d: %lf\n",
                              elem->nd[k]->eqno[l] + 1,
                              estifmPtr[cind] );
#             endif
            }
          }
        }

        for( l=0; l<dfel; l++ )
        {
          cind = dfcn*nnd + l;

          if( frind[cind] >= 0 )
          {
            rowEq[frind[cind]] += estifmPtr[cind];

#           ifdef kDebug
            fprintf( dbgId, "          col %d: %lf\n",
                            elem->eqno[l] + 1,
                            estifmPtr[cind] );
#           endif
          }
        }
      }
    }


    // -----------------------------------------------------------------------------------
    // eliminate complete equations

    for( ;; )
    {
      // search for PIVOT (diagonal), if any equation is complete ------------------------

      pivot  = 0.0;
      indPiv = -1;

      for( i=0; i<actualFW; i++ )
      {
        if( complete[i] )
        {
          double diag = fabs( frntEq[i][i] );

          if( diag > fabs(pivot) )
          {
            pivot  = frntEq[i][i];
            indPiv = i;
          }
        }
      }

      // any complete row found? (else continue with next element) -----------------------

      if( indPiv >= 0  &&  fabs(pivot) > kMinPivot )
      {
        // decrement the actual front width

        actualFW--;


        int elimEqno = frntEqno[indPiv];                  // elimination equation


        // exchange <pivot> row with <actualFW> row

        rowPiv           = frntEq[indPiv];
        frntEq[indPiv]   = frntEq[actualFW];
        frntEq[actualFW] = rowPiv;

        complete[indPiv] = complete[actualFW];
        frntEqno[indPiv] = frntEqno[actualFW];


        // eliminate the pivot

        rowPiv[indPiv] = rowPiv[actualFW];


        // Gauss elimination

        for( i=0; i<actualFW; i++ )                       // loop on remaining rows
        {
          rowI         = frntEq[i];                       // row <i>
          factor       = rowI[indPiv] / pivot;            // elimination factor
          rowI[indPiv] = rowI[actualFW];                  // eliminate column


          if( fabs(factor) > kMinPivot )
          {
            for( j=0; j<actualFW; j++ )
            {
              rowI[j] -= factor * rowPiv[j];
            }

            vec[frntEqno[i]] -= factor * vec[elimEqno];
          }

          rowI[actualFW] = 0.0;                           // initialize column
        }


        // save <pivot> equation to buffers

        frontWidth[eqCounter] = actualFW + 1;             // length of equation

        eqnoBuf[bufCounter] = elimEqno;                   // save the pivot
        eqBuf[bufCounter]   = 1.0;
        bufCounter++;

        vec[elimEqno] /= pivot;

        int*    eqnoBufPtr = &eqnoBuf[bufCounter];
        REALPR* eqBufPtr   = &eqBuf[bufCounter];


        for( i=0; i<actualFW; i++ )
        {
          eqnoBufPtr[i] = frntEqno[i];
          eqBufPtr[i]   = (REALPR) (rowPiv[i] / pivot);
        }


        // -------------------------------------------------------------------------------

#       ifdef kDebug
        {
          int i;
          fprintf( dbgId, "\nelimination equation is %d; pivot = %lf\n",
                          elimEqno+1, pivot );

          for( i=0; i<actualFW; i++ )
            fprintf( dbgId, "          col %6d: %lf\n",
                            frntEqno[i]+1, rowPiv[i] );
        }
#       endif


        bufCounter += actualFW;


        // initialize elimination equation -----------------------------------------------

        for( i=0; i<=actualFW; i++ )  rowPiv[i] = 0.0;


        // save buffers, if bufCounter exceeds bufsz - mfw -------------------------------

        if( bufCounter > bufsz - mfw )
        {
          WriteEq( bufCounter, bufsz, eqnoBuf, eqBuf );
          record++;
          bufCounter = 0;
        }

        eqCounter++;
      }

      else
      {
        break;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // check for singularity, other errors and give size message

  if( record > 0 )
  {
    sprintf(text, "\n%-25s%d %s %lu %s\n%-25s%d %s\n",
                  " (FRONT::Direct)",
                  record, "records of",
                  bufsz * (sizeof(int)+sizeof(REALPR)) + sizeof(long),
                  "bytes written",
                  " ", maximumFW,
                  "was maximum front width");
  }
  else
  {
    sprintf( text, "\n%-25s%d was maximum front width\n",
                   " (FRONT::Direct)", maximumFW );
  }
  REPORT::rpt.Output( text, 2 );


  if( actualFW > 0 )
  {
    REPORT::rpt.Error( kUnexpectedFault, "singular matrix - FRONT::Direct(8)" );
  }

  if( neq < eqCounter )
  {
    REPORT::rpt.Error( kUnexpectedFault, "bad number of equations - FRONT::Direct(9)" );
  }


  // -------------------------------------------------------------------------------------
  // backward substitution

  do
  {
    double* solve;
    int     eqnoPiv, eqno;

    if( bufCounter == 0 )
    {
      record--;
      ReadEq( record, bufsz, &bufCounter, eqnoBuf, eqBuf );
    }

    eqCounter--;

    actualFW    = frontWidth[eqCounter];
    bufCounter -= actualFW;


    // sum equation directly to <vec[eqnoPiv]>

    eqnoPiv = eqnoBuf[bufCounter];                       // index to pivot

    solve = &(vec[eqnoPiv]);                             // pointer to right side


    j = bufCounter + 1;


    for( i=1; i<actualFW; i++ )
    {
      eqno = eqnoBuf[j];

      *solve -= eqBuf[j] * vec[eqno];
      j++;
    }

    *solve /= eqBuf[bufCounter];                         // divide *solve by pivot

  } while( eqCounter );

  REPORT::rpt.Screen( 2, "\n" );


  // free temporary allocated memory -----------------------------------------------------

  delete[] frntEqBuf;
  delete[] frntEq;
  delete[] frntEqno;
  delete[] complete;
  delete[] eqnoBuf;
  delete[] eqBuf;
  delete[] frontWidth;
  delete[] frind;


# ifdef kDebug
  fclose( dbgId );
# endif
}


// ---------------------------------------------------------------------------------------
// read record containing bufCounter and contents of eqBuf from scratch
// unit first, the file pointer is set to the record number recno

void FRONT::ReadEq( int     recno,
                    size_t  bufsz,
                    long*   bufCounter,
                    int*    eqnoBuf,
                    REALPR* eqBuf )
{
  unsigned int items;
  long         position;


  // set file pointer position (count from beginning of file) ----------------------------

  position = recno * ( sizeof(long) + bufsz*sizeof(int)
                                    + bufsz*sizeof(REALPR) );
  if( fseek(id, position, 0) != 0)
    REPORT::rpt.Error( "can not position temporary file pointer - FRONT::ReadEq(1)" );


  // read record -------------------------------------------------------------------------

  items = 1;

  if( fread((char *)bufCounter, sizeof(long), items, id) != items )
    REPORT::rpt.Error( "reading temporary file - FRONT::ReadEq(2)" );


  items = (unsigned int) bufsz;

  if( fread((char *)eqnoBuf, sizeof(int), items, id) != items )
    REPORT::rpt.Error( "reading temporary file - FRONT::ReadEq(3)" );

  if( fread((char *)eqBuf, sizeof(REALPR), items, id) != items )
    REPORT::rpt.Error( "reading temporary file - FRONT::ReadEq(4)" );
}


// ---------------------------------------------------------------------------------------
// write sequential bufCounter and contents of eqBuf to scrach unit
// positioning of file pointer is not needed

void FRONT::WriteEq( long    bufCounter,
                     size_t  bufsz,
                     int*    eqnoBuf,
                     REALPR* eqBuf )
{
  unsigned int items;

  if(!id) REPORT::rpt.Error( "scratch unit not open - FRONT::WriteEq(1)" );


  items = (unsigned int) 1;

  if( fwrite ((char *)(&bufCounter), sizeof(long), items, id) != items )
    REPORT::rpt.Error( "writing temporary file - FRONT::WriteEq(2)" );


  items = (unsigned int) bufsz;

  if( fwrite((char *)eqnoBuf, sizeof(int), items, id) != items )
    REPORT::rpt.Error( "writing temporary file - FRONT::WriteEq(3)" );

  if( fwrite((char *)eqBuf, sizeof(REALPR), items, id) != items )
    REPORT::rpt.Error( "writing temporary file - FRONT::WriteEq(4)" );
}


//////////////////////////////////////////////////////////////////////////////////////////

void FRONT::OpenEquation()
{
  char theFile[200];

  if( size > 0 )
  {
    sprintf( theFile, "%stemp.scr", path );

    REPORT::rpt.Message( 2, "\n%-25s%s%s\n",
                            " (FRONT::OpenEquation)", "scratch file: ", theFile );


    id = fopen( theFile, "rb+");


    // try to create scratch file, if fopen has failed -----------------------------------

    if( !id )  id = fopen( theFile, "wb+");

    if( !id )  REPORT::rpt.Error( kOpenFileFault, "%s (EQS::openEquation - 1)",
                                  "can not open scratch files" );
  }
}


void FRONT::CloseEquation()
{
  if( id )  fclose( id );

  id = (FILE *) 0;
}
