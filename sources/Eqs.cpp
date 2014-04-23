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

#include "Defs.h"
#include "Report.h"
#include "CRSMat.h"
#include "Node.h"
#include "Elem.h"
#include "Type.h"
#include "Shape.h"
#include "Grid.h"
#include "Memory.h"
#include "Model.h"
#include "Subdom.h"
#include "Preco_ilu0.h"
#include "Solver.h"
#include "Front.h"
#include "Bicgstab.h"
#include "P_bcgstabd.h"
#include "P_fgmresd.h"
#include "Project.h"

#include "Eqs.h"


EQS::EQS( int dfcn, int dfmn, int dfel )
{
  this->dfcn = dfcn;
  this->dfmn = dfmn;
  this->dfel = dfel;


  // allocate memory for estifm and force ------------------------------------------------

  int df = (dfcn>=dfmn)? (dfcn):(dfmn);

  maxEleq = df * kMaxNodes2D  +  dfel;

  estifm = MEMORY::memo.Dmatrix( maxEleq, maxEleq );
  force  = new double [maxEleq];

  if( !force  || !estifm )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (EQS::EQS - 1)" );


  // further initializations -------------------------------------------------------------

  modelInit = 0;

  initStructure = true;

  iterCountCG = 0;

  nodeEqno = NULL;
  elemEqno = NULL;

  eqnoNode = NULL;

  crsm = NULL;
}


EQS::~EQS()
{
  if( nodeEqno )
  {
    delete[] nodeEqno[0];
    delete[] nodeEqno;
  }

  if( elemEqno )
  {
    delete[] elemEqno[0];
    delete[] elemEqno;
  }

  delete[] force;

  MEMORY::memo.Delete( estifm );
}


int EQS::Coefs( ELEM*    elem,
                PROJECT* project,
                double** estifm,
                double*  force )
{
  // -------------------------------------------------------------------------------------
  // initialization of arrays force[] and estifm[][]

  if( force )
  {
    for( int i=0; i<maxEleq; i++ )  force[i] = 0.0;
  }

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      double* estifmPtr = estifm[i];

      for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }

  return true;
}


int EQS::GetEqno( NODE* node, int no )
{
  return nodeEqno[no][node->Getno()];
}


int EQS::GetEqno( ELEM* elem, int no )
{
  if( isFS(elem->flag, ELEM::kBound) )  return -1;

  return elemEqno[no][elem->Getno()];
}


//////////////////////////////////////////////////////////////////////////////////////////

void EQS::ExportIM( const char* filename, CRSMAT* crsm, PROJECT* project )
{
  MODEL* model  = project->M2D;
  GRID*  region = model->region;

  FILE* id = fopen( filename, "w" );
  fprintf( id, "%d\n", neq );

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode( n );

    for( int e=0; e<dfcn; e++ )
    {
      int eqno = GetEqno( nd, e );

      if( eqno >= 0 )
      {
        fprintf( id, "%5d %1d (%5d) ", nd->Getname(), e, eqno );

        for( int i=0; i<crsm->Getwidth(eqno); i++ )
        {
          fprintf( id, " %5d", crsm->Getindex(eqno,i) );
        }

        fprintf( id, "\n" );
      }
    }
  }

  fclose( id );
}


void EQS::ExportEQS( const char* filename, CRSMAT* crsm, PROJECT* project )
{
  MODEL* model  = project->M2D;
  GRID*  region = model->region;

  FILE* id = fopen( filename, "w" );
  fprintf( id, "%d\n", neq );

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode( n );

    for( int e=0; e<dfcn; e++ )
    {
      int eqno = GetEqno( nd, e );

      if( eqno >= 0 )
      {
        fprintf( id, "%5d %1d (%5d) ", nd->Getname(), e, eqno );

        for( int i=0; i<crsm->Getwidth(eqno); i++ )
        {
          fprintf( id, " %12.4le", crsm->Get(eqno,i) );
        }

        fprintf( id, "\n" );
      }
    }
  }

  fclose( id );
}


void EQS::ExportVEC( const char* filename, double* vec, PROJECT* project )
{
  MODEL* model  = project->M2D;
  GRID*  region = model->region;

  FILE* id = fopen( filename, "w" );
  fprintf( id, "%d\n", neq );

   for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode( n );

    for( int e=0; e<dfcn; e++ )
    {
      int eqno = GetEqno( nd, e );

      fprintf( id, "%5d %1d (%5d) ", nd->Getname(), e, eqno );

      if( eqno >= 0 )  fprintf( id, "%14.8lf\n", vec[eqno] );
      else             fprintf( id, "%14.8lf\n", 0.0 );
    }
  }

  fclose( id );
}



void EQS::ScaleDiag( double* B, PROJECT* project )
{
  // set diagonal vector diag[] ----------------------------------------------------------
  double* diag = (double*) MEMORY::memo.Array_eq( neq );

  for( int i=0; i<neq; i++ )  diag[i] = crsm->m_A[i][0];

  Mpi_assemble( diag, project );


  // scale matrix ------------------------------------------------------------------------
  for( int i=0; i<neq; i++ )
  {
    double scale = diag[i];

    if( fabs(scale) > kZero )
    {
      B[i] /= scale;

      for( int j=0; j<crsm->m_width[i]; j++)
      {
        crsm->m_A[i][j] /= (REALPR) scale;
      }
    }
  }

  MEMORY::memo.Detach( diag );
}


//////////////////////////////////////////////////////////////////////////////////////////
// assemble vector "vec[]" across subdomains

void EQS::Mpi_assemble( double* vec, PROJECT* project )
{
# ifdef _MPI_
  if( project->subdom.npr > 1 )
  {
    SUBDOM* subdom = &project->subdom;
    INFACE* inface = subdom->inface;


    // copy vector "vec[]" to a temporary array ------------------------------------------

    double* tmp = (double*) MEMORY::memo.Array_eq( neq );

    memcpy( tmp, vec, neq*sizeof(double) );


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
              for( int e=0; e<dfcn; e++ )
              {
                int eqno = GetEqno( nd, e );

                if( eqno >= 0 )
                {
                  inface[s].send[cnt] = tmp[eqno];
                  cnt++;
                }
              }
            }
          }
        }

#       ifdef _MPI_DBG
        {
          int rcnt;
          MPI_Sendrecv( &cnt,  1, MPI_INT, s, 1,
                        &rcnt, 1, MPI_INT, s, 1,
                        MPI_COMM_WORLD, &status );

          if( cnt != rcnt )
          {
            REPORT::rpt.Warning( kUnexpectedFault, "%s (%d/%d) %s | sending from %d to %d",
                                 "Different number of values", cnt, rcnt,
                                 "- EQS::Mpi_assemble(1)",  subdom->pid+1, s+1 );

            REPORT::rpt.Output( "### Following list of equations for MPI_Sendrecv...\n" );

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
                  for( int e=0; e<dfcn; e++ )
                  {
                    int eqno = GetEqno( nd, e );

                    if( eqno >= 0 )
                    {
                      int bckind = 0;
                      if( nd->bc )  bckind = nd->bc->kind;

                      char text[200];
                      sprintf( text, "     NODE %7d | eq %d: flag=%d; bc=%d; dry(%d)=%d",
                                     nd->Getname(), e, nd->flag, bckind,
                                     subdom->pid+1, (nd->flag&NODE::kDry)/NODE::kDry );
                      REPORT::rpt.Output( text );

                      SUB* sub = nd->sub;
                      while( sub )
                      {
                        sprintf( text, "; dry(%d)=%d", sub->no+1, sub->dry );
                        REPORT::rpt.Output( text );
                        sub = sub->next;
                      }

                      REPORT::rpt.Output( "\n" );
                    }
                  }
                }
              }
            }

            project->M2D->Output( project, 0 );
            M2D->DetachOutput();
            REPORT::rpt.Error( kMPI_Error, "Error in Mpi_assemble()" );
          }
        }
#       endif

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
              for( int e=0; e<dfcn; e++ )
              {
                int eqno = GetEqno( nd, e );

                if( eqno >= 0 )
                {
                  vec[eqno] += inface[s].recv[cnt];
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
}


void EQS::DataOut( char* name, int step, char* time, int release, GRID* region, char* label[], ... )
{
  char fileName[80];
  sprintf( fileName, name, step );

  // open file ---------------------------------------------------------------------------
  FILE* id = fopen( fileName, "w" );

  if( !id )
  {
    REPORT::rpt.Error( kOpenFileFault,
                       "can not open data output file (EQS::DataOut - 1)" );
    return;
  }

  // count number of data ----------------------------------------------------------------
  int ndata = 0;
  while( label[ndata][0] && label[ndata][0] != '\n')  ndata++;

  int edata = 0;
  while( label[edata+ndata+1][0] )  edata++;

  // write file header -------------------------------------------------------------------
  fprintf( id, "%d  %d  %d  %d  %d\n", region->Getnp(), region->Getne(), ndata, edata, 0 );

  // write nodes -------------------------------------------------------------------------
  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);
    fprintf( id, "%7d  %12.6lf  %12.6lf  %12.6lf\n",
             nd->Getname(), nd->x, nd->y, nd->zor );
  }

  // write element connectivity --------------------------------------------------------
  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* el = region->Getelem(e);

    char elemShape[6];

    switch( el->GetshapeID() )
    {
      case kLine:      strcpy( elemShape, "line" );   break;
      case kTriangle:  strcpy( elemShape, "tri" );    break;
      case kSquare:    strcpy( elemShape, "quad" );   break;
    }

    int mat = TYPE::Getid(el->type)->no(el->type);
    if( mat <= 0 )  mat = abs(el->type);

    fprintf( id, "%7d  %5d  %5s ", el->Getname(), mat, elemShape );

    int nnd = el->Getnnd();
    for( int i=0; i<nnd; i++ )  fprintf( id, " %7d", el->nd[i]->Getname() );

    fprintf( id, "\n" );
  }

  // write labels and units --------------------------------------------------------------
  fprintf( id, "%d", ndata );
  for( int i=0; i<ndata; i++ )  fprintf( id, " 1" );
  fprintf( id, "\n" );

  for( int i=0; i<ndata; i++ )  fprintf( id, "%s,-\n", label[i] );

  // write nodal values ------------------------------------------------------------------
  va_list argPtr;

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);
    int   no = nd->Getno();

    fprintf( id, "%7d", nd->Getname() );

    va_start( argPtr, label );

    for( int i=0; i<ndata; i++ )
    {
      double* val = va_arg( argPtr, double* );
      fprintf( id, " %14.6le", val[no] );
    }

    fprintf( id, "\n" );

    va_end( argPtr );
  }

  // write labels and units --------------------------------------------------------------
  fprintf( id, "%d", edata );
  for( int i=0; i<edata; i++ )  fprintf( id, " 1" );
  fprintf( id, "\n" );

  for( int i=0; i<edata; i++ )  fprintf( id, "%s,-\n", label[i+ndata+1] );

  // write element values ----------------------------------------------------------------
  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* el = region->Getelem(e);
    int   no = el->Getno();

    fprintf( id, "%7d", el->Getname() );

    va_start( argPtr, label );

    for( int i=0; i<ndata; i++ )  va_arg( argPtr, double* );

    for( int i=0; i<edata; i++ )
    {
      double* val = va_arg( argPtr, double* );
      fprintf( id, " %14.6le", val[no] );
    }
    fprintf( id, "\n" );

    va_end( argPtr );
  }

  // close file --------------------------------------------------------------------------
  fclose( id );
}
