// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class PRECO_ILUT
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
#include "Node.h"
#include "Project.h"
#include "Model.h"
#include "CRSMat.h"
#include "Eqs.h"
#include "Fromat.h"

#include "Preco_ilut.h"

#define kMinPivot   1.0e-20
#define kZero       1.0e-40

#define kEqCounter

//#define kDebug


PRECO_ILUT::PRECO_ILUT()
{
}

PRECO_ILUT::~PRECO_ILUT()
{
}


//////////////////////////////////////////////////////////////////////////////////////////
// Incomplete LU factorization (ILU)
// Matrix A is assumed to be symmetric in structure (not in values)
//////////////////////////////////////////////////////////////////////////////////////////

void PRECO_ILUT::Factor( PROJECT* project, EQS* eqs, CRSMAT* crsm )
{
  REPORT::rpt.Message( 2, "\n (PRECO_ILUT::Factor)    preconditioning: ILUT\n" );


  ////////////////////////////////////////////////////////////////////////////////////////
  // 1.  initializations and memory allocation

  int*     width = crsm->m_width;
  int**    index = crsm->m_index;
  REALPR** A     = crsm->m_A;

  //int*     Iw    = crsi->m_width;
  //int**    Ii    = crsi->m_index;
  REALPR** ILU   = crsi->m_A;

  int ceq    = crsi->m_ceq;

  int neq    = crsm->m_neq;
  int neq_dn = crsm->m_neq_dn;
  int neq_up = crsm->m_neq_up;

  //for( int e=0; e<neq; e++ )  Iw[e] = 1;

  //crsi->Init();


  // allocate dynamic front memory -------------------------------------------------------

  int mfw = 4 * ceq;

  fromat = new FROMAT( mfw, neq );

  if( !fromat )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory - PRECO_ILUT::Factor(1)" );


  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.  LU factorization
  //
  // 2.1 MPI: incomplete LU factorization restricted to interior subdomain nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Factor)    starting with 2.1\n" );
# endif

  int maximumFW = 0;

  for( int e=0; e<neq_up; e++ )
  {
    if( fabs(A[e][0]) < kMinPivot )
    {
      REPORT::rpt.Error( kParameterFault, "%s - PRECO_ILUT::Factor(2)",
                         "factorization can not be applied to the matrix" );
    }

    // equation "e" is elimination equation ----------------------------------------------

    fromat->Insert( e, width, index, A );

    if( fromat->m_actFW > maximumFW )  maximumFW = fromat->m_actFW;

#   ifdef kEqCounter
    if( ((e+1)%1000) == 0 )
    {
      REPORT::rpt.Screen( 2, "\r (PRECO_ILUT::Factor)    %d (maximum FW = %d)", e+1, maximumFW );
    }
#   endif

    fromat->Eliminate( e );

    ILU[e][0] = (REALPR) fromat->m_U[fromat->m_actFW];      // diagonal element

    //if( fabs(ILU[e][0]) < kZero )
    //{
    //  REPORT::rpt.Error( kParameterFault, "%s - PRECO_ILUT::Factor(2)",
    //                     "factorization can not be applied to the matrix" );
    //}

    for( int i=1; i<width[e]; i++ )
    {
      int f = index[e][i];
      int j = fromat->m_eql[f].ind;

      if( f > e )
      {
        ILU[e][i] = (REALPR) fromat->m_U[j];

        for( int k=1; k<width[f]; k++ )
        {
          if( index[f][k] == e )
          {
            ILU[f][k] = (REALPR) fromat->m_L[j];
            break;
          }
        }
      }
    }

    //fromat->Sort();

    //Ii[e][0]  = e;
    //ILU[e][0] = (REALPR) fromat->m_eqtr[fromat->m_actFW].U;

    //EQTREE* eqtr = (EQTREE*) fromat->m_root.First();

    //while( eqtr )
    //{
    //  int i;
    //  int f = eqtr->e;

    //  i         = Iw[e];
    //  Ii[e][i]  = f;
    //  ILU[e][i] = (REALPR) eqtr->U;
    //  Iw[e]++;

    //  i         = Iw[f];
    //  Ii[f][i]  = e;
    //  ILU[f][i] = (REALPR) eqtr->L;
    //  Iw[f]++;

    //  if( Iw[e] >= ceq  ||  Iw[f] >= ceq/2 ) break;

    //  eqtr = (EQTREE*) eqtr->Next();
    //}


    // ...
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.2 MPI communication:
  //     receive interface equations from subdomains with smaller pid (2 <- 1)
  //     s < pid  ==> upstream interface nodes

//# ifdef _MPI_DBG
//  REPORT::rpt.Output( " (PRECO_ILUT::Factor)    starting with 2.2\n" );
//# endif

# ifdef _MPI_

  if( project->subdom.npr > 1 )
  {
    SUBDOM* subdom = &project->subdom;
    INFACE* inface = subdom->inface;

    MPI_Status status;

    int df = eqs->dfcn;

    for( int s=subdom->npr-1; s>=0; s-- )
    {
      int np = inface[s].np;                      // number of interface nodes

      if( np > 0  &&  s < subdom->pid )           // upstream interface nodes
      {
        for( int n=0; n<np; n++ )
        {
          NODE* ndr = inface[s].node[n];

          if( !isFS(ndr->flag, NODE::kDry) )      // nothing to do if node is dry...
          {
            SUB* subr = ndr->sub;
            while( subr )
            {
              if( subr->no == s )  break;
              subr = subr->next;
            }

            if( !subr->dry )                       // ...or if the node is dry in
            {                                     // any adjacent subdomain
              for( int e=0; e<df; e++ )
              {
                int req = eqs->GetEqno( ndr, e );

                if( req >= 0 )
                {
                  int cnt;

                  MPI_Recv( &cnt, 1, MPI_INT, s, 1, MPI_COMM_WORLD, &status );

                  if( cnt )
                  {
                    MPI_Recv( inface[s].ria1, cnt, MPI_CHAR,   s, 2, MPI_COMM_WORLD, &status );
                    MPI_Recv( inface[s].ria2, cnt, MPI_INT,    s, 3, MPI_COMM_WORLD, &status );
                    MPI_Recv( inface[s].recv, cnt, MPI_DOUBLE, s, 4, MPI_COMM_WORLD, &status );

                    for( int j=0; j<cnt; j++ )
                    {
                      int   eid  = inface[s].ria1[j];
                      int   name = inface[s].ria2[j];
                      NODE* ndc  = subdom->node[name - 1];

                      if( ndc )
                      {
                        int ceq = eqs->GetEqno( ndc, eid );

                        for( int c=0; c<width[req]; c++ )
                        {
                          if( index[req][c] == ceq )
                          {
                            ILU[req][c] += (REALPR)inface[s].recv[j];
                            break;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.3 MPI: factorization of interface nodes

//# ifdef _MPI_DBG
//  REPORT::rpt.Output( " (PRECO_ILUT::Factor)    starting with 2.3\n" );
//# endif

  for( int e=neq_up; e<neq_dn; e++ )
  {
    if( fabs(ILU[e][0]) < kMinPivot )
    {
      int    piv = -1;
      double max = 0.0;

      for( int j=1; j<width[e]; j++ )
      {
        int eq = index[e][j];

        if( eq > e )
        {
          // find column "e" in equation "eq" --------------------------------------------

          for( int k=0; k<width[eq]; k++ )
          {
            if( index[eq][k] == e )
            {
              double p = fabs( ILU[eq][k] );

              if( p > max )
              {
                piv = eq;
                max = p;
              }
            }
          }
        }
      }

      if( piv > 0 )
      {
        for( int j=0; j<width[e]; j++ )
        {
          for( int k=0; k<width[piv]; k++ )
          {
            if( index[e][j] == index[piv][k] )  ILU[e][j] += ILU[piv][k];
          }
        }
      }
      else
      {
        REPORT::rpt.Error( kParameterFault, "ILU-Factorization cannot be applied to the matrix" );
      }
    }


    // equation "e" is elimination equation ----------------------------------------------

    for( int j=1; j<width[e]; j++ )
    {
      int eq = index[e][j];

      if( eq > e )// &&  eq > neq_dn ) // Note: The elimination of just received upstream
      {                                // equations (2.2) was performed in prior subdomain
        double factor = 0.0;

        // find column "e" in equation "eq" and compute elimination factor ---------------

        int k;

        for( k=0; k<width[eq]; k++ )
        {
          if( index[eq][k] == e )
          {
            factor = -ILU[eq][k] / ILU[e][0];
            break;
          }
        }


        // save elimination factor (L-Matrix: index[eq][k] == e < eq) --------------------

        ILU[eq][k] = (REALPR)factor;


        // add elimination equation ILU[i] to ILU[eq] ------------------------------------

        for( k=1; k<width[e]; k++ )
        {
          if( index[e][k] >= e )
          {
            for( int l=0; l<width[eq]; l++ )
            {
              if( index[eq][l] == index[e][k] )
              {
                ILU[eq][l] += (REALPR)factor * ILU[e][k];
                break;
              }
            }
          }
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.4 MPI communication:
  //     send to subdomain with larger pid      (2 -> 3)
  //     s > pid  ==> downstream interface nodes

//# ifdef _MPI_DBG
//  REPORT::rpt.Output( " (PRECO_ILUT::Factor)    starting with 2.4\n" );
//# endif

  if( project->subdom.npr > 1 )
  {
    SUBDOM* subdom = &project->subdom;
    INFACE* inface = subdom->inface;

    int df = eqs->dfcn;

    for( int s=0; s<subdom->npr; s++ )
    {
      int np = inface[s].np;                      // number of interface nodes

      if( np > 0  &&  s > subdom->pid )           // downstream interface nodes
      {
        for( int n=0; n<np; n++ )
        {
          NODE* ndr = inface[s].node[n];

          if( !isFS(ndr->flag, NODE::kDry) )      // nothing to send if the node is dry...
          {
            SUB* subr = ndr->sub;
            while( subr )
            {
              if( subr->no == s )  break;
              subr = subr->next;
            }

            if( !subr->dry )                      // ...or if the node is dry in
            {                                     // any adjacent subdomain
              for( int e=0; e<df; e++ )
              {
                int req = eqs->GetEqno( ndr, e );

                if( req >= 0 )
                {
                  int cnt = 0;

                  for( int c=0; c<width[req]; c++ )
                  {
                    int   ceq = index[req][c];
                    NODE* ndc = eqs->eqnoNode[ceq];
                    int   eid = eqs->eqid[ceq];

                    // is ndc a point on interface "s"?
                    SUB* subc = ndc->sub;
                    while( subc )
                    {
                      if( subc->no == s )  break;
                      subc = subc->next;
                    }

                    if( subc )
                    {
                      inface[s].sia1[cnt] = eid;
                      inface[s].sia2[cnt] = ndc->Getname();
                      inface[s].send[cnt] = ILU[req][c];
                      cnt++;
                    }
                  }

                  MPI_Send( &cnt, 1, MPI_INT, s, 1, MPI_COMM_WORLD );

                  if( cnt )
                  {
                    MPI_Send( inface[s].sia1, cnt, MPI_CHAR,   s, 2, MPI_COMM_WORLD );
                    MPI_Send( inface[s].sia2, cnt, MPI_INT,    s, 3, MPI_COMM_WORLD );
                    MPI_Send( inface[s].send, cnt, MPI_DOUBLE, s, 4, MPI_COMM_WORLD );
                  }
                }
              }
            }
          }
        }
      }
    }
  }
# endif
  ////////////////////////////////////////////////////////////////////////////////////////

# ifdef kDebug
  {
    char filename[80];

    if( project->subdom.npr == 1 )
      sprintf( filename, "ilu_2.4.dbg" );
    else
      sprintf( filename, "ilu_2.4_%02d.dbg", project->subdom.pid+1 );

    eqs->ExportEQS( filename, crsi, project );
  }

# endif

//# ifdef _MPI_DBG
//  REPORT::rpt.Output( " (PRECO_ILUT::Factor)    ILU factorization finished\n" );
//# endif
}


//////////////////////////////////////////////////////////////////////////////////////////
// forward and backward solve; determine X from: U * X  =  L * B
//////////////////////////////////////////////////////////////////////////////////////////

void PRECO_ILUT::Solve( PROJECT* project, EQS* eqs, double* B, double* X )
{
  int      neq    = crsi->m_neq;
  int      neq_dn = crsi->m_neq_dn;
  int      neq_up = crsi->m_neq_up;

  int*     width  = crsi->m_width;
  int**    index  = crsi->m_index;
  REALPR** ILU    = crsi->m_A;


  ////////////////////////////////////////////////////////////////////////////////////////
  // 1.  forward solve: manipulate vector B with lower part of matrix ILU
  //
  // 1.1 MPI: forward solve for interior subdomain nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 1.1\n" );
# endif

  for( int i=0; i<neq_up; i++ )
  {
    X[i] = B[i];

    for( int j=1; j<width[i]; j++ )
    {
      int eq = index[i][j];

      if( eq < i )  // L-matrix
      {
        X[i] += ILU[i][j] * X[eq];
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 1.2 MPI communication:
  //     receive from subdomains with smaller pid  (2 <- 1)
  //     s < pid  ==> upstream interface nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 1.2\n" );
# endif

# ifdef _MPI_

  if( project->subdom.npr > 1 )
  {
    SUBDOM* subdom = &project->subdom;
    INFACE* inface = subdom->inface;

    MPI_Status status;

    int df = eqs->dfcn;

    for( int s=subdom->npr-1; s>=0; s-- )
    {
      int np = inface[s].np;                      // number of interface nodes

      if( np > 0  &&  s < subdom->pid )           // upstream interface nodes
      {
        int cnt = 0;

#       ifdef _MPI_DBG
        {
          char text[200];
          sprintf( text, " ### receiving from %d:", s+1 );
          REPORT::rpt.Output( text );
        }
#       endif

        MPI_Recv( &cnt, 1, MPI_INT, s, 1, MPI_COMM_WORLD, &status );

#       ifdef _MPI_DBG
        {
          char text[200];
          sprintf( text, " %d values\n", cnt );
          REPORT::rpt.Output( text );
        }
#       endif

        if( cnt )
        {
          MPI_Recv( inface[s].recv, cnt, MPI_DOUBLE, s, 2, MPI_COMM_WORLD, &status );

          cnt = 0;

          for( int n=0; n<np; n++ )
          {
            NODE* nd = inface[s].node[n];

            if( !isFS(nd->flag, NODE::kDry) )       // nothing received if node is dry...
            {
              SUB* sub = nd->sub;
              while( sub )
              {
                if( sub->no == s )  break;
                sub = sub->next;
              }

              if( !sub->dry )                       // ...or if the node is dry in
              {                                     // any adjacent subdomain
                for( int e=0; e<df; e++ )
                {
                  int eqno = eqs->GetEqno( nd, e );

                  if( eqno >= 0 )
                  {
                    X[eqno] = inface[s].recv[cnt];
                    cnt++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 1.3 MPI: forward solve for upstream interface nodes
  //          Note: X[i] is not initialized with B[i] since
  //                this was performed in prior subdomain

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 1.3\n" );
# endif

  for( int i=neq_up; i<neq_dn; i++ )
  {
    for( int j=1; j<width[i]; j++ )
    {
      int eq = index[i][j];

      if( eq < i )
      {
        X[i] += ILU[i][j] * X[eq];
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 1.4 MPI: faktorisation of downstream interface nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 1.4\n" );
# endif

  for( int i=neq_dn; i<neq; i++ )
  {
    X[i] = B[i];

    for( int j=1; j<width[i]; j++ )
    {
      int eq = index[i][j];

      if( eq < neq_dn )
      {
        X[i] += ILU[i][j] * X[eq];
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 1.5 MPI communication:
  //     send to subdomains with larger pid        (2 -> 3)
  //     s > pid  ==> downstream interface nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 1.5\n" );
# endif

  if( project->subdom.npr > 1 )
  {
    SUBDOM* subdom = &project->subdom;
    INFACE* inface = subdom->inface;

    int df = eqs->dfcn;

    for( int s=0; s<subdom->npr; s++ )
    {
      int np = inface[s].np;                      // number of interface nodes

      if( np > 0  &&  s > subdom->pid )           // downstream interface nodes
      {
        int cnt = 0;

        for( int n=0; n<np; n++ )
        {
          NODE* nd = inface[s].node[n];

          if( !isFS(nd->flag, NODE::kDry) )       // nothing to send if the node is dry...
          {
            SUB* sub = nd->sub;
            while( sub )
            {
              if( sub->no == s )  break;
              sub = sub->next;
            }

            if( !sub->dry )                       // ...or if the node is dry in
            {                                     // any adjacent subdomain
              for( int e=0; e<df; e++ )
              {
                int eqno = eqs->GetEqno( nd, e );

                if( eqno >= 0 )
                {
                  inface[s].send[cnt] = X[eqno];
                  cnt++;
                }
              }
            }
          }
        }

#       ifdef _MPI_DBG
        {
          char text[200];
          sprintf( text, " ### sending %d values to %d\n", cnt, s+1 );
          REPORT::rpt.Output( text );
        }
#       endif

        MPI_Send( &cnt, 1, MPI_INT, s, 1, MPI_COMM_WORLD );
        if( cnt ) MPI_Send( inface[s].send, cnt, MPI_DOUBLE, s, 2, MPI_COMM_WORLD );
      }
    }
  }

# endif
  ////////////////////////////////////////////////////////////////////////////////////////

# ifdef kDebug
  {
    FILE* dbg;
    char  filename[80];

    MODEL* model  = project->M2D;
    GRID*  region = model->region;

    if( project->subdom.npr == 1 )
      sprintf( filename, "forwX.dbg" );
    else
      sprintf( filename, "forwX_%02d.dbg", project->subdom.pid );

    dbg = fopen( filename, "w" );
    fprintf( dbg, "%d\n", region->Getnp() );

    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE* nd = region->Getnode( n );
      for( int e=0; e<eqs->dfcn; e++ )
      {
        int eqno = eqs->GetEqno( nd, e );
        fprintf( dbg, "%5d %1d ", nd->Getname(), e );
        if( eqno >= 0 )
          fprintf( dbg, "%14.6le\n", X[eqno] );
        else
          fprintf( dbg, "%14.6le\n", 0.0 );
      }
    }

    fclose( dbg );
  }
# endif

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     forward solve finished\n" );
# endif

  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.  solve for X with upper part of matrix ILU (backward substitution)
  //
  // 2.1 MPI communication:
  //     receive from subdomains with larger pid    (2 <- 3)
  //     s > pid  ==> downstream interface nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 2.1\n" );
# endif

# ifdef _MPI_

  if( project->subdom.npr > 1 )
  {
    SUBDOM* subdom = &project->subdom;
    INFACE* inface = subdom->inface;

    MPI_Status status;

    int df = eqs->dfcn;

    for( int s=subdom->npr-1; s>=0; s-- )
    {
      int np = inface[s].np;                // number of interface nodes

      if( np > 0  &&  s > subdom->pid )     // downstream interface nodes
      {
        int cnt = 0;

#       ifdef _MPI_DBG
        {
          char text[200];
          sprintf( text, " ### receiving from %d:", s+1 );
          REPORT::rpt.Output( text );
        }
#       endif

        MPI_Recv( &cnt, 1, MPI_INT, s, 1, MPI_COMM_WORLD, &status );

#       ifdef _MPI_DBG
        {
          char text[200];
          sprintf( text, " %d values\n", cnt );
          REPORT::rpt.Output( text );
        }
#       endif

        if( cnt )
        {
          MPI_Recv( inface[s].recv, cnt, MPI_DOUBLE, s, 2, MPI_COMM_WORLD, &status );

          cnt = 0;

          for( int n=0; n<np; n++ )
          {
            NODE* nd = inface[s].node[n];

            if( !isFS(nd->flag, NODE::kDry) )       // nothing received if node is dry...
            {
              SUB* sub = nd->sub;
              while( sub )
              {
                if( sub->no == s )  break;
                sub = sub->next;
              }

              if( !sub->dry )                       // ...or if the node is dry in
              {                                     // any adjacent subdomain
                for( int e=0; e<df; e++ )
                {
                  int req = eqs->GetEqno( nd, e );

                  if( req >= 0 )
                  {
                    X[req] = inface[s].recv[cnt];
                    cnt++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.2 MPI: solve for upstream interface nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 2.2\n" );
# endif

  for( int i=neq_dn-1; i>=neq_up; i-- )
  {
    for( int j=1; j<width[i]; j++ )
    {
      int eq = index[i][j];

      if( eq > i )
      {
        X[i] -= ILU[i][j] * X[eq];
      }
    }

    if( fabs(ILU[i][0]) < kZero )
      REPORT::rpt.Error( kParameterFault, "division by zero - EQS::ILU_solver(1)" );

    X[i] /= ILU[i][0];
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.3 MPI communication:
  //     send to subdomains with smaller pid        (2 -> 1)
  //     s < pid  ==> upstream interface nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 2.3\n" );
# endif

  if( project->subdom.npr > 1 )
  {
    SUBDOM* subdom = &project->subdom;
    INFACE* inface = subdom->inface;

    int df = eqs->dfcn;

    for( int s=0; s<subdom->npr; s++ )
    {
      int np = inface[s].np;                      // number of interface nodes

      if( np > 0  &&  s < subdom->pid )           // upstream interface nodes
      {
        int cnt = 0;

        for( int n=0; n<np; n++ )
        {
          NODE* nd = inface[s].node[n];

          if( !isFS(nd->flag, NODE::kDry) )       // nothing to send if node is dry...
          {
            SUB* sub = nd->sub;
            while( sub )
            {
              if( sub->no == s )  break;
              sub = sub->next;
            }

            if( !sub->dry )                       // ...or if the node is dry in
            {                                     // any adjacent subdomain
              for( int e=0; e<df; e++ )
              {
                int eqno = eqs->GetEqno( nd, e );

                if( eqno >= 0 )
                {
                  inface[s].send[cnt] = X[eqno];
                  cnt++;
                }
              }
            }
          }
        }

#       ifdef _MPI_DBG
        {
          char text[200];
          sprintf( text, " ### sending %d values to %d\n", cnt, s+1 );
          REPORT::rpt.Output( text );
        }
#       endif

        MPI_Send( &cnt, 1, MPI_INT, s, 1, MPI_COMM_WORLD );
        if( cnt )  MPI_Send( inface[s].send, cnt, MPI_DOUBLE, s, 2, MPI_COMM_WORLD );
      }
    }
  }

# endif
  ////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.4 MPI: backward solve for interior nodes

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     starting with 2.4\n" );
# endif

  for( int i=neq_up-1; i>=0; i-- )
  {
    for( int j=1; j<width[i]; j++ )
    {
      int eq = index[i][j];

      if( eq > i )
      {
        X[i] -= ILU[i][j] * X[eq];
      }
    }

    if( fabs(ILU[i][0]) < kZero )
      REPORT::rpt.Error( kParameterFault, "division by zero - EQS::ILU_solver(2)" );

    X[i] /= ILU[i][0];
  }

# ifdef kDebug
  {
    FILE* dbg;
    char  filename[80];

    MODEL* model  = project->M2D;
    GRID*  region = model->region;

    if( project->subdom.npr == 1 )
      sprintf( filename, "backX.dbg" );
    else
      sprintf( filename, "backX_%02d.dbg", project->subdom.pid );

    dbg = fopen( filename, "w" );
    fprintf( dbg, "%d\n", region->Getnp() );

    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE* nd = region->Getnode( n );
      for( int e=0; e<eqs->dfcn; e++ )
      {
        fprintf( dbg, "%5d %1d ", nd->Getname(), e );

        int eqno = eqs->GetEqno( nd, e );

        if( eqno >= 0 )
          fprintf( dbg, "%14.6le\n", X[eqno] );
        else
          fprintf( dbg, "%14.6le\n", 0.0 );
      }
    }

    fclose( dbg );
  }

  MPI_Barrier( MPI_COMM_WORLD );

# endif

# ifdef _MPI_DBG
  REPORT::rpt.Output( " (PRECO_ILUT::Solve)     backward solve finished\n" );
# endif
}
