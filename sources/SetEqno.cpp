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
// set equation numbers according to nodes and elements
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Memory.h"
#include "Elem.h"
#include "Node.h"
#include "Model.h"
#include "Project.h"

#include "Eqs.h"


//#define DEBUG 1


void EQS::SetEqno( MODEL*   model,
                   int      neqcn,       // equations at corner nodes
                   int      neqmn,       // equations at midside nodes
                   int      neqel,       // equations at elements
                   unsigned
                   int*     fix,
                   int      elemKind )
{
  char text[80];

  int    ne   = model->ne;
  ELEM** elem = model->elem;
  int    np   = model->np;
  NODE** node = model->node;

  int    rgnp = model->region->Getnp();
  int    rgne = model->region->Getne();


  // allocate memory for equation numbers ------------------------------------------------

  if( !nodeEqno  &&  dfcn > 0 )
  {
    nodeEqno = new int* [ dfcn ];
    int* buf = new int  [ dfcn * rgnp ];

    if( !nodeEqno || !buf )
      REPORT::rpt.Error( "can not allocate memory - EQS::SetEqno(1)" );

    for( int i=0; i<dfcn; i++ )  nodeEqno[i] = buf + i * rgnp;
  }

  if( !elemEqno  &&  dfel > 0 )
  {
    elemEqno = new int* [ dfel ];
    int* buf = new int  [ dfel * rgne ];

    if( !elemEqno || !buf )
      REPORT::rpt.Error( "can not allocate memory - EQS::SetEqno(2)" );

    for( int i=0; i<dfel; i++ )  elemEqno[i] = buf + i * rgne;
  }


  // initialization ----------------------------------------------------------------------

  for( int i=0; i<dfcn; i++ )
  {
    for( int j=0; j<rgnp; j++ )  nodeEqno[i][j] = -1;
  }

  for( int i=0; i<dfel; i++ )
  {
    for( int j=0; j<rgne; j++ )  elemEqno[i][j] = -1;
  }

  // determine equation numbers at nodes -------------------------------------------------
/*
 * MPI: neq_up and neq_dn is now set in method EQS:RestEqOrder()
 *
  neq_up = -1;
  neq_dn = -1;
*/

  neq = 0;

  for( int n=0; n<np; n++ )
  {
    BCON* bc = &node[n]->bc;
/*
    if( neq_up < 0  &&  isFS(node[n]->flag, NODE::kInface_UP) )  neq_up = neq;
    if( neq_dn < 0  &&  isFS(node[n]->flag, NODE::kInface_DN) )  neq_dn = neq;
*/
    if( isFS(node[n]->flag, NODE::kCornNode) )
    {
      for( int j=0; j<neqcn; j++ )
      {
        if( isFS(bc->kind,fix[j]) )
        {
          continue;
        }

        else
        {
          nodeEqno[j][node[n]->Getno()] = neq;
          neq++;
        }
      }
    }

    else
    {
      for( int j=0; j<neqmn; j++ )
      {
        if( bc  &&  isFS(bc->kind, fix[j]) )
        {
          continue;
        }

        else
        {
          nodeEqno[j][node[n]->Getno()] = neq;
          neq++;
        }
      }
    }
  }
/* MPI: obsolete
       if( neq_up < 0  &&  neq_dn < 0 )  neq_up = neq_dn = neq;
  else if( neq_up < 0 )                  neq_up = neq_dn;
  else if( neq_dn < 0 )                  neq_dn = neq;
*/

  // determine equation numbers at elements ----------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = elem[e];

    for( int j=0; j<neqel; j++ )
    {
      if( isFS(el->flag, ELEM::kBound) )
      {
        continue;
      }

      else
      {
        elemEqno[j][el->Getno()] = neq;
        neq++;
      }
    }
  }


  // report total number of simultaneous equations ---------------------------------------
  sprintf( text, "\n (EQS::SetEqno)          %s %d\n",
                 "number of equations is", neq );
  REPORT::rpt.Output( text, 3 );


  // -------------------------------------------------------------------------------------
  // determine last occurence of equations (frontal solver) an reset equation order
  // ResetEqOrder() is also necessary for proper working of ILU preconditioning

  LastEquation( model, elemKind );

  ResetEqOrder( model );


  // set up list of node pointers for equation numbers

  if( !eqnoNode || !eqid )
  {
    eqid     = new int  [kSimDF * rgnp];
    eqnoNode = new NODE*[kSimDF * rgnp];
    if( !eqid || !eqnoNode )
      REPORT::rpt.Error( "can not allocate memory - EQS::SetEqno(3)" );
  }

  for( int n=0; n<np; n++ )
  {
    for( int j=0; j<dfcn; j++ )
    {
      int e = nodeEqno[j][node[n]->Getno()];

      if( e >= 0 )
      {
        eqid[e]     = j;
        eqnoNode[e] = node[n];
      }
    }
  }
}


// ---------------------------------------------------------------------------------------
// determine last occurence of equations during element assembling

void EQS::LastEquation( MODEL* model, long elemKind )
{
  int    ne   = model->ne;
  ELEM** elem = model->elem;


  // allocate temporary memory for counter vector ----------------------------------------

  int* counter = (int*) MEMORY::memo.Array_eq( neq, "EQS::LastEquation(1)" );


  // initialize counter array ------------------------------------------------------------

  for( int i=0; i<neq; i++ ) counter[i] = 0;


  // initialize the elem->isLast array ---------------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = elem[e];

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      for( int j=0; j<dfcn; j++ )
      {
        el->isLast[j][i] = 0;
      }
    }
  }


  // determine the number of elements associated to equations ----------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = elem[e];

    int eqno, nnd;

    if( isFS(el->flag, elemKind)  &&  !isFS(el->flag, ELEM::kDry) )
    {
      nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        for( int j=0; j<dfcn; j++ )
        {
          eqno = GetEqno( el->nd[i], j );

          if( eqno >= 0 )
            counter[eqno]++;
        }
      }
    }
  }


  // determine the last occurence of equations -------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = elem[e];

    int eqno, nnd;

    if( isFS(el->flag, elemKind)  &&  !isFS(el->flag, ELEM::kDry) )
    {
      nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        for( int j=0; j<dfcn; j++ )
        {
          eqno = GetEqno( el->nd[i], j );

          if( eqno >= 0 )
          {
            counter[eqno]--;

            if( counter[eqno] <  0 )
              REPORT::rpt.Error ("unexpected: wrong count of equations - lastEquation (2)");

            if( counter[eqno] == 0 )
              el->isLast[j][i] = 1;
          }
        }
      }
    }
  }

  MEMORY::memo.Detach( counter );

  REPORT::rpt.Output("\n (EQS::LastEquation)     last occurence of equations determined\n", 4 );


# ifdef DEBUG
  {
  int   count;
  FILE* id;

  id = fopen( "lastEquation.Report", "w" );

  count = 0;

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = elem[e];

    int nnd = el->getnnd();

    if( isFS(el->flag, ELEM::kRegion) )
      fprintf( id, "--- region element %d\n", el->Getname() );

    else
      fprintf( id, "--- boundary element %d\n", el->Getname() );


    for( int i=0; i<nnd; i++ )
    {
      fprintf( id, "    node %6d .......................\n", el->nd[i]->Getname() );

      for( int j=0; j<dfcn; j++ )
      {
        int eqno = GetEqno( el->nd[i], j );

        if( eqno >= 0 )
        {
          if( el->isLast[j][i] )
          {
            fprintf( id, "               %6d - last occurence\n", eqno );
            count++;
          }

          else
          {
            fprintf( id, "               %6d\n", eqno );
          }
        }
      }
    }

    fprintf( id, "\n\n" );
  }

  fprintf( id, "total number of equations is: %d\n", count );

  fclose( id );
  }
# endif

}


void EQS::ResetEqOrder( MODEL* model )
{
  int* eqOrd = (int*) MEMORY::memo.Array_eq( neq, "EQS::ResetEqOrder(1)" );


  ////////////////////////////////////////////////////////////////////////////////////////
  // MPI: keep ordering of equations
  // ordering: 1. interior nodes             (nd->flag & kInface)
  //           2. upstream interface nodes   (nd->flag & kInface_UP)
  //           3. downstream interface nodes (nd->flag & kInface_DN)

  int n = 0;

  for( int e=0; e<model->ne; e++ )
  {
    ELEM* el = model->elem[e];

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      NODE* nd = el->nd[i];

      if( !isFS(nd->flag, NODE::kInface) )
      {
        for( int j=0; j<dfcn; j++ )
        {
          int eqno = GetEqno( nd, j );

          if( eqno >= 0  &&  el->isLast[j][i] )
          {
            if( n >= neq )
              REPORT::rpt.Error( kUnexpectedFault, "%s - EQS::resetEqOrder(2)",
                                         "wrong number of equations" );

            eqOrd[eqno] = n;
            n++;
          }
        }
      }
    }

    for( int j=0; j<dfel; j++ )
    {
      int eqno = GetEqno( el, j );

      if( eqno >= 0 )
      {
        if( n >= neq )
          REPORT::rpt.Error( kUnexpectedFault, "%s - EQS::resetEqOrder(3)",
                                     "wrong number of equations" );
        eqOrd[eqno] = n;
        n++;
      }
    }
  }

  neq_up = n;
  neq_dn = n;


  for( int e=0; e<model->ne; e++ )
  {
    ELEM* el = model->elem[e];

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      NODE* nd = el->nd[i];

      if( isFS(nd->flag, NODE::kInface) && isFS(nd->flag, NODE::kInface_UP) )
      {
        for( int j=0; j<dfcn; j++ )
        {
          int eqno = GetEqno( nd, j );

          if( eqno >= 0  &&  el->isLast[j][i] )
          {
            if( n >= neq )
              REPORT::rpt.Error( kUnexpectedFault, "%s - EQS::resetEqOrder(2)",
                                         "wrong number of equations" );

            eqOrd[eqno] = n;
            n++;
          }
        }
      }
    }
  }

  neq_dn = n;


  for( int e=0; e<model->ne; e++ )
  {
    ELEM* el = model->elem[e];

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      NODE* nd = el->nd[i];

      if( isFS(nd->flag, NODE::kInface) && isFS(nd->flag, NODE::kInface_DN) )
      {
        for( int j=0; j<dfcn; j++ )
        {
          int eqno = GetEqno( nd, j );

          if( eqno >= 0  &&  el->isLast[j][i] )
          {
            if( n >= neq )
              REPORT::rpt.Error( kUnexpectedFault, "%s - EQS::resetEqOrder(2)",
                                 "wrong number of equations" );

            eqOrd[eqno] = n;
            n++;
          }
        }
      }
    }
  }


  for( int i=0; i<model->np; i++ )
  {
    NODE* nd = model->node[i];

    for( int j=0; j<dfcn; j++ )
    {
      int eqno = GetEqno( nd, j );

      if( eqno >= 0 )  nodeEqno[j][nd->Getno()] = eqOrd[eqno];
    }
  }


  for( int e=0; e<model->ne; e++ )
  {
    ELEM* el = model->elem[e];

    for( int j=0; j<dfel; j++ )
    {
      int eqno = GetEqno( el, j );

      if( eqno >= 0 )  elemEqno[j][el->Getno()] = eqOrd[eqno];
    }
  }

  MEMORY::memo.Detach( eqOrd );

  //char text[100];
  //sprintf( text, "\n interface equations in SetEqno(): neq_up = %d / neq_dn = %d / neq = %d\n",
  //         neq_up, neq_dn, neq );
  //REPORT::rpt.Output( text );
}
