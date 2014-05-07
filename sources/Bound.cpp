// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software.
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================

// ---------------------------------------------------------------------------------------
// determine boundary elements and set up boundary conditions
// ---------------------------------------------------------------------------------------

#include "Bcon.h"
#include "Node.h"
#include "Shape.h"
#include "Elem.h"
#include "Type.h"
#include "Subdom.h"

#include "Model.h"

//#define kDebug


void MODEL::Boundary()
{
  char text[200];

  // initializations ---------------------------------------------------------------------

  if( !boundList )
  {
    boundList = new ELEM* [region->Getnp()];
    if( !boundList )
      REPORT::rpt.Error( kMemoryFault, "can not allocate memory - MODEL::Boundary(1)" );
  }

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    CF( nd->flag, NODE::kBound );

    nd->mark = false;


    // remove slip boundary conditions ---------------------------------------------------

    if( isFS(nd->bc.kind, BCON::kAutoSlip) )
    {
      CF( nd->bc.kind, BCON::kAutoSlip | BCON::kFixU | BCON::kFixV );
    }

    if( isFS(nd->bc.kind, BCON::kAutoKD) )
    {
      CF( nd->bc.kind, BCON::kAutoKD | BCON::kFixK | BCON::kFixD );
    }
  }


  // determine boundary midside nodes ----------------------------------------------------

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    // midside nodes which are connected to only one element are boundary nodes ----------

    if( nd->noel == 1  &&  isFS(nd->flag, NODE::kMidsNode) )
    {
      if( !isFS(nd->flag, NODE::kInface) )  SF( nd->flag, NODE::kBound );
    }
  }


  // determine number of boundary elements: nb -------------------------------------------

  int nb = 0;

  for( int n=0; n<region->Getnp(); n++ )
  {
    if( isFS(region->Getnode(n)->flag,NODE::kBound) )  nb++;
  }


  // allocate memory for boundary elements -----------------------------------------------

  bound->Free();
  bound->Alloc( 0, nb );


  // set up boundary elements ------------------------------------------------------------

  int be = 0;       // counter for boundary elements

  for( int re=0; re<region->Getne(); re++ )
  {
    ELEM* el = region->Getelem(re);

    if( isFS(el->flag, ELEM::kDry) )  continue;

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    for( int i=ncn; i<nnd; i++ )
    {
      // check, if el->nd[i] is a midside boundary node ----------------------------------

      if( isFS(el->nd[i]->flag, NODE::kBound) )
      {
        ELEM* bd = bound->Getelem(be);

        boundList[el->nd[i]->Getno()] = bd;

        int left = i - ncn;
        int rght = (left + 1) % ncn;

        bd->nd[0] = el->nd[left];           // corner nodes
        bd->nd[1] = el->nd[rght];
        bd->nd[2] = el->nd[i];              // midside node

        SF( bd->nd[0]->flag, NODE::kBound );
        SF( bd->nd[1]->flag, NODE::kBound );


        // set shape specifications ------------------------------------------------------

        bd->Setshape( kLine );
        bd->Setname( el->Getname() );

        SF( bd->flag, ELEM::kBound );

        bd->type     = el->type;
        bd->areaFact = 1.0;

        be++;
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // communicate boundary nodes

//# ifdef _MPI_DBG
//  REPORT::rpt.Output( " (MODEL::Boundary)       communication of boundary nodes", 1 );
//# endif

# ifdef _MPI_
  if( subdom->npr > 1 )
  {
    INFACE* inface = subdom->inface;

    // loop on all interfaces: exchange bound flag ---------------------------------------
    for( int s=0; s<subdom->npr; s++ )
    {
      MPI_Status status;

      int npinf = inface[s].np;

      if( npinf > 0 )
      {
        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];

          if( isFS(nd->flag, NODE::kBound) )  inface[s].sia1[n] = true;
          else                                inface[s].sia1[n] = false;
        }

        MPI_Sendrecv( inface[s].sia1, npinf, MPI_CHAR, s, 1,
                      inface[s].ria1, npinf, MPI_CHAR, s, 1,
                      MPI_COMM_WORLD, &status );

        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];
          if( inface[s].ria1[n] )  SF( nd->flag, NODE::kBound );
        }
      }
    }
  }
# endif
  ////////////////////////////////////////////////////////////////////////////////////////


  // -------------------------------------------------------------------------------------
  // count for newly required boundary conditions
  // note: (sc, 30.10.2004)
  // a boundary condition is needed for marsh-nodes in case of dry-rewet-method 3

  int nbc = 0;

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    if( isFS(nd->flag, NODE::kBound) )  nbc++;
  }

  sprintf( text,"\n (MODEL::Boundary)       number of boundary elements: %d\n", nb );
  REPORT::rpt.Output( text, 3 );


# ifdef kDebug
{
  int pid;
  MPI_Comm_rank( MPI_COMM_WORLD, &pid );

  char fname[40];
  sprintf( fname, "bound_%02d.inp", pid+1 );

  FILE* id = fopen( fname, "w" );

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    CF( nd->flag, NODE::kMarker );
  }


  for( int e=0; e<bound->Getne(); e++ )
  {
    ELEM* el = bound->Getelem(e);

    for( int i=0; i<el->getnnd(); i++ )
    {
      SF( el->nd[i]->flag, NODE::kMarker );
    }
  }


  int j = 0;

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    if( isFS(nd->flag, NODE::kMarker) )  j++;
  }


  fprintf( id, "%6d  %6d   0   0   0\n", j, nb );


  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    if( isFS(nd->flag, NODE::kMarker) )
    {
      fprintf( id, "%6d  %17.9le  %17.9le  %17.9le\n",
                   nd->Getname(), nd->x, nd->y, nd->z );
    }
  }

  for( int e=0; e<bound->Getne(); e++ )
  {
    ELEM* el = bound->Getelem(e);

    fprintf( id, "%6d  %3d  line   %6d  %6d  %6d\n",
                 el->Getname(), TYPE::getid(el->type),
                 el->nd[0]->Getname(), el->nd[1]->Getname(), el->nd[2]->Getname() );
  }

  fclose( id );
}
# endif
}


ELEM* MODEL::Getbound( int no )
{
  return boundList[no];
}
