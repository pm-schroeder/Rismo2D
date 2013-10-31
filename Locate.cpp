// =================================================================================================
//
// Copyright (C) 1992-2013  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// =================================================================================================
//
// -------------------------------------------------------------------------------------------------
// set location of boundary elements and nodes
// -------------------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Bcon.h"
#include "Node.h"
#include "Elem.h"

#include "Model.h"


void MODEL::SetLocation()
{
  // -----------------------------------------------------------------------------------------------
  // set location of inflow, outflow and openflow elements

  for( int e=0; e<bound->Getne(); e++ )
  {
    ELEM* el = bound->Getelem(e);

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    int inlet   = true;
    int outlet  = true;
    int openBnd = true;

    for( int i=0; i<ncn; i++ )
    {
      BCON* bc = &el->nd[i]->bc;

      if( !isFS(bc->kind, BCON::kInlet) )    inlet   = false;
      if( !isFS(bc->kind, BCON::kOutlet) )   outlet  = false;
      if( !isFS(bc->kind, BCON::kOpenBnd) )  openBnd = false;
    }

    if( inlet )    SF( el->flag, ELEM::kInlet );
    if( outlet )   SF( el->flag, ELEM::kOutlet );
    if( openBnd )  SF( el->flag, ELEM::kOpenBnd );

    for( int i=0; i<nnd; i++ )
    {
      if( inlet )    SF( el->nd[i]->flag, NODE::kInlet );
      if( outlet )   SF( el->nd[i]->flag, NODE::kOutlet );
      if( openBnd )  SF( el->nd[i]->flag, NODE::kOpenBnd );
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // communicate location of inlet/outlet nodes
  // !!! to do: check if this is necessary !!! (sc, 19.07.2005)

//# ifdef _MPI_
//  if( subdom->npr > 1 )
//  {
//    INFACE* inface = subdom->inface;
//
//    // loop on all interfaces: exchange NODE::flag -----------------------------------------------
//    for( int s=0; s<subdom->npr; s++ )
//    {
//      MPI_Status status;
//
//      int npinf = inface[s].np;
//
//      if( npinf > 0 )
//      {
//        for( int n=0; n<npinf; n++ )
//        {
//          NODE* nd = inface[s].node[n];
//
//          inface[s].sia1[n] = 0;
//
//          if( isFS(nd->flag, NODE::kInlet) )   inface[s].sia1[n] |= 1;
//          if( isFS(nd->flag, NODE::kOutlet) )  inface[s].sia1[n] |= 2;
//        }
//
//        MPI_Sendrecv( inface[s].sia1, npinf, MPI_CHAR, s, 1,
//                      inface[s].ria1, npinf, MPI_CHAR, s, 1,
//                      MPI_COMM_WORLD, &status );
//
//        for( int n=0; n<npinf; n++ )
//        {
//          NODE* nd = inface[s].node[n];
//
//          if( inface[s].ria1[n] & 1 )  SF( nd->flag, NODE::kInlet );
//          if( inface[s].ria1[n] & 2 )  SF( nd->flag, NODE::kOutlet );
//        }
//      }
//    }
//  }
//# endif
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // -----------------------------------------------------------------------------------------------
  // Set kNoMoment flag for boundary corner nodes connected to one element only (convex nodes).
  // Nodes on open boundary (inflow/outflow) are not convex.
  // MPI: nodes on interfaces are not convex.

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    CF( nd->flag, NODE::kNoMoment );

    // ### test - 09.01.2008 #######################################################################
    // ### set flag "kNoMoment" for marsh nodes
    // if( isFS(nd->flag, NODE::kMarsh) ) SF( nd->flag, NODE::kNoMoment );
    // ### end of test - 09.01.2008 ################################################################
    //
    // ### new test - 15.02.2010 ###################################################################
    // ### with some restrictions concerning in-, out- and openflow-boundaries
    // ### -> convergence of Newton-Raphson and BiCGSTAB gets worse
    // if( isFS(nd->flag, NODE::kMarsh)
    //     && !isFS(nd->flag, NODE::kInlet)
    //     && !isFS(nd->flag, NODE::kOutlet)
    //     && !isFS(nd->flag, NODE::kOpenBnd) )
    // {
    //   SF( nd->flag, NODE::kNoMoment );
    // }
    // ### end of test - 15.02.2010 ################################################################
    //
    // ### new test - 26.02.2013 ###################################################################
    // ### a problem occurs in elements downstrem to the inlet, where the boundary condition
    // kNoMoment inhibits flowing and retains the water leading to an increasing water elevation.
    // if( isFS(nd->flag, NODE::kMarsh) && !isFS(nd->flag, NODE::kInlet) )
    // {
    //   SF( nd->flag, NODE::kNoMoment );
    // }
    // ### end of test - 26.02.2013 ################################################################
  }


  for( int e=0; e<bound->Getne(); e++ )
  {
    ELEM* el = bound->Getelem(e);

    int ncn = el->Getncn();

    for( int i=0; i<ncn; i++ )
    {
      // 1. treat convex nodes ---------------------------------------------------------------------
      if( el->nd[i]->noel == 1
          && !isFS(el->nd[i]->flag, NODE::kInface)
          && !isFS(el->nd[i]->flag, NODE::kInlet)
          && !isFS(el->nd[i]->flag, NODE::kOutlet)
          && !isFS(el->nd[i]->flag, NODE::kOpenBnd) )
      {
        SF( el->nd[i]->flag, NODE::kNoMoment );
      }
    }
  }

  char text[200];
  sprintf( text, "\n (MODEL::SetLocation)    %s\n",
                 "location of boundary nodes (inlet, outlet, ...)" );
  REPORT::rpt.Output( text, 4 );
}
