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
// dry/rewet and divergence free flow field cycle
// -------------------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Shape.h"
#include "Memory.h"
#include "Model.h"
#include "Subdom.h"

#include "Project.h"

//#define kDebug_1
//#define kDebug_2


void MODEL::DoDryRewet( PROJECT* project, int* dried, int* wetted )
{
  DRYREW *dryRew = &region->dryRew;

  region->Connection( 0L );

  int del = 0;
  int wel = 0;

  if( dryRew->method == 1 )
  {
    // mark nodes and elements to be rewetted ------------------------------------------------------
    wel = region->Rewet( dryRew->rewetLimit, dryRew->rewetPasses, project );

    // mark dry nodes and elements -----------------------------------------------------------------
    del = region->Dry( dryRew->dryLimit, dryRew->countDown );
  }
  else if( dryRew->method == 2 )
  {
    region->DryRewet( dryRew->dryLimit, dryRew->rewetLimit, dryRew->countDown, &del, &wel );
  }
  else if( dryRew->method == 3 )
  {
    region->RewetDry( dryRew->dryLimit, dryRew->rewetLimit, dryRew->countDown, &del, &wel );
  }
/*
  // future work ...
  else if( dryRew->method == 4 )
  {
    // determine dynamic boundary ------------------------------------------------------------------
    region.DynamicBound( np, node, *elem, dryRew->dryLimit, project );
  }
*/

  //////////////////////////////////////////////////////////////////////////////////////////////////

# ifdef _MPI_
  del = project->subdom.Mpi_sum( del );
  wel = project->subdom.Mpi_sum( wel );
# endif

  char text [200];

  sprintf( text, "\n (MODEL::DoDryRewet)     %d elements have got dry\n", del );
  REPORT::rpt.Output( text, 3 );

  sprintf( text, "\n (MODEL::DoDryRewet)     %d elements have got wet\n", wel );
  REPORT::rpt.Output( text, 3 );

  int dryRewFlag = del + wel;

  if( dryRewFlag )
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // exchange information on dry nodes on interfaces
#   ifdef _MPI_
    if( project->subdom.npr > 1 )
    {
      MPI_Comm_Dry( project, true );

      //////////////////////////////////////////////////////////////////////////////////////////////
      // reset wetted area, in case of dry nodes that are not dry in adjacent subdomains
      // added on 20.04.2006, sc

      if( dryRew->method == 2  ||  dryRew->method == 3 )
      {
        int wetted = 0;

        for( int n=0; n<region->Getnp(); n++ )
        {
          NODE* nd = region->Getnode(n);

          nd->mark = false;

          double H = nd->v.S - nd->zor;

          if( H < dryRew->dryLimit )
          {
            SF( nd->flag, NODE::kDry );
          }

          else if( isFS(nd->flag, NODE::kDry) )
          {
            nd->mark = true;

            CF( nd->flag, NODE::kDry );
            wetted++;
          }

          CF( nd->flag, NODE::kMarsh );
          nd->z = nd->zor;
        }


        if( dryRew->method == 2 )
        {
          region->DryRewet( dryRew->dryLimit, dryRew->rewetLimit, dryRew->countDown, &del, &wel );
        }
        else if( dryRew->method == 3 )
        {
          region->RewetDry( dryRew->dryLimit, dryRew->rewetLimit, dryRew->countDown, &del, &wel );
        }

        ////////////////////////////////////////////////////////////////////////////////////////////

        MPI_Comm_Dry( project, false );

        wetted = project->subdom.Mpi_sum( wetted );

        sprintf( text, "\n (MODEL::DoDryRewet)     %d interface nodes have got wet\n", wetted );
        REPORT::rpt.Output( text, 3 );
      }
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // reset interface flags: kInface, kInface_DN and kInface_UP; this is
    // necessary for the correct ordering of equations in EQS::ResetEqOrder()

    project->subdom.SetInface( region );

#   endif  //  #ifdef _MPI_
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // determine 1D boundary elements and ----------------------------------------------------------
    // set up slip velocity boundary conditions

    Initialize();

    SetNormal();
    SetRotation();

    region->SetSlipFlow();

    REPORT::rpt.PrintTime( 3 );

//  ================================================================================================
#   ifdef kDebug_1
    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE* ndbg = region->Getnode(n);

      switch( ndbg->Getname() )
      {
        case 3654:
          REPORT::rpt.Message( "\n" );
          REPORT::rpt.Message( "### NODE %6d: inface    = %d\n", ndbg->Getname(), ndbg->flag&NODE::kInface );
          REPORT::rpt.Message( "###            : inface_dn = %d\n", ndbg->flag&NODE::kInface_DN );
          REPORT::rpt.Message( "###            : inface_up = %d\n", ndbg->flag&NODE::kInface_UP );
          REPORT::rpt.Message( "\n" );
          REPORT::rpt.Message( "###            : dry       = %d\n", ndbg->flag&NODE::kDry );
          REPORT::rpt.Message( "###            : marsh     = %d\n", ndbg->flag&NODE::kMarsh );
          REPORT::rpt.Message( "###            : H         = %f\n", ndbg->v.S - ndbg->zor );
          REPORT::rpt.Message( "\n" );
          REPORT::rpt.Message( "###            : bound     = %d\n", ndbg->flag&NODE::kBound );
          REPORT::rpt.Message( "###            : inflow    = %d\n", ndbg->flag&NODE::kInlet );
          REPORT::rpt.Message( "###            : outflow   = %d\n", ndbg->flag&NODE::kOutlet );
          REPORT::rpt.Message( "###            : rotat     = %d\n", ndbg->flag&NODE::kRotat );
          REPORT::rpt.Message( "###            : noMoment  = %d\n", ndbg->flag&NODE::kNoMoment );
          REPORT::rpt.Message( "\n" );
          REPORT::rpt.Message( "###            : countDown = %d\n", ndbg->countDown );

          {
            SUB* sub = ndbg->sub;
            while( sub )
            {
              REPORT::rpt.Message( "### SUBDOM %4d: dry       = %d\n", sub->no+1, sub->dry );
              sub = sub->next;
            }
          }
          break;
      }
    }
#   endif
//  ================================================================================================
  }

  else
  {
    // set up structure for connection of nodes to elements ----------------------------------------
    // dry elements are not taken into consideration

    region->Connection( ELEM::kDry );
  }

  if( dried )  *dried  = del;
  if( wetted ) *wetted = wel;

  region->firstDryRew = false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void MODEL::MPI_Comm_Dry( PROJECT* project, int initialize )
{
# ifdef _MPI_

  DRYREW* dryRew = &region->dryRew;

  SUBDOM* subdom = &project->subdom;
  INFACE* inface = subdom->inface;


  // loop on all interfaces: exchange dry flag -----------------------------------------------------
  // set sub->dry flag for dry interface nodes in adjacent subdomains

  for( int s=0; s<subdom->npr; s++ )
  {
    MPI_Status status;

    int np = inface[s].np;

    if( np > 0 )
    {
      // initialize the flag "nd->sub->dry" --------------------------------------------------------

      for( int n=0; n<np; n++ )
      {
        NODE* nd = inface[s].node[n];

        SUB* sub = nd->sub;
        while( sub )
        {
          if( sub->no == s )  sub->dry = false;
          sub = sub->next;
        }

        if( isFS(nd->flag, NODE::kDry) )  inface[s].sia1[n] = true;
        else                              inface[s].sia1[n] = false;
      }


      // exchange data on interface nodes ----------------------------------------------------------

      MPI_Sendrecv( inface[s].sia1, np, MPI_CHAR, s, 1,
                    inface[s].ria1, np, MPI_CHAR, s, 1,
                    MPI_COMM_WORLD, &status );


      // set the flag "nd->sub->dry" for dry nodes in adjacent subdomain and -----------------------
      // adjust water elevation (necessary for just wetted interface nodes)

      for( int n=0; n<np; n++ )
      {
        NODE* nd = inface[s].node[n];

        if( inface[s].ria1[n] )
        {
          SUB* sub = nd->sub;
          while( sub )
          {
            if( sub->no == s )  sub->dry = true;
            sub = sub->next;
          }
        }
      }
    }
  }


  // average data across domains (necessary for wetting interfaces) --------------------------------

  if( initialize )
  {
    MPI_Status status;

    int  rgnp = region->Getnp();
    int* cnt  = (int*) MEMORY::memo.Array_nd( rgnp );

    for( int n=0; n<rgnp; n++ )  cnt[n] = 1;

    for( int s=0; s<subdom->npr; s++ )
    {
      int npinf = inface[s].np;

      if( npinf > 0 )
      {
        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];

          inface[s].send[3*n]   = nd->v.U;
          inface[s].send[3*n+1] = nd->v.V;
          inface[s].send[3*n+2] = nd->v.S;
        }
      }
    }

    for( int s=0; s<subdom->npr; s++ )
    {
      int npinf = inface[s].np;

      if( npinf > 0 )
      {
        MPI_Sendrecv( inface[s].send, 3*npinf, MPI_DOUBLE, s, 1,
                      inface[s].recv, 3*npinf, MPI_DOUBLE, s, 1,
                      MPI_COMM_WORLD, &status );

        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];
          SUB* sub = nd->sub;

          while( sub )
          {
            if( sub->no == s )
            {
              if( !sub->dry )
              {
                nd->v.U += inface[s].recv[3*n];
                nd->v.V += inface[s].recv[3*n+1];
                nd->v.S += inface[s].recv[3*n+2];
                cnt[nd->Getno()]++;
              }
            }
            sub = sub->next;
          }
        }
      }
    }

    for( int n=0; n<rgnp; n++ )
    {
      if( cnt[n] > 1 )
      {
        NODE* nd = region->Getnode(n);

        nd->v.U /= cnt[n];
        nd->v.V /= cnt[n];
        nd->v.S /= cnt[n];
      }
    }

    MEMORY::memo.Detach( cnt );
  }

# ifdef kDebug_2
  {
    char text[200];

    sprintf( text, "### DEBUGGING: List of dry interface nodes...\n\n" );
    REPORT::rpt.Output( text, 1 );

    for( int s=0; s<subdom->npr; s++ )
    {
      int np = inface[s].np;

      if( np > 0 )
      {
        sprintf( text, "interface to domain %d\n", s+1 );
        REPORT::rpt.Output( text, 1 );

        for( int n=0; n<np; n++ )
        {
          NODE* nd = inface[s].node[n];

          SUB* sub = nd->sub;
          while( sub )
          {
            if( sub->no == s )
            {
              if( sub->dry )
              {
                if( isFS(nd->flag, NODE::kDry) )
                {
                  sprintf( text, "dry interface node %5d: dry\n", nd->Getname() );
                  REPORT::rpt.Output( text, 1 );
                }
                else
                {
                  sprintf( text, "wet interface node %5d: dry\n", nd->Getname() );
                  REPORT::rpt.Output( text, 1 );
                }
              }
              else
              {
                if( isFS(nd->flag, NODE::kDry) )
                {
                  sprintf( text, "dry interface node %5d: wet\n", nd->Getname() );
                  REPORT::rpt.Output( text, 1 );
                }
                else
                {
                  sprintf( text, "wet interface node %5d: wet\n", nd->Getname() );
                  REPORT::rpt.Output( text, 1 );
                }
              }
            }

            sub = sub->next;
          }
        }
      }
    }
  }
# endif  // #ifdef kDebug_2
# endif  // #ifdef _MPI_
}
