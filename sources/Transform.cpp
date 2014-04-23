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
// methods in this module:
//
//   SetNormal()   -  determine normal direction on boundary nodes
//
//   SetRotation() -  determine rotation matrix for boundary nodes
//
// -------------------------------------------------------------------------------------------------

#include "Defs.h"

#include "Bcon.h"
#include "Shape.h"
#include "Node.h"
#include "Elem.h"
#include "Subdom.h"

#include "Model.h"

// defines for debugging output in method MODEL::SetNormal()
// #define kDebug_1_1
// #define kDebug_1_2

// defines for debugging output in method MODEL::SetRotation()
// #define kDebug_2


void MODEL::SetNormal()
{
  // -----------------------------------------------------------------------------------------------
  // initialization

  for( int i=0; i<region->Getnp(); i++ )
  {
    BCON* bc = &region->Getnode(i)->bc;

    bc->nx   = bc->ny   = 0.0;    // normal direction on closed boundaries
    bc->niox = bc->nioy = 0.0;    // normal direction on open boundaries (inlet, outlet, ...)
  }


  // -----------------------------------------------------------------------------------------------
  // initialize boundary conditions, loop on all boundary elements

  for( int e=0; e<bound->Getne(); e++ )
  {
    ELEM* el = bound->Getelem(e);

    SHAPE* qShape = el->GetQShape();

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      double* dndx = qShape->dndx[i];

      BCON* bc = &el->nd[i]->bc;

      double nx = 0.0;
      double ny = 0.0;

      for( int j=0; j<nnd; j++ )
      {
        nx += 2.0 * (dndx[j] * el->nd[j]->y);
        ny -= 2.0 * (dndx[j] * el->nd[j]->x);
      }

      if(   isFS(el->flag,ELEM::kInlet)
         || isFS(el->flag,ELEM::kOutlet)
         || isFS(el->flag,ELEM::kOpenBnd) )
      {
        bc->niox += nx;
        bc->nioy += ny;
      }
      else
      {
        bc->nx += nx;
        bc->ny += ny;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
# ifdef _MPI_
  if( subdom->npr > 1 )
  {
    MPI_Status status;
    INFACE*    inface = subdom->inface;

    for( int s=0; s<subdom->npr; s++ )
    {
      int npinf = inface[s].np;

      for( int n=0; n<npinf; n++ )
      {
        NODE* nd = inface[s].node[n];

        if( isFS(nd->flag, NODE::kBound) )
        {
          inface[s].send[4*n]   = nd->bc.nx;
          inface[s].send[4*n+1] = nd->bc.ny;
          inface[s].send[4*n+2] = nd->bc.niox;
          inface[s].send[4*n+3] = nd->bc.nioy;
        }
        else
        {
          inface[s].send[4*n]   = 0.0;
          inface[s].send[4*n+1] = 0.0;
          inface[s].send[4*n+2] = 0.0;
          inface[s].send[4*n+3] = 0.0;
        }

      }
    }

    for( int s=0; s<subdom->npr; s++ )
    {
      int npinf = inface[s].np;

      MPI_Sendrecv( inface[s].send, 4*npinf, MPI_DOUBLE, s, 1,
                    inface[s].recv, 4*npinf, MPI_DOUBLE, s, 1,
                    MPI_COMM_WORLD, &status );

      for( int n=0; n<npinf; n++ )
      {
        NODE* nd = inface[s].node[n];

        if( isFS(nd->flag, NODE::kBound) )
        {
          nd->bc.nx   += inface[s].recv[4*n];
          nd->bc.ny   += inface[s].recv[4*n+1];
          nd->bc.niox += inface[s].recv[4*n+2];
          nd->bc.nioy += inface[s].recv[4*n+3];

          //////////////////////////////////////////////////////////////////////////////////////////
#         ifdef kDebug_1_1
          REPORT::rpt.Message( "\n" );
          REPORT::rpt.Message( "communication of normal directions...\n" );
          REPORT::rpt.Message( "### NODE %6d: {n_%03d} = {%12.6lf; %12.6lf}\n",
                               nd->Getname(), subdom->pid+1, nd->bc.nx, nd->bc.ny );
          REPORT::rpt.Message( "         %6s: {n_%03d} = {%12.6lf; %12.6lf}\n",
                               "  SEND",      s+1, inface[s].send[0], inface[s].send[1] );
          REPORT::rpt.Message( "         %6s: {n_%03d} = {%12.6lf; %12.6lf}\n",
                               "  RECV",      s+1, inface[s].recv[0], inface[s].recv[1] );
#         endif
          //////////////////////////////////////////////////////////////////////////////////////////
        }
      }
    }
  }
# endif
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // set normal vector to unit length --------------------------------------------------------------

  for( int i=0; i<region->Getnp(); i++ )
  {
    NODE* nd = region->Getnode(i);

    double length;

    if( isFS(nd->flag, NODE::kDry) )  continue;


    BCON* bc = &nd->bc;

    if( isFS(nd->flag, NODE::kBound) )
    {
      length = sqrt( bc->nx * bc->nx  +  bc->ny * bc->ny );

      if( fabs(length) < 1.0e-12 )
      {
        if(    !isFS(nd->flag, NODE::kInlet)
            && !isFS(nd->flag, NODE::kOutlet)
            && !isFS(nd->flag, NODE::kOpenBnd) )
        {
          char text[200];
          sprintf( text, "\n (MODEL::SetNormal)     %s %d\n",
                         "length of normal direction is zero at node",
                         nd->Getname() );
          REPORT::rpt.Output( text, 2 );
        }

        bc->nx = 0.0;
        bc->ny = 0.0;
      }

      else
      {
        bc->nx /= length;
        bc->ny /= length;
      }


      length = sqrt( bc->niox * bc->niox  +  bc->nioy * bc->nioy );

      if( fabs( length ) < 1.0e-12 )
      {
        bc->niox = 0.0;
        bc->nioy = 0.0;
      }

      else
      {
        bc->niox /= length;
        bc->nioy /= length;
      }
    }
  }

  char text[200];
  sprintf( text, "\n (MODEL::SetNormal)      %s\n",
                 "normal direction on boundary elements determined" );
  REPORT::rpt.Output( text, 4 );


# ifdef kDebug_1_2
  {
    char  name[80];

    sprintf( name, "SetNormal_%02d.report", subdom->pid+1 );
    FILE* id = fopen( name, "w" );

    fprintf( id, "Normal directions\n" );
    fprintf( id, "%d\n", region->Getnp() );

    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE* nd = region->Getnode(n);
      BCON* bc = &nd->bc;


      fprintf( id, "%6d  %12.5lf %12.5lf %12.5lf  %12.5lf\n",
                   nd->Getname(), bc->nx,   bc->ny,
                                  bc->niox, bc->nioy );
    }

    fclose( id );
  }
# endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void MODEL::SetRotation()
{
  // reset automatic slip flow conditions (BCON::kAutoSlip) for the new boundary -------------------

  for( int i=0; i<region->Getnp(); i++ )
  {
    NODE* nd = region->Getnode(i);

    CF( nd->flag, NODE::kRotat );

    BCON* bc = &nd->bc;

    bc->rot = 0;

    if( isFS(bc->kind, BCON::kInlet) )
    {
      CF( bc->kind, BCON::kFixU | BCON::kFixV );
    }

    if( isFS(bc->kind, BCON::kAutoSlip) )
    {
      CF( bc->kind, BCON::kAutoSlip | BCON::kFixU | BCON::kFixV );
    }

    // initialization of rotation matrices ---------------------------------------------------------

   bool bound = isFS(nd->flag, NODE::kBound);

//   if( region->dryRew.method >= 2 )
//   {
//     bound = isFS(nd->flag, NODE::kGridBnd)
//         || (isFS(nd->flag, NODE::kBound) && nd->z < nd->zor-1.0e-9 );
//         //|| (isFS(nd->flag, NODE::kBound) && isFS(nd->flag, NODE::kMarsh) );
//   }
//   else

   if(  bound  &&  !isFS(bc->kind, BCON::kInlet)
               &&  !isFS(bc->kind, BCON::kOutlet)
               &&  !isFS(bc->kind, BCON::kOpenBnd)
               &&  !isFS(bc->kind, BCON::kSlip)
               &&  !isFS(bc->kind, BCON::kSetUV) )
    {
      bc->rot = BCON::kTangentFlowRot;

      // set boundary condition type kAutoSlip and delete Y-momentum equation ----------------------
      SF( bc->kind, BCON::kAutoSlip | BCON::kFixV );
    }
  }


  // loop on all boundary elements -----------------------------------------------------------------
  // set boundary condition for nodes in outflow corner (for users convenience)

  for( int e=0; e<bound->Getne(); e++ )
  {
    ELEM* bd = bound->Getelem(e);

    // regard only those elements which do not lie on an outlet boundary
    if(    !isFS(bd->flag, ELEM::kOutlet)
        && !isFS(bd->flag, ELEM::kOpenBnd) )
    {
      // initialization of rotation matrices -------------------------------------------------------

      int ncn = bd->Getncn();

      for( int i=0; i<ncn; i++ )
      {
        BCON* bc = &bd->nd[i]->bc;

        // if one of the corner nodes belongs to an outlet boundary,
        // the flow direction will be set to the normal direction for the corner node
        // note: the normal direction has been computed in method MODEL:SetNormal()
        if( isFS(bc->kind, BCON::kOutlet) || isFS(bc->kind, BCON::kOpenBnd) )
        {
          if( !isFS(bc->kind, BCON::kSlip) && !isFS(bc->kind, BCON::kSetUV) )
          {
            bc->rot = BCON::kNormalFlowRot;

            // set boundary condition type kAutoSlip and delete Y-momentum equation ----------------
            SF( bc->kind, BCON::kAutoSlip | BCON::kFixV );
          }
        }
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // communication of BCON::kAutoSlip
# ifdef _MPI_
  if( subdom->npr > 1 )
  {
    INFACE* inface = subdom->inface;

    for( int s=0; s<subdom->npr; s++ )
    {
      MPI_Status status;

      int npinf = inface[s].np;

      for( int n=0; n<npinf; n++ )
      {
        NODE* nd = inface[s].node[n];

        if( isFS(nd->bc.kind, BCON::kAutoSlip) )  inface[s].sia1[n] = true;
        else                                      inface[s].sia1[n] = false;
      }

      MPI_Sendrecv( inface[s].sia1, npinf, MPI_CHAR, s, 1,
                    inface[s].ria1, npinf, MPI_CHAR, s, 1,
                    MPI_COMM_WORLD, &status );

      for( int n=0; n<npinf; n++ )
      {
        NODE* nd = inface[s].node[n];

        if( inface[s].ria1[n] )
        {
          SF( nd->bc.kind, BCON::kAutoSlip | BCON::kFixV );
        }
      }
    }
  }
# endif
  //////////////////////////////////////////////////////////////////////////////////////////////////


  // -----------------------------------------------------------------------------------------------
  // Set the kNoMoment flag for boundary corner nodes with a zero normal direction.
  // This may occure if two marsh elements are connected over one corner node.

  for( int i=0; i<region->Getnp(); i++ )
  {
    NODE* nd = region->Getnode(i);
    BCON* bc = &nd->bc;

    if( isFS(bc->kind, BCON::kAutoSlip) )
    {
      double length = sqrt( bc->nx * bc->nx  +  bc->ny * bc->ny );

      if( length < 1.0e-12 ) SF( nd->flag, NODE::kNoMoment );
    }
  }

  // -----------------------------------------------------------------------------------------------
  // determine rotation matrix from normal direction

  for( int i=0; i<region->Getnp(); i++ )
  {
    NODE* nd = region->Getnode(i);
    BCON* bc = &nd->bc;


    // set rotation flag for experimental boundaries -----------------------------------------------

    if( isFS(bc->kind, BCON::kSlip) )
    {
      SF( nd->flag, NODE::kRotat );
      bc->rot = BCON::kSlipFlowRot;
    }

    else if( isFS(bc->kind, BCON::kInlet) )
    {
      SF( nd->flag, NODE::kRotat );
      SF( bc->kind, BCON::kFixV );

      bc->rot = BCON::kNormalFlowRot;
    }


    // determine rotation matrix for automatic boundary conditions ---------------------------------

    if( isFS(nd->flag, NODE::kNoMoment)
        &&  !isFS(bc->kind, BCON::kSetUV)
        &&  !isFS(bc->kind, BCON::kOutlet)
        &&  !isFS(bc->kind, BCON::kOpenBnd) )
    {
      SF( bc->kind, BCON::kAutoSlip | BCON::kFixU | BCON::kFixV );

      nd->v.U = nd->v.V = 0.0;
    }

    else if( isFS(bc->kind, BCON::kFixV)  &&  isFS(bc->kind, BCON::kAutoSlip) )
    {
      SF( nd->flag, NODE::kRotat );

      // if the node is also on the outlet boundary the flow direction is set
      // normal to the outlet boundary
      if( isFS(bc->kind, BCON::kOutlet) || isFS(bc->kind, BCON::kOpenBnd) )
      {
        bc->rot = BCON::kNormalFlowRot;
      }
      else
      {
        bc->rot = BCON::kTangentFlowRot;
      }
    }
  }

  char text[200];
  sprintf( text, "\n (MODEL::SetRotation)    %s\n",
                 "transformation matrix for boundary elements set up" );
  REPORT::rpt.Output( text, 4 );


#ifdef kDebug_2
{
  char  name[80];

  sprintf( name, "SetRotation_%02d.report", subdom->pid+1 );
  FILE* id = fopen( name, "w" );

  fprintf( id, "  Node:  ((rot[i,j], i=1,2), j=1,2)\n" );

  for( int i=0; i<region->Getnp(); i++ )
  {
    NODE* nd = region->Getnode(i);
    BCON* bc = &nd->bc;

    fprintf( id, "%6d   %7.4lf  %7.4lf",
                 nd->Getname(),
                 bc->Getrot(0,0), bc->Getrot(1,0) );
    fprintf( id, "     %7.4lf  %7.4lf\n",
                 bc->Getrot(0,1), bc->Getrot(1,1) );
  }

  fclose( id );
}
#endif
}
