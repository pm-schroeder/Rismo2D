// ======================================================================================
//                            BCONSET | GetBcon | InitBcon
// ======================================================================================
//
// Copyright (C) 1992-2010  by  P.M. SCHROEDER
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
//
// methods:              description
// -------------------   ----------------------------------------------------------------
// BCONSET::GetBcon()
// BCONSET::InitBcon()   initialize specified boundary conditions
//
// ======================================================================================
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.09.1992     sc     Tirst implementation.
// 02.01.2010     ko     Implementation of diffuse Q at Nodes (Source/Sink).
// 31.03.2010     sc     Rismo 4.01.00: The implementation of boundary conditions has
//                       fundamentally changed. The two methods BCONSET::SetInflow()
//                       and BCONSET::SetOutflow() have vanished and are replaced by a
//                       call to the method BCONSET::InitBcon(), which should most
//                       likely produce identical results.
//                       For further information on changes see also the header file
//                       Bcon.h and the input method TIMEINT::Input_40100().
// 13.12.2012     sc     Rismo 4.03.00: BCONSET::stat[] array introduced (see TIMEINT)
//
// ======================================================================================

#include "Defs.h"
#include "Report.h"
#include "Node.h"
#include "Elem.h"
#include "Grid.h"
#include "Model.h"
#include "Memory.h"
#include "Time.h"

#include "Bcon.h"
#include "Project.h"

//#define kDebug_1

//////////////////////////////////////////////////////////////////////////////////////////

BCONSET::BCONSET()
{
  for( int i=0; i<kMaxCycles; i++ )
  {
    cycle[i]  = 0;
    stat[i]   = 0;
    turb[i]   = 0;
    disp[i]   = 0;
    cycit[i]  = 0;
    solver[i] = 0;
  }

  nofLines    = 0;
  nofNodes    = 0;

  nbc         = 0;
  nbcbuf      = 0;
  bc          = NULL;
  bcval       = NULL;

  ngct = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////

BCONSET::~BCONSET()
{
}

//////////////////////////////////////////////////////////////////////////////////////////

BCON* BCONSET::GetBcon( int no )
{
  for( int i=0; i<nbc; i++ )
  {
    if( bc[i].no == no )  return &bc[i];
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////

void BCONSET::InitBcon( PROJECT* project, TIME* actualTime, double* preQ )
{
  char text[200];

  MODEL* model = project->M2D;
  GRID*  ct    = model->control;
  GRID*  rg    = model->region;

  // class members:
  // nofLines - number of boundary-conditions at lines
  // nofNodes - number of boundary-conditions at nodes

  REPORT::rpt.Output( "\n (BCONSET::InitBcon)     initializing boundary conditions\n", 4 );


  // -------------------------------------------------------------------------------------
  // set up boundary conditions for lines and nodes
  // -------------------------------------------------------------------------------------

  // count number of nodes with boundary conditions to allocate memory -------------------

  // preparing task: initialize the marker for nodes
  for( int i=0; i<rg->Getnp(); i++ )
  {
    NODE* nd = rg->Getnode(i);
    nd->mark = false;
  }

  // loop over all boundary condition lines of actuall set ($TM_BOUND_LINE)
  for( int i=0; i<nofLines; i++ )
  {
    // loop over all (1D-)elements on control lines
    for( int c=0; c<ct->Getne(); c++ )
    {
      ELEM* el = ct->Getelem(c); // element c on the control line

      // (Control)lines are specified by their MAT-ID
      // To get only lines with a boundary condition
      // we are looking for identical number
      // in MAT-ID (el->type) and $TM_BOUND_LINE (bcLine[i].no)

      if( el->type == bcLine[i].no )
      {
        int nnd = el->Getnnd();
        // mark nodes of the element with boundary condition
        for( int j=0; j<nnd; j++ ) el->nd[j]->mark = true;
      }
    }
  }

  // loop on all nodes
  // count nodes with nd->mark == TRUE
  // nbc is the MAXIMUM number of nodes with BCs out of control-lines
  nbc = 0;

  for( int i=0; i<rg->Getnp(); i++ )
  {
    NODE* nd = rg->Getnode(i);
    if( nd->mark ) nbc++;
  }

  // add all nodes with BCs out of nodes
  // nbc - is now the MAXIMUM number of nodes with boundary conditions that can occur
  nbc += nofNodes;

  // allocate memory for boundary conditions
  if( nbc > nbcbuf )
  {
    if( nbcbuf )
    {
      delete[] bc;
      delete[] bcval;
    }

    nbcbuf = nbc;

    bc    = new BCON [nbcbuf];
    bcval = new BCVAL[nbcbuf];

    if( !bc || !bcval )
      REPORT::rpt.Error( kMemoryFault, "can not allocate memory - BCONSET::InitBcon() #1" );

    // set pointers to boundary values class BCVAL
    for( int i=0; i<nbcbuf; i++ )  bc[i].val = &bcval[i];
  }

  // -----------------------------------------------------------------------------------------------
  // set up boundary conditions for lines
  nbc = 0;     // reset nbc to zero

  for( int i=0; i<nofLines; i++ )
  {
    // transfer boundary conditions from lines to nodes
    // nbc counts the number of nodes with boundary in bc
    nbc = bcLine[i].GenBcon( model, actualTime, nbc, bc, preQ );
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // communicate generated boundary conditions for lines
  // !!! to check if this is necessary !!! (sc, 19.07.2005)
  //
  //# ifdef _MPI_
  //  if( subdom->npr > 1 )
  //  {
  //    INFACE* inface = subdom->inface;
  //
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
  //          ...
  //        }
  //
  //        MPI_Sendrecv( inface[s].sia1, npinf, MPI_CHAR, s, 1,
  //                      inface[s].ria1, npinf, MPI_CHAR, s, 1,
  //                      MPI_COMM_WORLD, &status );
  //
  //        for( int n=0; n<npinf; n++ )
  //        {
  //          NODE* nd = inface[s].node[n];
  //          ...
  //        }
  //      }
  //    }
  //  }
  //# endif
  ////////////////////////////////////////////////////////////////////////////////////////


  // -------------------------------------------------------------------------------------
  // append boundary conditions for nodes
  // -------------------------------------------------------------------------------------

  for( int i=0; i<nofNodes; i++ )
  {
    // The following loop over the current list of boundary conditions is searching
    // for boundary conditions at nodes which already have been set by boundary lines.

    int b;

    for( b=0; b<nbc; b++ )
    {
      if( bcNode[i].no == bc[b].no )  break;
    }

    // no prior bc found: copy specified BCs from bcnode to bc[] and increase k
    if( b == nbc )
    {
      // NOTE !!!
      // In the case of parallel computation (MPI) the specified node number of the
      // boundary condition "bcNode[i].no" is a global number and may be not localized
      // in this subdomain. "nd->Getno()" is the local index of the node.

      NODE* nd = NULL;

      if( project->subdom.npr > 1 )  nd = project->subdom.node[bcNode[i].no];
      else                           nd = rg->Getnode(bcNode[i].no);

      if( nd )
      {
        bc[b].no         = nd->Getno();
        bc[b].kind       = bcNode[i].kind;

        bc[b].val->U     = bcNode[i].val->U;
        bc[b].val->V     = bcNode[i].val->V;
        bc[b].val->Q     = bcNode[i].val->Q;
        bc[b].val->S     = bcNode[i].val->S;
        bc[b].val->K     = bcNode[i].val->K;
        bc[b].val->D     = bcNode[i].val->D;
        bc[b].val->C[0]  = bcNode[i].val->C[0];
        bc[b].val->Qb[0] = bcNode[i].val->Qb[0];
        bc[b].val->gct   = bcNode[i].val->gct;

        nbc++;
      }
    }

    // prior bc found: overwrite previously specified boundary conditions
    else
    {
      if( bc[b].Consistent( bcNode[i].kind ) )
      {
        if(     isFS(bcNode[i].kind, BCON::kInlet)
            ||  isFS(bcNode[i].kind, BCON::kSlip)
            ||  isFS(bcNode[i].kind, BCON::kSetUV) )
        {
          bc[b].kind  |= bcNode[i].kind;
          bc[b].val->U = bcNode[i].val->U;
          bc[b].val->V = bcNode[i].val->V;
        }

        if(     isFS(bcNode[i].kind, BCON::kOutlet)
            ||  isFS(bcNode[i].kind, BCON::kSetS) )
        {
          bc[b].kind   |= bcNode[i].kind;
          bc[b].val->S  = bcNode[i].val->S;
        }

        if(isFS(bcNode[i].kind, BCON::kSource) )
        {
          bc[b].kind  |= bcNode[i].kind;
          bc[b].val->Q = bcNode[i].val->Q;
        }

        if( isFS(bcNode[i].kind, BCON::kSetKD) )
        {
          bc[b].kind  |= bcNode[i].kind;
          bc[b].val->K = bcNode[i].val->K;
          bc[b].val->D = bcNode[i].val->D;
        }

        if( isFS(bcNode[i].kind, BCON::kSetC) )
        {
          bc[b].kind     |= bcNode[i].kind;
          bc[b].val->C[0] = bcNode[i].val->C[0];
        }

        if( isFS(bcNode[i].kind, BCON::kRateC) )
        {
          bc[b].kind      |= bcNode[i].kind;
          bc[b].val->Qb[0] = bcNode[i].val->Qb[0];
        }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  // remove boundary conditions from previous time step and initialize some values

  for( int i=0; i<rg->Getnp(); i++ )
  {
    NODE* nd = rg->Getnode(i);

    nd->cfw     = 0;
    nd->bc.kind = 0;
    nd->bc.val  = NULL;
  }

  // -------------------------------------------------------------------------------------
  // determine the area of node patches and store it into the array "patchArea"

  double* patchArea = (double*) MEMORY::memo.Array_nd( rg->Getnp() );

  for( int i=0; i<rg->Getnp(); i++ )  patchArea[i] = 0.0;

  for( int e=0; e<rg->Getne(); e++ )
  {
    ELEM* el  = rg->Getelem(e);
    int   nnd = el->Getnnd();

    if( !isFS(el->flag, ELEM::kDry) )
    {
      for( int j=0; j<nnd; j++ )
      {
        patchArea[el->nd[j]->Getno()] += el->area();
      }
    }
  }

# ifdef _MPI_
  project->subdom.Mpi_assemble( patchArea );
# endif

  // -------------------------------------------------------------------------------------
  // copy boundary conditions to nodes
  // -------------------------------------------------------------------------------------

  for( int i=0; i<nbc; i++ )
  {
    NODE* nd = rg->Getnode( bc[i].no );

    if( i == 0 )
    {
      sprintf( text, "\n (BCONSET::InitBcon)     %s kind %d at node %d %s\n",
                     "boundary condition", bc[i].kind, nd->Getname(), "set" );
      REPORT::rpt.Output( text, 5 );
    }
    else
    {
      sprintf( text, "                         %s kind %d at node %d %s\n",
                     "boundary condition", bc[i].kind, nd->Getname(), "set" );
      REPORT::rpt.Output( text, 5 );
    }

    // -----------------------------------------------------------------------------------

    if( isFS(bc[i].kind, BCON::kInlet) )
    {
      // avoid errors: remove kOutlet, kSlip, kFixVelo -----------------------------------
      CF( bc[i].kind, BCON::kOutlet | BCON::kSlip | BCON::kSetUV );
    }

    else if( isFS(bc[i].kind, BCON::kSlip) )
    {
      // avoid errors: remove kSetUV -----------------------------------------------------
      CF( bc[i].kind, BCON::kSetUV );

      double U = bc[i].val->U;
      double V = bc[i].val->V;

      if( fabs(U*U + V*V) > 1.0e-10 )
      {
        double Us = sqrt( U*U + V*V );

        bc[i].val->U /= Us;
        bc[i].val->V /= Us;

        bc[i].rot = BCON::kSlipFlowRot;

        SF( bc[i].kind, BCON::kFixV );      // delete Y-momentum-equation
      }

      else
      {
        sprintf( text, "\n (BCONSET::InitBcon)     %s %d\n",
                       "no flow direction specified at node", nd->Getname() );
        REPORT::rpt.Warning( kUserFault, text );
      }
    }

    if( isFS(bc[i].kind, BCON::kSource) )
    {
      bc[i].val->A = patchArea[nd->Getno()];
    }

    nd->bc.no   = nd->Getno();
    nd->bc.kind = bc[i].kind;
    nd->bc.val  = bc[i].val;


    // -----------------------------------------------------------------------------------
    // initialize flow parameters at nodes and boundaries ...

    if( isFS(bc[i].kind, BCON::kSetUV) )
    {
      nd->v.U = bc[i].val->U;
      nd->v.V = bc[i].val->V;

      SF( nd->bc.kind, BCON::kFixU | BCON::kFixV );
    }

    if( isFS(bc[i].kind, BCON::kSetS) )
    {
      nd->v.S = bc[i].val->S;
    }

    if( isFS(bc[i].kind, BCON::kLoglaw) )
    {
      double U = nd->v.U;
      double V = nd->v.V;
      double Us = sqrt( U*U + V*V );

      double dw = bc[i].val->dw;        // wall distance
      double kw = bc[i].val->kw;        // roughness height

      nd->cfw = Loglaw( Us, dw, kw, project->kappa, project->vk, project->g );
    }

    if( isFS(bc[i].kind, BCON::kSetKD) )
    {
      nd->v.K = bc[i].val->K;
      nd->v.D = bc[i].val->D;

      SF( nd->bc.kind, BCON::kFixK | BCON::kFixD );
    }

    if( isFS(bc[i].kind, BCON::kSetC) )
    {
      nd->v.C  = bc[i].val->C[0];
    }

    if( isFS(bc[i].kind, BCON::kRateC) )
    {
      nd->v.Qb = bc[i].val->Qb[0];
    }
  }

  // -------------------------------------------------------------------------------------
  // loop on all gauges to adapt the natural outflow boundary condition

  // set the targeted water elevation at gauges in PROJECT::gaugeSo[]
  // and copy the current water elevation at gauges to PROJECT::gaugeS[]

  for( int i=0; i<project->ngauge; i++ )
  {
    int gauge = 0;

    for( int j=0; j<nbc; j++ )
    {
      if( isFS(bc[j].kind, BCON::kOutlet)  &&  bc[j].val->gct.nocg > 0 )
      {
        gauge = bc[j].val->gct.nocg;

        if( project->gauge[i] == gauge )
        {
          project->gaugeSo[i] = bc[j].val->gct.So;
          break;
        }
      }
    }

#   ifdef kDebug_1
    REPORT::rpt.Message( 5, "%-25s%s%d%s%.6lf%s%.6lf\n",
                         " (BCONSET::InitBcon)", "Mark #1: gauge = ", gauge,
                         " | S[i] = ", project->gaugeS[i], " | So[i] = ", project->gaugeSo[i] );
#   endif

    ////////////////////////////////////////////////////////////////////////////////////////
#   ifdef _MPI_

    if( project->subdom.npr > 1 )
    {
      if( project->subdom.pid == 0 )
      {
        MPI_Status status;

        int    no;
        double So = 0.0;

        for( int s=1; s<project->subdom.npr; s++ )
        {
          MPI_Recv( &no, 1, MPI_INT,    s, 1, MPI_COMM_WORLD, &status );
          MPI_Recv( &So, 1, MPI_DOUBLE, s, 2, MPI_COMM_WORLD, &status );

#         ifdef kDebug_1
          REPORT::rpt.Message( 5, "%-25s%s%d%s%d%s%.6lf\n",
                                  " (BCONSET::InitBcon)",
                                  "Mark #1.1 (", s, "): Received no = ", no, " | So = ", So );
#         endif

          if( no > 0 )
          {
            project->gaugeSo[i] = So;
          }
        }
      }
      else
      {
        int    no;
        double So = project->gaugeSo[i];

        if( gauge == project->gauge[i] ) no = gauge;
        else                             no = 0;

        MPI_Send( &no, 1, MPI_INT,    0, 1, MPI_COMM_WORLD );
        MPI_Send( &So, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD );

#       ifdef kDebug_1
        REPORT::rpt.Message( 5, "%-25s%s%d%s%.6lf\n",
                                " (BCONSET::InitBcon)",
                                "Mark #1.2: Sent no = ", no, " | So = ", So );
#       endif
      }

      MPI_Barrier( MPI_COMM_WORLD );

//    MPI_Bcast( &project->gaugeSo[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

      if( project->subdom.pid == 0 )
      {
        for( int s=1; s<project->subdom.npr; s++ )
        {
          MPI_Send( &project->gaugeSo[i], 1, MPI_DOUBLE, s, 1, MPI_COMM_WORLD );
        }
      }
      else
      {
        MPI_Status status;
        MPI_Recv( &project->gaugeSo[i], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status );
      }
    }

#   endif
    ////////////////////////////////////////////////////////////////////////////////////////

#   ifdef kDebug_1
    REPORT::rpt.Message( 5, "%-25s%s%d\n",
                            " (BCONSET::InitBcon)", "Mark #2: gauge[i] = ", project->gauge[i] );
#   endif

    NODE *ndg = NULL;
    if( project->subdom.npr > 1 ) ndg = project->subdom.node[project->gauge[i] - 1];
    else                          ndg = rg->Getnode(project->gauge[i] - 1);

    if( ndg )
    {
      project->gaugeS[i] = ndg->v.S;
    }

    ////////////////////////////////////////////////////////////////////////////////////////
#   ifdef _MPI_

    if( project->subdom.npr > 1 )
    {
      if( project->subdom.pid == 0 )
      {
        MPI_Status status;

        int    no;
        double S;

        for( int s=1; s<project->subdom.npr; s++ )
        {
          MPI_Recv( &no, 1, MPI_INT,    s, 1, MPI_COMM_WORLD, &status );
          MPI_Recv( &S,  1, MPI_DOUBLE, s, 2, MPI_COMM_WORLD, &status );

#         ifdef kDebug_1
          REPORT::rpt.Message( 5, "%-25s%s%d%s%d%s%.6lf\n",
                                  " (BCONSET::InitBcon)",
                                  "Mark #2.1 (", s, "): Received no = ", no, " | S = ", S );
#         endif

          if( no > 0 )
          {
            project->gaugeS[i] = S;
          }
        }
      }
      else
      {
        int no;
        if( ndg ) no = project->gauge[i];
        else      no = 0;

        double S = project->gaugeS[i];

        MPI_Send( &no, 1, MPI_INT,    0, 1, MPI_COMM_WORLD );
        MPI_Send( &S,  1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD );

#       ifdef kDebug_1
        REPORT::rpt.Message( 5, "%-25s%s%d%s%.6lf\n",
                                " (BCONSET::InitBcon)",
                                "Mark #2.2: Sent no = ", no, " | S = ", S );
#       endif
      }

      MPI_Barrier( MPI_COMM_WORLD );

//    MPI_Bcast( &project->gaugeS[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

      if( project->subdom.pid == 0 )
      {
        for( int s=1; s<project->subdom.npr; s++ )
        {
          MPI_Send( &project->gaugeS[i], 1, MPI_DOUBLE, s, 1, MPI_COMM_WORLD );
        }
      }
      else
      {
        MPI_Status status;
        MPI_Recv( &project->gaugeS[i], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status );
      }
    }

#   endif
    ////////////////////////////////////////////////////////////////////////////////////////
  }

  // determine NODE::bc.val->gct.dS, which is the iterated difference to the specified
  // water level on outlet boundaries; the previous iteration step is project->gct[j].dS
  if( project->ngauge > 0 )
  {
    // allocate memory for the list of gauge controlled nodes
    if( ngct > 0 ) delete[] gct;

    ngct = 0;       // reset ngct to count the current number of controlled nodes

    for( int i=0; i<nbc; i++ )
    {
      NODE *nd = rg->Getnode( bc[i].no );
      if( isFS(nd->bc.kind, BCON::kOutlet)  &&  nd->bc.val->gct.nocg > 0 ) ngct++;
    }

    if( ngct > 0 )
    {
      REPORT::rpt.Message( 5, "\n\n%-25s%s\n",
                              " (BCONSET::InitBcon)", "adaption of outlet boundary conditions..." );
    }

#   ifdef kDebug_1
    REPORT::rpt.Message( 5, "%-25s%s%d\n",
                            " (BCONSET::InitBcon)", "Mark #3: ngct = ", ngct );
#   endif

    gct  = new BCVAL::GAUGECT[ngct];
    ngct = 0;       // reset ngct again

    int first = true;

    for( int i=0; i<nbc; i++ )
    {
      NODE *nd = rg->Getnode( bc[i].no );

      if( isFS(nd->bc.kind, BCON::kOutlet)  &&  nd->bc.val->gct.nocg > 0 )
      {
        gct[ngct]    = nd->bc.val->gct;
        gct[ngct].dS = 0.0;

        for( int j=0; j<project->ngct; j++ )
        {
          if( project->gct[j].node == gct[ngct].node )
          {
            gct[ngct].dS = project->gct[j].dS;
          }
        }

        double dS = 0.0;

        for( int j=0; j<project->ngauge; j++ )
        {
          if( project->gauge[j] == gct[ngct].nocg )
          {
            double S  = project->gaugeS[j];
            double So = project->gaugeSo[j];
            dS = So - S;
            break;
          }
        }

        gct[ngct].dS  += dS;
        nd->bc.val->S += gct[ngct].dS;
        ngct++;

        if( first )
        {
          first = false;
          REPORT::rpt.Message( 5, "%-25s%s%8d%s%.5lf\n",
                                  " ", "node = ", nd->Getname(), "  |  S = ", nd->bc.val->S );
        }
        else
        {
          REPORT::rpt.Message( 6, "%-25s%s%8d%s%.5lf\n",
                                  " ", "node = ", nd->Getname(), "  |  S = ", nd->bc.val->S );
        }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  // detach the temporary used memory

  MEMORY::memo.Detach( patchArea );
}

//////////////////////////////////////////////////////////////////////////////////////////
// compute friction coefficient for wall roughness (logarithmic law of the wall)
//////////////////////////////////////////////////////////////////////////////////////////

double BCONSET::Loglaw( double Us,           // scalar velocity
                        double dw,           // wall distance
                        double kw,           // roughness height
                        double ka,           // von Karman's constant ( = 0.41 )
                        double vk,           // kinematic viscosity
                        double g )           // gravity acceleration ( = 9.81 )
{
  int    iter;
  double Re, Ust, dwPlus;
  double cf = 100.0;

  if( Us < 1.0e-9 )
  {
    cf = 0.10;
  }

  else
  {
    double oldUst;

    if( kw < 1.0e-9 )
    {
      // iterative computation of friction velocity Ust

      iter = 0;
      Ust  = 100.0 * vk / dw;

      do
      {
        dwPlus = dw * Ust / vk;
        if( dwPlus < 11.0 ) dwPlus = 11.0;

        oldUst = Ust;
        Ust    = ka * Us / log(9.0 * dwPlus);

        iter++;

      } while( fabs((Ust-oldUst)/Ust) > 1.0e-6  &&  iter < 50 );
    }

    else
    {
      Ust = ka * Us / ( log(dw/kw) + 8.5 );
      Re  = kw * Ust / vk;

      if( Re < 70.0 )
      {
        // iterative computation of friction velocity Ust

        iter = 0;

        do
        {
          dwPlus = dw * Ust / vk;
          if( dwPlus < 11.0 ) dwPlus = 11.0;

          Re     = kw * Ust / vk;
          oldUst = Ust;

          if( Re < 3.32 )
          {
            Ust = ka * Us / log(9.0 * dwPlus);
          }
          else
          {
            Ust = ka * Us / ( log(dw/kw)  +  3.32*log(Re)/Re  +  ka * (8.5 - 9.96/Re) );
          }

          iter++;

        } while( fabs((Ust-oldUst)/Ust) > 1.0e-6  &&  iter < 50 );
      }
    }


    dwPlus = dw * Ust / vk;

    if( dwPlus < 11.0 )
    {
      Ust = 11.0 * vk / dw;             // constant
      Ust = sqrt (Us * vk / dw);        // linear
    }

    cf = Ust * Ust / Us / Us;
  }

  return cf;
}
