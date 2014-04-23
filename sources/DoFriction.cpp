// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================
//
// ---------------------------------------------------------------------------------------
// Compute a friction coefficient for each element; the influence of bottom roughness
// and single roughness elements are superimposed.
// Roughness elements are computed with the formula of LINDNER/PASCHE, bottom roughness
// may be characterized by different roughness coefficients, like:
//
//     (1) kLogLaw   logarithmic law of the wall
//     (2) kFricLaw  equivalent sand roughness coefficient 'ks' in [m]
//     (3) kChezy    roughness coefficient 'C' of Chezy
//     (4) kManning  roughness coefficient Manning's 'n' (= 1/kst)
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Type.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Project.h"

#include "Model.h"


void MODEL::DoFriction( PROJECT* project )
{
  int    i, cnt;

  double cf, h, U, V, Vres;
  double Ust, dwPlus, dwMax, dwMin, dwAve;
  char   text[500];


  GRID* rg = region;
  GRID* bd = bound;

  int   np = rg->Getnp();
  int   ne = rg->Getne();
  int   nb = bd->Getne();


  // allocations and initializations -----------------------------------------------------

  int* counter  = (int*) MEMORY::memo.Array_nd( np );

  for( i=0; i<np; i++ )
  {
    rg->Getnode(i)->cf = 0.0;
    counter[i] = 0;
  }


  // loop on all elements: compute friction coefficient at nodes -------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* el   = rg->Getelem(e);

    int   nnd  = el->Getnnd();
    TYPE* type = TYPE::Getid( el->type );

    if( type->rtype > 0 )
    {
      for( i=0; i<nnd; i++ )
      {
        h = el->nd[i]->v.S - el->nd[i]->z;

        if( h < project->hmin )  h = project->hmin;

        U = el->nd[i]->v.U;
        V = el->nd[i]->v.V;

        Vres = sqrt( U*U + V*V );

        // ### test - 10.01.2008 #########################################################
        // ### compute laminar roughness coefficient for marsh nodes
        //if( isFS(el->nd[i]->flag, NODE::kMarsh)  &&  Vres > 1.0e-6 )
        //{
        //  double Re = Vres * 4.0 * h / project->vk;
        //  cf = 8.0 / Re;
        //}
        //else
        //{
        //  cf = type->bottom( Vres, h, project->kappa,
        //                     project->vk, project->dw, project->g );
        //}
        // ### end of test - 10.01.2008 ##################################################

        switch( type->kslaw )
        {
          case 0:
            cf = type->bottom( Vres, h, project->kappa, project->vk, project->g,
                               project->rho, project->sed.rhob, type->d50, type->d90 );
            break;

          case 1:
            cf = type->bottom( Vres, h, project->kappa, project->vk, project->g,
                               project->rho, project->sed.rhob, type->d50, type->d90 );
            break;

          case 2:
            cf = type->bottom( Vres, h, project->kappa, project->vk, project->g,
                               project->rho, project->sed.rhob, project->sed.d50, project->sed.d90 );
            break;
        }

        el->nd[i]->cf += el->areaFact * cf;
        counter[ el->nd[i]->Getno() ]++;
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // MPI: assemble friction coefficient cf across interfaces

# ifdef _MPI_

  if( project->subdom.npr > 1 )
  {
    SUBDOM* subdom = &project->subdom;
    INFACE* inface = subdom->inface;

    // loop on all interfaces: exchange vector data --------------------------------------

    for( int s=0; s<subdom->npr; s++ )
    {
      MPI_Status status;

      int npinf = inface[s].np;

      if( npinf > 0 )
      {
        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];

          inface[s].sia2[n] = counter[nd->Getno()];
          inface[s].send[n] = nd->cf;
        }

        MPI_Sendrecv( inface[s].sia2, npinf, MPI_INT, s, 1,
                      inface[s].ria2, npinf, MPI_INT, s, 1,
                      MPI_COMM_WORLD, &status );

        MPI_Sendrecv( inface[s].send, npinf, MPI_DOUBLE, s, 1,
                      inface[s].recv, npinf, MPI_DOUBLE, s, 1,
                      MPI_COMM_WORLD, &status );

        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];

          counter[nd->Getno()] += inface[s].ria2[n];
          nd->cf               += inface[s].recv[n];
        }
      }
    }
  }

  MPI_Barrier( MPI_COMM_WORLD );
# endif
  ////////////////////////////////////////////////////////////////////////////////////////

  // divide cf-array by counter ----------------------------------------------------------

  for( int i=0; i<np; i++ )
  {
    NODE* nd = rg->Getnode(i);

    if( counter[i] ) nd->cf /= counter[i];

    counter[i] = 0;
  }

  REPORT::rpt.Output( "\n (MODEL::DoFriction)     friction coefficients determined\n", 4 );

  TYPE* type = TYPE::Getid( rg->Getelem(0)->type );

  type->statisLindner();


  MEMORY::memo.Detach( counter );
}
