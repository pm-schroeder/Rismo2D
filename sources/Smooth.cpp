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
#include "Defs.h"
#include "Report.h"
#include "Memory.h"
#include "Bcon.h"
#include "Node.h"
#include "Elem.h"

#include "Grid.h"


void GRID::SmoothS( int passes )
{
  int*    counter = (int*)    MEMORY::memo.Array_nd( np );
  double* Save    = (double*) MEMORY::memo.Array_el( ne );
  double* UH      = (double*) MEMORY::memo.Array_nd( np );
  double* VH      = (double*) MEMORY::memo.Array_nd( np );


  // compute local discharge qx = U * H and qy = V * H

  for( int n=0; n<np; n++ )
  {
    if( !isFS(node[n].flag, NODE::kDry) )
    {
      double H = node[n].v.S - node[n].z;

      UH[n] = node[n].v.U * H;
      VH[n] = node[n].v.V * H;
    }
  }


  for( int p=0; p<passes; p++ )
  {
    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        int ncn = el->Getncn();

        double S = 0.0;

        for( int j=0; j<ncn; j++ )
        {
          S += el->nd[j]->v.S;
        }

        Save[el->Getno()] = S / ncn;
      }
    }

    for( int n=0; n<np; n++ )  counter[n] = 1;

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        int nnd = el->Getnnd();

        for( int j=0; j<nnd; j++ )
        {
          BCON* bc = &el->nd[j]->bc;

          if(     !isFS(bc->kind, BCON::kSetS)
              &&  !isFS(bc->kind, BCON::kOutlet) )
          {
            el->nd[j]->v.S += Save[el->Getno()];

            counter[el->nd[j]->Getno()]++;
          }
        }
      }
    }

    for( int n=0; n<np; n++ )  node[n].v.S /= counter[n];
  }


  // reestablish local discharge
  for( int n=0; n<np; n++ )
  {
    if( !isFS(node[n].flag, NODE::kDry) )
    {
      double H = node[n].v.S - node[n].z;

      node[n].v.U = UH[n] / H;
      node[n].v.V = VH[n] / H;
    }
  }

  MEMORY::memo.Detach( counter );
  MEMORY::memo.Detach( Save );
  MEMORY::memo.Detach( UH );
  MEMORY::memo.Detach( VH );
}


void GRID::SmoothKD( int passes )
{
  int*    Kcounter = (int*)    MEMORY::memo.Array_nd( np );
  int*    Dcounter = (int*)    MEMORY::memo.Array_nd( np );
  double* Kave     = (double*) MEMORY::memo.Array_el( ne );
  double* Dave     = (double*) MEMORY::memo.Array_el( ne );


  for( int p=0; p<passes; p++ )
  {
    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        int nnd = el->Getnnd();

        double K = 0.0;
        double D = 0.0;

        for( int j=0; j<nnd; j++ )
        {
          K += el->nd[j]->v.K;
          D += el->nd[j]->v.D;
        }

        Kave[el->Getno()] = K / nnd;
        Dave[el->Getno()] = D / nnd;
      }
    }


    for( int n=0; n<np; n++ )  Kcounter[n] = Dcounter[n] = 1;


    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry)  )
      {
        int nnd = el->Getnnd();

        for( int j=0; j<nnd; j++ )
        {
          BCON* bc = &el->nd[j]->bc;

          if( !isFS(bc->kind, BCON::kFixK) )
          {
            el->nd[j]->v.K += Kave[el->Getno()];

            Kcounter[el->nd[j]->Getno()]++;
          }

          if( !bc  ||  (bc && !isFS(bc->kind, BCON::kFixD)) )
          {
            el->nd[j]->v.D += Dave[el->Getno()];

            Dcounter[el->nd[j]->Getno()]++;
          }
        }
      }
    }

    for( int n=0; n<np; n++ )
    {
      node[n].v.K /= Kcounter[n];
      node[n].v.D /= Dcounter[n];
    }
  }

  MEMORY::memo.Detach( Kcounter );
  MEMORY::memo.Detach( Dcounter );
  MEMORY::memo.Detach( Kave );
  MEMORY::memo.Detach( Dave );
}
