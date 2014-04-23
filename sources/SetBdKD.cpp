// ======================================================================================
//
// Copyright (C) 1992-2010  by  P.M. SCHROEDER
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
// -----------------------------------------------------------------------------------------
// Function to determine values for the turbulent kinetic energy (K) and
// the dissipation (D) on boundaries.
//
//
//  kind            description                              boundary
//  ----  ----------------------------------  -----------  ------------
//   1    symmetrie plane; normal gradients     Neumann      outlet
//        (diffusive flow) set to zero
//
//   2    fixed K and D are computed with      Dirichlet     side, inlet
//        algebraic eddy viscosity model;
//        K and D related to flow depth
//        and bottom shear stresses;
//
//   3    diffusive flow of K set to zero       Neumann      side
//        fixed D computed from wall shear     Dirichlet
//        stress at solid boundaries
//        local equilibrium:
//        production equals dissipation
//
//   4    fixed K and D are computed from      Dirichlet     side
//        wall shear stress at solid
//        boundaries (local equilibrium:
//        production equals dissipation)
//
//   5    fixed K and D, leave K and D         Dirichlet     inlet
//        unchanged from initialization
//
// -----------------------------------------------------------------------------------------
// Michael Schroeder, November   1992
//                    March      1994
//                    September  1999
// -----------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Type.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Project.h"

#include "Model.h"


void MODEL::SetBoundKD( PROJECT* project )
{
  GRID* rg = region;


  // interpolate dimensionless diffusivity to nodes ----------------------------------------

  double* stEst   = (double*) MEMORY::memo.Array_nd( rg->Getnp() );
  int*    counter = (int*)    MEMORY::memo.Array_nd( rg->Getnp() );

  for( int i=0; i<rg->Getnp(); i++ )
  {
    counter[i] = 0;
    stEst[i]   = 0.0;
  }


  for( int e=0; e<rg->Getne(); e++ )
  {
    ELEM* el = rg->Getelem(e);

    TYPE* type = TYPE::Getid( el->type );

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      int no = el->nd[i]->Getno();

      stEst[no] += type->st * type->KDest;
      counter[no]++;
    }
  }


  // compute production of turbulence energy ---------------------------------------------

  double* phi = NULL;

  if( isFS(project->actualTurb, BCONSET::kVtAlgShear) )
  {
    phi = project->M2D->Phi2D();
  }


  // set boundary flags ------------------------------------------------------------------

  for( int i=0; i<rg->Getnp(); i++ )
  {
    if( counter[i] ) stEst[i] /= counter[i];

    BCON* bcon = &rg->Getnode(i)->bc;

    if( bcon  &&  isFS(bcon->kind, BCON::kAutoKD) )
      CF( bcon->kind, BCON::kAutoKD | BCON::kFixK | BCON::kFixD );
  }


  for( int e=0; e<bound->Getne(); e++ )
  {
    ELEM* el = bound->Getelem(e);

    int nnd = el->Getnnd();

    for( int i=0; i<nnd; i++ )
    {
      SF( el->nd[i]->bc.kind, BCON::kAutoKD );
    }
  }


  for( int i=0; i<rg->Getnp(); i++ )
  {
    NODE* nd = rg->Getnode(i);

    if( isFS(nd->flag, NODE::kDry) ) continue;

    BCON* bcon = &nd->bc;

    if( isFS(bcon->kind, BCON::kSetKD) )
    {
      // set specified boundary conditions for K and D
      SF( bcon->kind, BCON::kFixK | BCON::kFixD );
      nd->v.K = bcon->val->K;
      nd->v.D = bcon->val->D;
    }


    else if( isFS(bcon->kind, BCON::kAutoKD) )
    {
      // compute K and D according to the specified model:
      // 1. set inflow and outflow boundaries
      // 2. overwrite values at wall and surface boundaries
      int index = -1;

           if( isFS(bcon->kind, BCON::kInlet) )  index = kIndInflow;
      else if( isFS(bcon->kind, BCON::kOutlet) ) index = kIndOutflow;
      else                                       index = kIndSide;

      if( index >= 0  &&  phi )
      {
        SetNodeKD(  project, nd, project->KDBcon[index], index, phi[i], stEst[i] );
      }
      else
      {
        SetNodeKD(  project, nd, project->KDBcon[index], index, 0.0, stEst[i] );
      }
    }
  }


  REPORT::rpt.Output( "\n (MODEL::SetBoundKD)     KD boundary condition set up\n", 4 );

  MEMORY::memo.Detach( stEst );
  MEMORY::memo.Detach( counter );
  if( phi ) MEMORY::memo.Detach(phi );

  //SmoothKDBound( project->smoothPassesBC );
}


void MODEL::SetNodeKD( PROJECT* project,
                       NODE*    node,
                       int      KDType,
                       int      index,
                       double   Fi,
                       double   stEst )
{
  KDCONST* KD    = &project->KD;
  double   hmin  =  project->hmin;
  double   kappa =  project->kappa;
  double   vk    =  project->vk;
  int      turb  =  project->actualTurb;

  double cf, U, V, H, Utau;

  double cm  = KD->cm;             // constants of the k-epsilon model
  double cd  = KD->cd;
  double c1D = KD->c1D;
  double c2D = KD->c2D;

  BCON* bcon = &node->bc;

  switch( KDType )
  {
    // -------------------------------------------------------------------------------------
    // diffusive flow set to zero for K and D (Neumann): nothing to do here
    case 1:
      break;

    // -------------------------------------------------------------------------------------
    // apply the algebraic eddy viscosity model (Dirichlet)
    case 2:
      {
        SF( bcon->kind, BCON::kAutoKD | BCON::kFixK | BCON::kFixD );
        cf   = node->cf;
        U    = node->v.U;
        V    = node->v.V;
        Utau = sqrt( cf * (U*U + V*V) );

        H  = node->v.S - node->z;
        if( H <= hmin ) H = hmin;

        double cK = 1.0 / sqrt( cf );

        node->v.K = sqrt( cK * stEst / cm / cd ) * Utau*Utau;
        node->v.D = cK * Utau*Utau*Utau / H;
      }

      if( isFS(turb, BCONSET::kVtAlgebraic) )
      {
        double cK = 1.0 / sqrt( cf );

        node->v.K = sqrt( cK * stEst / cm / cd ) * Utau*Utau;
        node->v.D = cK * Utau*Utau*Utau / H;
      }
      else if( isFS(turb, BCONSET::kVtAlgShear) )
      {
        double cK   = 1.0 / sqrt( cf );
        double PKv  = cK * Utau*Utau*Utau / H;

        double cD   = c2D * sqrt( cm/stEst / sqrt(cf*cf*cf) );
        double PDv  = cD * (Utau*Utau/H) * (Utau*Utau/H);

        double alfa = PDv*PDv / Fi / cm / cd;
        double beta = c1D * PKv;
        double delc = c2D - c1D;

        double ac = delc * delc;
        double bc = 2.0 * delc * beta;
        double cc = beta * beta - alfa;
        double dc = alfa * PKv;

        double x[3];
        int real = project->Cube( ac, bc, cc, dc, x );

        if( real == 1 )
        {
          node->v.D = x[0];
        }
        else
        {
          node->v.D = x[2];

          if( node->v.D < PKv )
          {
            node->v.D = x[0];
            if( node->v.D < PKv )  node->v.D = x[1];
          }
        }

        if( node->v.D < PKv )
        {
          node->v.K = sqrt( cK * stEst / cm / cd ) * Utau*Utau;
          node->v.D = cK * Utau*Utau*Utau / H;
        }
        else
        {
          node->v.K = (c2D*node->v.D - c1D*(node->v.D-PKv)) * node->v.D / PDv;
        }
      }
      break;

    // -------------------------------------------------------------------------------------
    // diffusive flow set to zero for K (Neumann)
    // D fixed value computed from wall distance (Dirichlet)
    case 3:
      if( index == kIndSide  &&  node->cfw > 1.0e-9  &&  node->bc.val )
      {
        SF( bcon->kind, BCON::kAutoKD | BCON::kFixD );
        cf   = node->cfw;
        U    = node->v.U;
        V    = node->v.V;
        Utau = sqrt( cf * (U*U + V*V) );

        node->v.D = Utau * Utau * Utau / kappa / node->bc.val->dw;
      }
      break;

    // -------------------------------------------------------------------------------------
    // local equilibrium, wall boundaries (Dirichlet)
    case 4:
      if( index == kIndSide  &&  node->cfw > 1.0e-9  &&  node->bc.val )
      {
        SF( bcon->kind, BCON::kAutoKD | BCON::kFixK | BCON::kFixD );
        cf   = node->cfw;
        U    = node->v.U;
        V    = node->v.V;
        Utau = sqrt( cf * (U*U + V*V) );

        node->v.K = Utau * Utau / sqrt(cm);
        node->v.D = Utau * Utau * Utau / kappa / node->bc.val->dw;
      }
      break;

    // -------------------------------------------------------------------------------------
    // error detected
    default:
      REPORT::rpt.Error( "KD boundary model not implemented - setNodeKD(3)" );
  }
}


// -----------------------------------------------------------------------------------------
// smooth transition of KD distribution on inflow, outflow and surface

void MODEL::SmoothKDBound( int passes )
{
  GRID* rg = region;
  GRID* bd = bound;


  // allocate temporary used memory --------------------------------------------------------

  int*    counter = (int*)    MEMORY::memo.Array_nd( rg->Getnp() );
  double* Kave    = (double*) MEMORY::memo.Array_el( bd->Getne() );
  double* Dave    = (double*) MEMORY::memo.Array_el( bd->Getne() );

  for( int i=0; i<passes; i++ )
  {
    for( int e=0; e<bd->Getne(); e++ )
    {
      ELEM* el = bd->Getelem(e);

      if(     isFS(el->flag, ELEM::kInlet)
          ||  isFS(el->flag, ELEM::kOutlet) )
      {
        int nnd = el->Getnnd();

        double K = 0.0;
        double D = 0.0;

        for( int j=0; j<nnd; j++ )
        {
          K += el->nd[j]->v.K;
          D += el->nd[j]->v.D;
        }

        Kave[e] = K / nnd;
        Dave[e] = D / nnd;
      }
    }

    for( int n=0; n<rg->Getnp(); n++ )  counter[n] = 1;

    for( int e=0; e<bd->Getne(); e++ )
    {
      ELEM* el = bd->Getelem(e);

      if(     isFS(el->flag, ELEM::kInlet)
          ||  isFS(el->flag, ELEM::kOutlet) )
      {
        int nnd = el->Getnnd();

        for( int j=0; j<nnd; j++ )
        {
          NODE* nd = el->nd[j];
          BCON* bc = &nd->bc;

          if( isFS(bc->kind, BCON::kAutoKD) )
          {
            nd->v.K += Kave[e];
            nd->v.D += Dave[e];

            counter[nd->Getno()]++;
          }
        }
      }
    }


    for( int n=0; n<rg->Getnp(); n++ )
    {
      if( counter[n] )
      {
        rg->Getnode(n)->v.K /= counter[n];
        rg->Getnode(n)->v.D /= counter[n];
      }
    }
  }

  MEMORY::memo.Detach( Dave );
  MEMORY::memo.Detach( Kave );
  MEMORY::memo.Detach( counter );
}
