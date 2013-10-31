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
// --------------------------------------------------------------------------------------
// The method Turbulence() computes the eddy viscosity depending on
// the value of turb = project->actualTurb:
//
//     turb == ##1: specified constant eddy viscosity
//                  vt = vto
//
//          == ##2: algebraic model for eddy viscosity
//                  vt = est*st * H * Utau
//
//          == ##3: algebraic shear stress model for eddy viscosity
//                  considering turbulence due to velocity gradients
//                  vt = cm * cd * Ke*Ke/De
//                  where: De = vt * Fi + PKv
//                         Ke = D * (c2D*D - c1D) / PDv
//
//          == ##4: Prandtl's mixing length model
//                  vt = lm * lm * sqrt( Fi )
//
//          == ##5: large eddy model
//                  vt = ls * ls * sqrt( Fi )
//
//          == ##6: Prandtl-Kolmogorov-Equation
//                  where K,D has been computed in a preceding cycle
//                       e.g.: 6x = Prandtl's mixing length model
//                             6x = large eddy model
//                             7x = KL-model
//                             8x = KD-model
//                  vt = cm * cd * K * K / D
//
//          == ##7: large eddy model + algebraic model for eddy viscosity
//                  vt = ls * ls * sqrt( Fi )  +  est*st * H * Utau
//
//     ----------------------------------------------------------------------------------
//     Further the following values of turb are recognized
//
//     turb == #1#: use anisotropic viscosity (ELDER-model)
//
//     ----------------------------------------------------------------------------------
//
//     turb == 1##: limit the minimum of eddy viscosity to the specified value
//          == 2##: do not change eddy viscosity during NR-iterations (stabilized)
//          == 3##: do both 1## + 2##
//
// ======================================================================================
//
#include "Defs.h"
#include "Report.h"
#include "Type.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Project.h"

#include "Grid.h"
#include "Model.h"

//#define kDebug_1


void GRID::Turbulence( PROJECT* project )
{
  char text[150];

  KDCONST* KD    = &project->KD;
  int      turb  = project->actualTurb;

  double   cm    = KD->cm;               // constants of the k-epsilon model
  double   cd    = KD->cd;
  double   c1D   = KD->c1D;
  double   c2D   = KD->c2D;


  ////////////////////////////////////////////////////////////////////////////////////////
  // ##1 use constant eddy viscosity

  if( isFS(turb, BCONSET::kVtConstant) )
  {
    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "using specified constant eddy viscosity" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
    double* vt  = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n]     = 0;
      vt[n]      = 0.0;
      node[n].uu = 0.0;
      node[n].uv = 0.0;
      node[n].vv = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          vt[nd[i]->Getno()] += type->vt;
          cnt[nd[i]->Getno()]++;
        }
      }
    }

    project->subdom.Mpi_assemble( cnt );
    project->subdom.Mpi_assemble( vt );


    // average values for vt -------------------------------------------------------------
    for( int n=0; n<np; n++ )
    {
      if( cnt[n] ) node[n].vt = vt[n] / cnt[n];
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( vt );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // ##7 large eddy model with algebraic bottom shear stresses

  else if( isFS(turb, BCONSET::kVtLES)  &&  isFS(turb, BCONSET::kVtAlgebraic) )
  {
    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "eddy viscosity with Smagorinsky's model and bottom shear stresses" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
    double* ste = (double*) MEMORY::memo.Array_nd( np );
    double* ls  = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n] = 0;
      ste[n] = 0.0;
      ls[n]  = 0.0;

      node[n].vt = 0.0;
      node[n].uu = 0.0;
      node[n].uv = 0.0;
      node[n].vv = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        double xmin = nd[0]->x;
        double xmax = nd[0]->x;
        double ymin = nd[0]->y;
        double ymax = nd[0]->y;

        for( int i=1; i<nnd; i++ )
        {
          if( nd[i]->x < xmin )  xmin = nd[i]->x;
          if( nd[i]->x > xmax )  xmax = nd[i]->x;

          if( nd[i]->y < ymin )  ymin = nd[i]->y;
          if( nd[i]->y > ymax )  ymax = nd[i]->y;
        }

        double dl = sqrt( (xmax-xmin)*(ymax-ymin) );

        for( int i=0; i<nnd; i++ )
        {
          int no = nd[i]->Getno();

          ls[no]  += type->lm * dl;
          ste[no] += type->st * type->KDest;
          cnt[no]++;
        }
      }
    }

    project->subdom.Mpi_assemble( cnt );
    project->subdom.Mpi_assemble( ls );
    project->subdom.Mpi_assemble( ste );

    double* phi = project->M2D->Phi2D();

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        // algebraic shear stresses
        ste[n] /= cnt[n];

        double U = node[n].v.U;
        double V = node[n].v.V;
        double H = node[n].v.S - node[n].z;

        double Utau = sqrt( node[n].cf * (U*U + V*V) );
        double cK   = 1.0 / sqrt( node[n].cf );
//        double PKv  = cK * Utau*Utau*Utau / H;

        // LES
        double L  = ls[n] / cnt[n];

        double vt = ste[n] * H * Utau  +  L * L * sqrt( phi[n] );

//        double D  = vt * phi[n];// + PKv;
//        double K  = sqrt( vt * D / cm / cd );

        double K = (L * L * phi[n]  +  Utau*Utau * sqrt(cK*ste[n]) ) / sqrt(cm/cd);

        node[n].vt = vt;
        node[n].uu =
        node[n].vv = -2.0/3.0 * K;
      }
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( ls );
    MEMORY::memo.Detach( ste );
    MEMORY::memo.Detach( phi );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // ##2 compute eddy viscosity with algebraic model

  else if( isFS(turb, BCONSET::kVtAlgebraic) )
  {
    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "eddy viscosity with algebraic model" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
    double* ste = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n] = 0;
      ste[n] = 0.0;

      node[n].vt = 0.0;
      node[n].uu = 0.0;
      node[n].uv = 0.0;
      node[n].vv = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          int no = nd[i]->Getno();

          ste[no] += type->st * type->KDest;
          cnt[no]++;
        }
      }
    }

    project->subdom.Mpi_assemble( cnt );
    project->subdom.Mpi_assemble( ste );

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        ste[n] /= cnt[n];

        double U = node[n].v.U;
        double V = node[n].v.V;
        double H = node[n].v.S - node[n].z;

        double Utau = sqrt( node[n].cf * (U*U + V*V) );
        double cK   = 1.0 / sqrt( node[n].cf );

        node[n].vt = ste[n] * H * Utau;

        node[n].uu =
        node[n].vv = -2.0/3.0 * Utau*Utau * sqrt(cK*ste[n]/cm/cd);
      }
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( ste );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // ##3 compute eddy viscosity with algebraic shear stress model

  else if( isFS(turb, BCONSET::kVtAlgShear) )
  {
    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "eddy viscosity with algebraic shear stress model #1#" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
    double* ste = (double*) MEMORY::memo.Array_nd( np );
    double* lm  = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n] = 0;
      ste[n] = 0.0;
      lm[n]  = 0.0;

      node[n].vt = 0.0;
      node[n].uu = 0.0;
      node[n].uv = 0.0;
      node[n].vv = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          int no = nd[i]->Getno();

          lm[no]  += type->lm;
          ste[no] += type->st * type->KDest;
          cnt[no]++;
        }
      }
    }

    project->subdom.Mpi_assemble( cnt );
    project->subdom.Mpi_assemble( lm );
    project->subdom.Mpi_assemble( ste );

    double* phi = project->M2D->Phi2D();

#   ifdef kDebug_1
    FILE *dbg = fopen( "Turbulence.dat", "w" );
    fprintf( dbg, "NODE\tcf\tPKv\tlm\tH\tUtau\tphi[n]\tNS\tx[0]\tx[1]\tx[2]\n" );
#   endif

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        double U  = node[n].v.U;
        double V  = node[n].v.V;
        double Us = sqrt( U*U + V*V );
        double H  = node[n].v.S - node[n].z;
        double cf = node[n].cf;

        double Utau = sqrt( cf ) * Us;
        double cK   = 1.0 / sqrt( cf );

        if( H < project->hmin ) H = project->hmin;

        ste[n] /= cnt[n];
        lm[n]  /= cnt[n];
        lm[n]  *= H;

        double PKv  = cK * Utau*Utau*Utau / H;

        double lm2 = lm[n] * lm[n];
        double lm4 = lm2 * lm2;

        double ac = 1.0;
        double bc = 0.0;
        double cc = -lm4 * phi[n];
        double dc = -lm4 * PKv;

        double x[3];
        int real = project->Cube( ac, bc, cc, dc, x );

        double vt = x[0];

        if( real == 3 )
        {
          if( x[1] > vt ) vt = x[1];
          if( x[2] > vt ) vt = x[2];
        }

        if( vt < 0.0 ) vt = ste[n] * H * Utau;

        double D = vt * phi[n]  +  PKv;
        double K = sqrt( vt * D / cm / cd );

#         ifdef kDebug_1
        if( real == 1 )
        {
          fprintf( dbg, "%d\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%d\t%12.4le\n",
                   node[n].Getname(), cf, PKv, lm[n], H, Utau, phi[n], real, x[0] );
        }
        else
        {
          fprintf( dbg, "%d\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%d\t%12.4le\t%12.4le\t%12.4le\n",
                   node[n].Getname(), cf, PKv, lm[n], H, Utau, phi[n], real, x[0], x[1], x[2] );
        }
#         endif

        node[n].vt = vt;

        node[n].uu =
        node[n].vv = -2.0/3.0 * K;
      }
    }

#   ifdef kDebug_1
    fclose( dbg );
#   endif

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( lm );
    MEMORY::memo.Detach( ste );
    MEMORY::memo.Detach( phi );
  }

//  else if( isFS(turb, BCONSET::kVtAlgShear) ) // EXPERIMENTAL #2# - not working
//  {                                           //                    in common cases
//    sprintf( text, "\n (GRID::Turbulence)      %s\n",
//                   "eddy viscosity with algebraic shear stress model #2#" );
//    REPORT::rpt.Output( text, 4 );
//
//    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
//    double* ste = (double*) MEMORY::memo.Array_nd( np );
//
//    for( int n=0; n<np; n++ )
//    {
//      cnt[n] = 0;
//      ste[n] = 0.0;
//
//      node[n].vt = 0.0;
//      node[n].uu = 0.0;
//      node[n].uv = 0.0;
//      node[n].vv = 0.0;
//    }
//
//    for( int e=0; e<ne; e++ )
//    {
//      ELEM* el = &elem[e];
//
//      if( !isFS(el->flag, ELEM::kDry) )
//      {
//        TYPE* type = TYPE::Getid( el->type );
//
//        NODE** nd  = el->nd;
//        int    nnd = el->Getnnd();
//
//        for( int i=0; i<nnd; i++ )
//        {
//          int no = nd[i]->Getno();
//          ste[no] += type->KDest * type->st;
//          cnt[no]++;
//        }
//      }
//    }
//
//    project->subdom.Mpi_assemble( cnt );
//    project->subdom.Mpi_assemble( ste );
//
//    double* phi = project->M2D->Phi2D();
//
//#   ifdef kDebug_1
//    FILE *dbg = fopen( "Turbulence.dat", "w" );
//    fprintf( dbg, "NODE\tcf\tPKv\tPDv\tH\tUtau\tphi[n]\tNS\tx[0]\tx[1]\tx[2]\n" );
//#   endif
//
//    for( int n=0; n<np; n++ )
//    {
//      if( cnt[n] )
//      {
//        ste[n] /= cnt[n];
//
//        double U  = node[n].v.U;
//        double V  = node[n].v.V;
//        double H  = node[n].v.S - node[n].z;
//        double cf = node[n].cf;
//
//        double Utau = sqrt( cf * (U*U + V*V) );
//        double cK   = 1.0 / sqrt( cf );
//
//        if( H < project->hmin  )
//        {
//          node[n].vt = ste[n] * H * Utau;
//
//          node[n].uu =
//          node[n].vv = -2.0/3.0 * Utau*Utau * sqrt(cK*ste[n]/cm/cd);
//        }
//        else
//        {
//          double PKv  = cK * Utau*Utau*Utau / H;
//
//          double cD   = c2D * sqrt( cm/ste[n] / sqrt(cf*cf*cf) );
//          double PDv  = cD * (Utau*Utau/H) * (Utau*Utau/H);
//
//          double alfa = PDv*PDv / phi[n] / cm / cd;
//          double beta = c1D * PKv;
//          double delc = c2D - c1D;
//
//          double ac = delc * delc;
//          double bc = 2.0 * delc * beta;
//          double cc = beta * beta - alfa;
//          double dc = alfa * PKv;
//
//          double x[3];
//          int real = project->Cube( ac, bc, cc, dc, x );
//
//#         ifdef kDebug_1
//          if( real == 1 )
//          {
//            fprintf( dbg, "%d\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%d\t%12.4le\n",
//                     node[n].Getname(), cf, PKv, PDv, H, Utau, phi[n], real, x[0] );
//          }
//          else
//          {
//            fprintf( dbg, "%d\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%12.4le\t%d\t%12.4le\t%12.4le\t%12.4le\n",
//                     node[n].Getname(), cf, PKv, PDv, H, Utau, phi[n], real, x[0], x[1], x[2] );
//          }
//#         endif
//
//          double K, D;
//
//          if( real == 1 )
//          {
//            D = x[0];
//          }
//          else
//          {
//            D = x[2];
//
//            if( D < PKv )
//            {
//              D = x[0];
//              if( D < PKv )  D = x[1];
//            }
//          }
//
//          if( D < PKv )
//          {
//            K = sqrt( cK * ste[n] / cm / cd ) * Utau*Utau;
//            D = cK * Utau*Utau*Utau / H;
//          }
//          else
//          {
//            //K = sqrt( D * (D - PKv) / Fi[n] / cm / cd );
//            K = (c2D*D - c1D*(D-PKv)) * D / PDv;
//          }
//
//          node[n].vt = cm * cd * K * K / D;
//
//          node[n].uu =
//          node[n].vv = -2.0/3.0 * K;
//        }
//      }
//    }
//
//#   ifdef kDebug_1
//    fclose( dbg );
//#   endif
//
//    MEMORY::memo.Detach( cnt );
//    MEMORY::memo.Detach( ste );
//    MEMORY::memo.Detach( phi );
//  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // Prandtl's mixing length model

  else if( isFS(turb, BCONSET::kVtMixingLength) )
  {
    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "eddy viscosity with Prandtl's mixing length model" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
    double* lm  = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n] = 0;
      lm[n]  = 0.0;

      node[n].vt = 0.0;
      node[n].uu = 0.0;
      node[n].uv = 0.0;
      node[n].vv = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          lm[nd[i]->Getno()] += type->lm;
          cnt[nd[i]->Getno()]++;
        }
      }
    }

    project->subdom.Mpi_assemble( cnt );
    project->subdom.Mpi_assemble( lm );

    double* phi = project->M2D->Phi2D();

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        double H = node[n].v.S - node[n].z;
        double L = H * lm[n] / cnt[n];

        node[n].vt = L * L * sqrt( phi[n] );
      }
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( lm );
    MEMORY::memo.Detach( phi );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // large eddy model

  else if( isFS(turb, BCONSET::kVtLES) )
  {
    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "eddy viscosity with Smagorinsky's model" );
    REPORT::rpt.Output( text, 4 );

    int*    cnt = (int*) MEMORY::memo.Array_nd( np );
    double* ls  = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n]     = 0;
      ls[n]      = 0.0;

      node[n].vt = 0.0;
      node[n].uu = 0.0;
      node[n].uv = 0.0;
      node[n].vv = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        NODE** nd  = el->nd;
        int    nnd = el->Getnnd();

        double xmin = nd[0]->x;
        double xmax = nd[0]->x;
        double ymin = nd[0]->y;
        double ymax = nd[0]->y;

        for( int i=1; i<nnd; i++ )
        {
          if( nd[i]->x < xmin )  xmin = nd[i]->x;
          if( nd[i]->x > xmax )  xmax = nd[i]->x;

          if( nd[i]->y < ymin )  ymin = nd[i]->y;
          if( nd[i]->y > ymax )  ymax = nd[i]->y;
        }

        double dl = sqrt( (xmax-xmin)*(ymax-ymin) );

        for( int i=0; i<nnd; i++ )
        {
          int no = nd[i]->Getno();

          ls[no] += type->lm * dl;
          cnt[no]++;
        }
      }
    }

    project->subdom.Mpi_assemble( cnt );
    project->subdom.Mpi_assemble( ls );

    double* phi = project->M2D->Phi2D();

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        double L  = ls[n] / cnt[n];

        double vt = L * L * sqrt( phi[n] );
        double D  = vt * phi[n];
        double K  = sqrt( vt * D / cm / cd );

        node[n].vt = vt;

        node[n].uu =
        node[n].vv = -2.0/3.0 * K;
      }
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( ls );
    MEMORY::memo.Detach( phi );
  }


  // -------------------------------------------------------------------------------------

  else if( isFS(turb, BCONSET::kVtPrandtlKol) )
  {
    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "eddy viscosity from K,D-values" );
    REPORT::rpt.Output( text, 4 );

    for( int n=0; n<np; n++ )
    {
      if( !isFS(node[n].flag, NODE::kDry) )
      {
        double K = node[n].v.K;
        double D = node[n].v.D;

        if( D < project->minD )  D = project->minD;
        node[n].vt = KD->cm * KD->cd * K * K / D;
      }
    }
  }


  // -------------------------------------------------------------------------------------

  else
  {
    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "initializing eddy viscosity with zero" );
    REPORT::rpt.Output( text, 4 );

    for( int n=0; n<np; n++ )  node[n].vt = 0.0;
  }


  // check eddy viscosity for minimum values ---------------------------------------------

  if( isFS(turb,BCONSET::kVtMin) )
  {
    double min   = -1.0;
    double minel = 0.0;
    int    jMin  = 0;

    for( int e=0; e<ne; e++ )
    {
      ELEM *el = &elem[e];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        TYPE* type = TYPE::Getid( el->type );

        for( int i=0; i<el->Getnnd(); i++ )
        {
          NODE *nd = el->Getnode(i);

          // check for minimum viscosity
          if( min < 0.0  ||  nd->vt < min )  min = nd->vt;

          if( nd->vt < type->vt )
          {
            // minVt = type->vt will be hold in EQS::RegionUVS2D(), ... 02-03-2010, SC
            // nd->vt = type->vt;
            minel = type->vt;
            jMin++;
          }
        }
      }
    }
//    for( int n=0; n<np; n++ )
//    {
//      if( !isFS(node[n].flag, NODE::kDry) )
//      {
//        // check for minimum viscosity
//        if( min < 0.0  ||  node[n].vt < min )  min = node[n].vt;

//        if( node[n].vt < minVt )
//        {
////        minVt will be hold only in EQS::RegionUVS2D(), ... 02-03-2010, SC
////        node[n].vt = minVt;
//          jMin++;
//        }
//      }
//    }

    if( jMin )
    {
      sprintf( text, "\n (GRID::Turbulence)      %s (%11.4le) set to %11.4le\n",
                     "minimum of eddy viscosity", min, minel );
      REPORT::rpt.Output( text, 4 );
    }
  }


  // check Peclet-number for maximum values ----------------------------------------------
/*
  double maxPe = 50.0;

  if( isFS(turb,BCONSET::kPeMax) )
  {
    double* Lx = (double*) MEMORY::memo.Array_nd( np );
    double* Ly = (double*) MEMORY::memo.Array_nd( np );

    int* cnt = (int*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      Lx[n]  = 0.0;
      Ly[n]  = 0.0;
      cnt[n] = 0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      if( isFS( el->flag, ELEM::kRegion ) )
      {
        int ncn = el->getncn();
        int nnd = el->getnnd();

        double xe = 0.0;
        double ye = 0.0;

        for( int i=1; i<ncn; i++ )
        {
          double x = fabs( el->nd[i]->x - el->nd[0]->x );
          double y = fabs( el->nd[i]->y - el->nd[0]->y );

          if( x > xe ) xe = x;
          if( y > ye ) ye = y;
        }

        for( int i=0; i<nnd; i++ )
        {
          int no = el->nd[i]->Getno();
          cnt[no]++;

          Lx[no] += xe;
          Ly[no] += ye;
        }
      }
    }

    for( int n=0; n<np; n++ )
    {
      if( cnt[n] > 0 )
      {
        Lx[n] /= cnt[n];
        Ly[n] /= cnt[n];
      }
    }

    MEMORY::memo.Detach( cnt );


    double max  = -1.0;
    int    kMax = 0;


    for( int n=0; n<np; n++ )
    {
      if( !isFS(node[n].flag, NODE::kDry) )
      {
        int no = node[n].Getno();

        double U = node[n].v.U;
        double V = node[n].v.V;

        if( node[n].vt > 0.0 )
        {
          double Pe = sqrt( U*Lx[no]*U*Lx[no] + V*Ly[no]*V*Ly[no] ) / node[n].vt;

          if( max < 0.0  ||  Pe > max )  max = Pe;

          if( Pe > maxPe )
          {
            node[n].vt *= Pe / maxPe;
            kMax++;
          }
        }
      }
    }

    MEMORY::memo.Detach( Lx );
    MEMORY::memo.Detach( Ly );

    if( kMax )
    {
      sprintf( text, "\n (GRID::Turbulence)      %s (%11.4le) set to %11.4le\n",
                     "maximum Peclet-number", max, maxPe );
      REPORT::rpt.Output( text, 4 );
    }
  }
*/

  // =====================================================================================
  //
  // compute eddy viscosity anisotrop - depending on flow direction (Elder)
  // see method: eddyDisp() for advection-diffusion-equation
  //
  // =====================================================================================

  if( isFS(turb, BCONSET::kVtAnisotrop) )
  {
    // initialization --------------------------------------------------------------------

    int*    cnt = (int*)    MEMORY::memo.Array_nd( np );
    double* exx = (double*) MEMORY::memo.Array_nd( np );
    double* eyy = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )
    {
      cnt[n] = 0;
      exx[n] = 0.0;
      eyy[n] = 0.0;
    }

    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      TYPE* type = TYPE::Getid( el->type );

      double exxInit = type->estx / type->KDest;
      double eyyInit = type->esty / type->KDest;

      NODE** nd  = el->nd;
      int    nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        exx[nd[i]->Getno()] += exxInit;
        eyy[nd[i]->Getno()] += eyyInit;

        cnt[nd[i]->Getno()]++;
      }
    }

    project->subdom.Mpi_assemble( cnt );
    project->subdom.Mpi_assemble( exx );
    project->subdom.Mpi_assemble( eyy );

    sprintf( text, "\n (GRID::Turbulence)      %s\n",
                   "anisotropic eddy viscosity (ELDER model)" );
    REPORT::rpt.Output( text, 4 );


    for( int n=0; n<np; n++ )
    {
      if( cnt[n] )
      {
        node[n].exx = exx[n] / cnt[n];
        node[n].eyy = eyy[n] / cnt[n];
      }

//    changed: multiplication with eddy viscosity is now performed in coefs()
//    29. April 2004, SC
//    node[n].Exx *= node[n].vt;
//    node[n].Eyy *= node[n].vt;
    }


    // transform local Exx,Exy,Eyy into global Exx,Exy,Eyy -------------------------------

    for( int n=0; n<np; n++ )
    {
      // angle between Ures and global x-Axis --------------------------------------------

      double U  = node[n].v.U;
      double V  = node[n].v.V;
      double Us = sqrt( U*U + V*V );

      double sa = 0.0;
      double ca = 1.0;

      if( Us > 1.0e-6 )
      {
        sa = V / Us;
        ca = U / Us;
      }


      // transformation ------------------------------------------------------------------

      double exx = node[n].exx;
      double eyy = node[n].eyy;

      node[n].exx =  exx * ca * ca  +  eyy * sa * sa;
      node[n].exy = (exx - eyy) * sa * ca;
      node[n].eyy =  exx * sa * sa  +  eyy * ca * ca;
    }

    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( exx );
    MEMORY::memo.Detach( eyy );
  }
}
