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
#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Sed.h"
#include "Project.h"


// ---------------------------------------------------------------------------------------
// method to initialize sediment parameters as qbe, Ls, sx, sy, ...
// ---------------------------------------------------------------------------------------

void SED::Initialize( PROJECT* project )
{
  isinit = true;

  MODEL*   model     = project->M2D;
  GRID*    rg        = model->region;
  int      rgnp      = rg->Getnp();

  // allocate memory for qbc, PLs, flow direction sx,sy and gradients --------------------

  qbc = (double*) MEMORY::memo.Array_nd( rgnp );
  for( int n=0; n<rgnp; n++ )  qbc[n]  = 0.0;

  if( lsType > 1 )
  {
    PLs = (double*) MEMORY::memo.Array_nd( rgnp );
    for( int n=0; n<rgnp; n++ )  PLs[n] = 1.0;
  }

  sx   = (double*) MEMORY::memo.Array_nd( rgnp );
  sy   = (double*) MEMORY::memo.Array_nd( rgnp );

  dzds = (double*) MEMORY::memo.Array_nd( rgnp );
  dzmx = (double*) MEMORY::memo.Array_nd( rgnp );
  dhds = (double*) MEMORY::memo.Array_nd( rgnp );

  // compute transport direction ---------------------------------------------------------
  Direction( project, model, sx, sy );

  // determine slope of bed in flow direction --------------------------------------------
  Slope( project, model, dzds, NULL, dzmx, dhds, NULL, NULL );

  // compute bed load transport capacity -------------------------------------------------
  Capacity( project, model, qbc, PLs, dzds, dhds );
}


void SED::Detach()
{
  // detach memory -----------------------------------------------------------------------

  if( qbc )   MEMORY::memo.Detach( qbc );
  if( PLs )   MEMORY::memo.Detach( PLs );
  if( sx )    MEMORY::memo.Detach( sx );
  if( sy )    MEMORY::memo.Detach( sy );
  if( dzds )  MEMORY::memo.Detach( dzds );
  if( dzmx )  MEMORY::memo.Detach( dzmx );
  if( dhds )  MEMORY::memo.Detach( dhds );

  qbc  = NULL;
  PLs  = NULL;
  sx   = NULL;
  sy   = NULL;
  dzds = NULL;
  dzmx = NULL;
  dhds = NULL;

  isinit = false;
}


// ---------------------------------------------------------------------------------------
// return available erosion of bed at node nd
// ---------------------------------------------------------------------------------------

double SED::Erode( NODE* nd )
{
  int no = nd->Getno();

  if( isFS(slope,kSLOPE_Max)  &&  dzmx[no]>maxSlope )
    return 0.0;

  else
    return nd->zor - nd->zero;
}


// ---------------------------------------------------------------------------------------
// determine bed slope at nodes
// ---------------------------------------------------------------------------------------

void SED::Slope( PROJECT* project,
                 MODEL*   model,
                 double*  dzds,      // slope of bed in flow direction (+/-)
                 double*  dzdn,      // slope of bed normal to flow direction (averaged)
                 double*  dzmx,      // maximum slope regarding only scours
                 double*  dhds,      // absolute change of flow depth
                 double*  dzdx,      // fe-averaged slope in x-direction
                 double*  dzdy )     // fe-averaged and y-direction
{
  GRID* rg = model->region;
  int   np = rg->Getnp();
  int   ne = rg->Getne();

  if( dzds )  for( int n=0; n<np; n++ )  dzds[n] = 0.0;
  if( dzdn )  for( int n=0; n<np; n++ )  dzdn[n] = 0.0;
  if( dzmx )  for( int n=0; n<np; n++ )  dzmx[n] = 0.0;
  if( dhds )  for( int n=0; n<np; n++ )  dhds[n] = 0.0;
  if( dzdx )  for( int n=0; n<np; n++ )  dzdx[n] = 0.0;
  if( dzdy )  for( int n=0; n<np; n++ )  dzdy[n] = 0.0;

  double* wds  = (double*) MEMORY::memo.Array_nd( np );
  double* wdn  = (double*) MEMORY::memo.Array_nd( np );
  double* wdh  = (double*) MEMORY::memo.Array_nd( np );
  double* wdxy = (double*) MEMORY::memo.Array_nd( np );

  for( int n=0; n<np; n++ )
  {
    wds[n]  = 0.0;
    wdn[n]  = 0.0;
    wdh[n]  = 0.0;
    wdxy[n] = 0.0;
  }

  // -------------------------------------------------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    ELEM* elem = rg->Getelem(e);

    if( isFS(elem->flag, ELEM::kBound) )  continue;

    SHAPE* shape = elem->GetLShape();

    int ncn = elem->Getncn();
    int nnd = elem->Getnnd();


    // -----------------------------------------------------------------------------------
    // local coordinates

    double x[kMaxNodes2D], y[kMaxNodes2D];

    x[0] = elem->nd[0]->x;
    y[0] = elem->nd[0]->y;

    for( int n=1; n<nnd; n++ )
    {
      x[n] = elem->nd[n]->x - x[0];
      y[n] = elem->nd[n]->y - y[0];
    }

    x[0] = y[0] = 0.0;


    // -----------------------------------------------------------------------------------
    // loop on nodes

    for( int n=0; n<nnd; n++ )
    {
      int no = elem->nd[n]->Getno();

      // ---------------------------------------------------------------------------------
      // transformation from local shape to global coordinates

      // form JACOBIAN transformation matrix with shape functions ------------------------

      double  trafo[2][2];

      double* dfdxPtr = shape->dndx[n];
      double* dfdyPtr = shape->dndy[n];

      double detj = shape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );

      double weight = elem->area();


      // ---------------------------------------------------------------------------------
      // compute values of linear shape functions at node i

      double dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

      for( int i=0; i<ncn; i++ )
      {
        dmdx[i] = trafo[0][0] * dfdxPtr[i]  +  trafo[0][1] * dfdyPtr[i];
        dmdy[i] = trafo[1][0] * dfdxPtr[i]  +  trafo[1][1] * dfdyPtr[i];
      }


      // ---------------------------------------------------------------------------------
      // compute bed slope with shape function at node "n"

      double dzx = 0.0;
      double dzy = 0.0;

      for( int i=0; i<ncn; i++ )
      {
        dzx += dmdx[i] * elem->nd[i]->zor;
        dzy += dmdy[i] * elem->nd[i]->zor;
      }

      // ---------------------------------------------------------------------------------
      // bed slope in flow direction

      if( dzds )
      {
        double U  = elem->nd[n]->v.U;
        double V  = elem->nd[n]->v.V;
        double Us = sqrt( U*U + V*V );

        //if( Us < 1.0e-6 )
        //{
        //  dzds[no] = 0.0;
        //}
        //else
        {
          //double tx = U / Us;
          //double ty = V / Us;
          double tx  = sx[no];  // approximate flow direction near the bed
          double ty  = sy[no];

          double xc, yc;
          elem->center( &xc, &yc );

          double dx = xc - elem->nd[n]->x;
          double dy = yc - elem->nd[n]->y;
          double Lk = sqrt( dx*dx + dy*dy );

          double cosa = dx*tx/Lk + dy*ty/Lk;

          if( cosa > 0.0 )
          {
            dzds[no] += weight * cosa * (dzx*tx + dzy*ty);
            wds[no]  += weight * cosa;
          }
        }
      }

      // ---------------------------------------------------------------------------------
      // bed slope normal to flow direction

      if( dzdn )
      {
        double U  = elem->nd[n]->v.U;
        double V  = elem->nd[n]->v.V;
        double Us = sqrt( U*U + V*V );

        if( Us < 1.0e-6 )
        {
          dzdn[no] = 0.0;
        }
        else
        {
          double nx = -V / Us;
          double ny =  U / Us;

          dzdn[no] += weight * (dzx*nx + dzy*ny);
          wdn[no]  += weight;
        }
      }

      // ---------------------------------------------------------------------------------

      if( dzmx )
      {
        double U  = elem->nd[n]->v.U;
        double V  = elem->nd[n]->v.V;
        double Us = sqrt( U*U + V*V );

        //if( Us < 1.0e-6 )
        //{
        //  dzmx[no] = 0.0;
        //}
        //else
        {
          //double tx = U / Us;
          //double ty = V / Us;
          double tx  = sx[no];  // approximate flow direction near the bed
          double ty  = sy[no];

          double xc, yc;
          elem->center( &xc, &yc );

          double dx = xc - elem->nd[n]->x;
          double dy = yc - elem->nd[n]->y;
          double Lk = sqrt( dx*dx + dy*dy );

          double cosa = dx*tx/Lk + dy*ty/Lk;

          double d = dzx*tx + dzy*ty;

          if( d * cosa > 0.0  &&  fabs(d) > dzmx[no] )
          {
            dzmx[no] = fabs(d);
          }
        }
      }

      // ---------------------------------------------------------------------------------

      if( dhds )
      {
        double dhx = 0.0;
        double dhy = 0.0;

        for( int i=0; i<ncn; i++ )
        {
          double H = elem->nd[i]->v.S - elem->nd[i]->zor;

          dhx += dmdx[i] * H;
          dhy += dmdy[i] * H;
        }

        double U  = elem->nd[n]->v.U;
        double V  = elem->nd[n]->v.V;
        double Us = sqrt( U*U + V*V );

        if( Us < 1.0e-6 )
        {
          dhds[no] = 0.0;
        }
        else
        {
          double tx = U / Us;
          double ty = V / Us;

          dhds[no] += weight * (dhx*tx + dhy*ty);
          wdh[no]  += weight;
        }
      }

      // ---------------------------------------------------------------------------------

      if( dzdx && dzdy )
      {
        dzdx[no] += weight * dzx;
        dzdy[no] += weight * dzy;
        wdxy[no] += weight;
      }
    }
  }


  if( dzds )
  {
    //////////////////////////////////////////////////////////////////////////////////////
#   ifdef _MPI_
    project->subdom.Mpi_assemble( dzds );
    project->subdom.Mpi_assemble( wds );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    for( int n=0; n<np; n++ )
    {
      if( wds[n] > 0.0 )  dzds[n] /= wds[n];
    }
  }


  if( dzdn )
  {
    //////////////////////////////////////////////////////////////////////////////////////
#   ifdef _MPI_
    project->subdom.Mpi_assemble( dzdn );
    project->subdom.Mpi_assemble( wdn );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    for( int n=0; n<np; n++ )
    {
      if( wdn[n] > 0.0 )  dzdn[n] /= wdn[n];
    }
  }


  if( dhds )
  {
    //////////////////////////////////////////////////////////////////////////////////////
#   ifdef _MPI_
    project->subdom.Mpi_assemble( dhds );
    project->subdom.Mpi_assemble( wdh );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    for( int n=0; n<np; n++ )
    {
      if( wdh[n] > 0.0 )  dhds[n] /= wdh[n];
    }
  }


  if( dzdx && dzdy )
  {
    //////////////////////////////////////////////////////////////////////////////////////
#   ifdef _MPI_
    project->subdom.Mpi_assemble( dzdx );
    project->subdom.Mpi_assemble( dzdy );
    project->subdom.Mpi_assemble( wdxy );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    for( int n=0; n<np; n++ )
    {
      if( wdxy[n] > 0.0 )
      {
        dzdx[n] /= wdxy[n];
        dzdy[n] /= wdxy[n];
      }
    }
  }


  // -------------------------------------------------------------------------------------

  MEMORY::memo.Detach( wds );
  MEMORY::memo.Detach( wdn );
  MEMORY::memo.Detach( wdh );
  MEMORY::memo.Detach( wdxy );
}


// ---------------------------------------------------------------------------------------
// determine transport direction
// ---------------------------------------------------------------------------------------

void SED::Direction( PROJECT* project, MODEL* model, double* sx, double* sy )
{
  GRID* rg = model->region;
  int   np = rg->Getnp();
  int   ne = rg->Getne();

  for( int n=0; n<np; n++ )
  {
    NODE* nd = rg->Getnode(n);
    int   no = nd->Getno();

    double U  = nd->v.U;
    double V  = nd->v.V;
    double Us = sqrt( U*U + V*V );

    if( Us > 1.0e-6 )
    {
      if( project->actualDisp > 0 )     // change flow direction due to secondary flow
      {
        U -= nd->Vsec * V/Us;
        V += nd->Vsec * U/Us;
        Us = sqrt( U*U + V*V );
      }

      sx[no] = U / Us;
      sy[no] = V / Us;
    }
    else
    {
      sx[no] = 0.0;
      sy[no] = 0.0;
    }
  }


  // -------------------------------------------------------------------------------------
  // compute change of flow direction due to gravity effects
  // an extended formulation of SEKINE & PARKER (1992) is used

  if( isFS(slope, kSLOPE_Angle) )
  {
    // initialization --------------------------------------------------------------------

    double grav  = project->g;
    double kappa = project->kappa;
    double vk    = project->vk;
    double rho   = project->rho;

    double rrg   = grav * (rhob/rho - 1.0);
    double rrg50 = rrg * d50;
    double Dst   = d50 * pow( rrg/vk/vk, 0.3333 );

    double* tand = (double*) MEMORY::memo.Array_nd( np );
    double* lmas = (double*) MEMORY::memo.Array_nd( np );

    for( int n=0; n<np; n++ )  tand[n] = lmas[n] = 0.0;

    // loop on elements ------------------------------------------------------------------

    for( int e=0; e<ne; e++ )
    {
      ELEM* elem = rg->Getelem(e);

      SHAPE* shape = elem->GetQShape();

      int ngp = shape->ngp;             // number of GAUSS points
      int nnd = shape->nnd;             // number of nodes

      // compute coordinates relative to first node --------------------------------------

      double x[kMaxNodes2D], y[kMaxNodes2D];

      x[0] = elem->nd[0]->x;
      y[0] = elem->nd[0]->y;

      for( int i=1; i<nnd; i++ )
      {
        x[i] = elem->nd[i]->x - *x;
        y[i] = elem->nd[i]->y - *y;
      }
      x[0] = y[0] = 0.0;

      // GAUSS point integration ---------------------------------------------------------

      for( int g=0; g<ngp; g++ )
      {
        // form JACOBIAN transformation matrix -------------------------------------------

        double  trafo[2][2];

        double* dfdxPtr = shape->dfdx[g];
        double* dfdyPtr = shape->dfdy[g];

        double  detj    = shape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );
        double  weight  = detj * shape->weight[g];

        // compute values of shape functions at GP g -------------------------------------

        double* n = shape->f[g];

        double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

        dfdxPtr = shape->dfdx[g];
        dfdyPtr = shape->dfdy[g];

        for( int i=0; i<nnd; i++ )
        {
          dndx[i] = trafo[0][0] * dfdxPtr[i] + trafo[0][1] * dfdyPtr[i];
          dndy[i] = trafo[1][0] * dfdxPtr[i] + trafo[1][1] * dfdyPtr[i];
        }

        // compute flow parameters and their derivatives at GP g -------------------------

        double dZdx = 0.0;
        double dZdy = 0.0;
        double U    = 0.0;
        double V    = 0.0;
        double H    = 0.0;
        double bx   = 0.0;
        double by   = 0.0;

        for( int i=0; i<nnd; i++ )
        {
          NODE* node = elem->nd[i];
          int   no   = node->Getno();

          double ndZ = node->zor;
          double ndU = node->v.U;
          double ndV = node->v.V;
          double ndH = node->v.S - ndZ;

          dZdx += dndx[i] * ndZ;
          dZdy += dndy[i] * ndZ;

          U    +=    n[i] * ndU;
          V    +=    n[i] * ndV;
          H    +=    n[i] * ndH;

          bx   +=    n[i] * sx[no];
          by   +=    n[i] * sy[no];
        }


        // -------------------------------------------------------------------------------

        double Us = sqrt( U*U + V*V );

        if( Us > 1.0e-6 )
        {
          double dzds    = bx * dZdx  +  by * dZdy;
          double dzdn    = by * dZdx  -  bx * dZdy;

          double sqrt_cf = kappa / log( 12.0 * H / d90 );

          double Utau    = sqrt_cf * Us;                    // friction velocity
          double theta   = Utau * Utau / rrg50;             // Shields-Parameter

          double thetacr = Shields_crit( Dst, dzds );
          double fs      = gammaSlope * pow( thetacr/theta, deltaSlope );

          // -----------------------------------------------------------------------------

          double wfdzdn  = weight * fs * dzdn;

          for( int i=0; i<nnd; i++ )
          {
            NODE* nd = elem->nd[i];
            int   no = nd->Getno();

            if( !isFS(nd->bc.kind,BCON::kFixU)  &&  !isFS(nd->bc.kind,BCON::kFixV) )
            {
            //tand[no] += n[i] * wfdzdn;
              tand[no] += wfdzdn;
            }

          //lmas[no] += n[i] * weight;
            lmas[no] += weight;
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
#   ifdef _MPI_
    project->subdom.Mpi_assemble( tand );
    project->subdom.Mpi_assemble( lmas );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    for( int n=0; n<np; n++ )
    {
      NODE* nd = rg->Getnode(n);
      int   no = nd->Getno();

      if( !isFS(nd->bc.kind,BCON::kFixU)  &&  !isFS(nd->bc.kind,BCON::kFixV) )
      {
        double t    = tand[no] / lmas[no];
        double tt   = sqrt( 1.0 + t*t );
        double sind = t / tt;
        double cosd = 1.0 / tt;

        double sinb = sx[no];
        double cosb = sy[no];

        sx[no] = sinb*cosd - cosb*sind;      // = sin(b - d)
        sy[no] = cosb*cosd + sinb*sind;      // = cos(b - d)

        // normalization of flow direction
        double r = sqrt( sx[no]*sx[no] + sy[no]*sy[no] );

        if( r > 1.0e-6 )
        {
          sx[no] /= r;
          sy[no] /= r;
        }
        else
        {
          sx[no] = 0.0;
          sy[no] = 0.0;
        }
      }
    }

    MEMORY::memo.Detach( tand );
    MEMORY::memo.Detach( lmas );
  }
}


// ---------------------------------------------------------------------------------------
// function to determine equilibrium bed load
// ---------------------------------------------------------------------------------------

void SED::Capacity( PROJECT* project,
                    MODEL*   model,
                    double*  qbc,
                    double*  PLs,
                    double*  dzds,
                    double*  dhds )
{
  GRID* rg = model->region;
  int   np = rg->Getnp();
  int   ne = rg->Getne();

  // determine equilibrium transport -----------------------------------------------------
  for( int n=0; n<np; n++ )
  {
    NODE* nd = rg->Getnode(n);
    int   no = nd->Getno();

    // determine scalar velocity, flow direction and bed slope ---------------------------
    double U  = nd->v.U;
    double V  = nd->v.V;
    double H  = nd->v.S - nd->z;

    double Us = sqrt( U*U + V*V );

    qbc[no] = Equilib( Us, H, dzds[no], dhds[no], project );

    if( PLs )  PLs[no] = 1.0 / GetLs( Us, H, project );
  }

  char text[200];
  sprintf( text, "\n%-25s%s\n",
                 " (SED::GetBLequi)",
                 "equilibrium bedload transport rates determined" );
  REPORT::rpt.Output( text, 4 );
}


// ---------------------------------------------------------------------------------------
// Bcon :  determine equilibrium values for qb at nodes
//         where equilibrium boundary condition is specified
// ---------------------------------------------------------------------------------------

void SED::Bcon( PROJECT* project, MODEL* model )
{
  GRID* rg = model->region;
  int   np = rg->Getnp();

  for( int n=0; n<np; n++ )
  {
    NODE* nd = rg->Getnode(n);
    int   no = nd->Getno();

    BCON* bc = &nd->bc;

    if( isFS(bc->kind,BCON::kRateC) )
    {
      // compute equilibrium transport capacity ------------------------------------------

      double U = nd->v.U;
      double V = nd->v.V;
      double H = nd->v.S - nd->z;

      double Us = sqrt( U*U + V*V );

      if( Us < 1.0e-6 )
      {
        qbc[no] = 0.0;
      }
      else
      {
        qbc[no] = Equilib( Us, H, dzds[no], dhds[no], project );
      }


      // set boundary condition ----------------------------------------------------------

      if( bc->val->Qb[0] < 0.0 )
      {
        nd->v.Qb = fabs(bc->val->Qb[0]) * qbc[no];
      }
      else
      {
        if( bc->val->Qb[0] <= qbc[no] )
        {
          nd->v.Qb = bc->val->Qb[0];
        }
        else
        {
          nd->v.Qb = qbc[no];
        }
      }
    }
  }
}


// ---------------------------------------------------------------------------------------
// Equilib: get equilibrium value of bedload transport
// ---------------------------------------------------------------------------------------

double SED::Equilib( double   Us,
                     double   H,
                     double   dzds,
                     double   dhds,
                     PROJECT* project )
{
  double qbe = 0.0;

  switch( project->sed.loadeq )
   {
     case 1:
       qbe = Vanrijn( Us, H, dzds, dhds, project );

     case 2:
     case 3:
       qbe = Meyerpm( Us, H, dzds, dhds, project );
   }


   if( isFS(slope, kSLOPE_Qbe) )                  // gravitation effect by WANG (1998)
   {
     double f = 1.0 - alfaSlope * dzds;

     if( f >= 0.0 )  qbe *= f;
     else            qbe  = 0.0;
   }

   return qbe;
}


// ---------------------------------------------------------------------------------------
// Meyerpm: formula of Meyer-Peter and Mueller
// ---------------------------------------------------------------------------------------

double SED::Meyerpm( double   Us,
                     double   H,
                     double   dzds,
                     double   dhds,
                     PROJECT* project )
{
  double g    = project->g;
  double rho  = project->rho;
  double vk   = project->vk;

  double rr   = rhob/rho - 1.0;

  // -------------------------------------------------------------------------------------

  double Dst     = 0.0;       // particle diameter parameter
  double thetacr = 0.0;

  if( loadeq == 2 )
  {
    thetacr = 0.047;
  }
  else
  {
    Dst     = d50 * pow( rr*g/vk/vk, 0.3333 );
    thetacr = Shields_crit( Dst, dzds );
  }

  // -------------------------------------------------------------------------------------

  double Ubseff = GetUtau( Us, H, dhds, project );

  double theta  = Ubseff*Ubseff / rr / g / d50;

  double qbe = 0.0;

  if( theta > thetacr )
  {
    qbe = 8.0 * Ubseff * d50 * theta * pow( (1.0 - thetacr/theta), 1.5 );
  }

  return qbe;
}


// ---------------------------------------------------------------------------------------
// Vanrijn: van Rijn formula
// ---------------------------------------------------------------------------------------

double SED::Vanrijn( double   Us,
                     double   H,
                     double   dzds,
                     double   dhds,
                     PROJECT* project )
{
  double  Dst     = 0.0;      // particle diameter parameter
  double  thetacr = 0.0;
  double  tstage  = 0.0;      // transport stage parameter
  double  Ubscr   = 0.0;      // critical bed-shear velocity

  double  g    = project->g;
  double  rho  = project->rho;
  double  vk   = project->vk;

  double  rr   = rhob/rho - 1.0;

  Dst     = d50 * pow( rr*g/vk/vk, 0.3333 );
  thetacr = Shields_crit( Dst, dzds );
  Ubscr   = sqrt( thetacr * rr * g * d50 );


  // -------------------------------------------------------------------------------------

  double Ubseff = GetUtau( Us, H, dhds, project );

  double qbe = 0.0;

  if( Ubseff > Ubscr  &&  Ubscr > 1.0e-6 )
  {
    tstage = ( Ubseff*Ubseff - Ubscr*Ubscr ) / Ubscr / Ubscr;

    qbe = 0.053  * sqrt( rr*g ) * pow( d50, 1.5 )
                 * pow( tstage, 2.1 ) * pow( Dst, -0.3 );
  }

  return qbe;
}


// ---------------------------------------------------------------------------------------
// GetLs: compute value for nonequilibrium parameter Ls
// ---------------------------------------------------------------------------------------

double SED::GetLs( double   Us,
                   double   H,
                   PROJECT* project )
{
  double Ls = 0.0;

  switch( lsType )
  {
    case 2:
      Ls = minLs;
      break;

    case 3:
      Ls = factLs * d50;
      break;

    case 4:
      Ls = Ls_vanrijn( Us, H, project );
      break;

    case 5:
      Ls = Ls_phillips( Us, H, project );
      break;
  }

  if( Ls < minLs )  Ls = minLs;

  return Ls;
}


// ---------------------------------------------------------------------------------------
// Ls_vanrijn: compute Ls by formula of van Rijn
// ---------------------------------------------------------------------------------------

double SED::Ls_vanrijn( double   Us,
                        double   H,
                        PROJECT* project )
{
  double  Dst     = 0.0;      // sedimentological particle diameter
  double  thetacr = 0.0;
  double  tstage  = 0.0;      // transport stage parameter
  double  Ubscr   = 0.0;      // critical bed-shear velocity


  double g   = project->g;
  double rho = project->rho;
  double vk  = project->vk;
  double rr  = rhob/rho - 1.0;

  Dst     = d50 * pow( rr*g/vk/vk, 0.3333 );
  thetacr = Shields_crit( Dst, 0.0 );
  Ubscr   = sqrt( thetacr*rr*g*d50 );


  // -------------------------------------------------------------------------------------

  double Ubseff = GetUtau( Us, H, 0.0, project );

  double Ls = 0.0;

  if( Ubseff > Ubscr  &&  Ubscr > 1.0e-6 )
  {
    tstage = ( Ubseff*Ubseff - Ubscr*Ubscr ) / Ubscr / Ubscr;

    Ls = 3.0 * d50 * pow( Dst, 0.6 ) * pow( tstage, 0.9 );
  }

  return Ls;
}


// ---------------------------------------------------------------------------------------
// Ls_phillips: compute Ls by formula of Phillips and Sutherland
// ---------------------------------------------------------------------------------------

double SED::Ls_phillips( double   Us,
                         double   H,
                         PROJECT* project )
{
  double alfaLs   = project->sed.alfaLs;

  double theta   = 0.0;
  double thetacr = 0.0;
  double Dst     = 0.0;
  double Ls      = 0.0;

  double g   = project->g;
  double rho = project->rho;
  double vk  = project->vk;

  double rr  = rhob/rho - 1.0;

  Dst     = d50 * pow( rr*g/vk/vk, 0.3333 );
  thetacr = Shields_crit( Dst, 0.0 );


  // -------------------------------------------------------------------------------------

  double Ubseff = GetUtau( Us, H, 0.0, project );

  theta  = Ubseff * Ubseff / rr / g / d50;

  Ls = alfaLs * ( theta - thetacr ) * d50;

  return Ls;
}


// ---------------------------------------------------------------------------------------
//  Shields_crit :  get critical value for shields-parameter
// ---------------------------------------------------------------------------------------

double SED::Shields_crit( double Dst,
                          double dzds )
{
  double thetacr = 0.055;

  if(                     Dst <=   6.0 ) thetacr = 0.109 * pow( Dst, -0.50 );
  else if(   6.0 < Dst && Dst <=  10.0 ) thetacr = 0.14  * pow( Dst, -0.64 );
  else if(  10.0 < Dst && Dst <=  20.0 ) thetacr = 0.04  * pow( Dst, -0.10 );
  else if(  20.0 < Dst && Dst <= 150.0 ) thetacr = 0.013 * pow( Dst,  0.29 );
  else if( 150.0 < Dst )                 thetacr = 0.055;

  // -------------------------------------------------------------------------------------
  // compute reduction/increase of thetacr due to graviational effects (slope)
  //   - dzdx is slope of bed elevation (gradient in flow direction)
  //   - phir is the tangens of friction angle for bed material

  if( phir > 0.001 && isFS(slope,kSLOPE_Shields) )        // KOCH, 1980
  {
    thetacr *= 1.0 + dzds/phir;
    if( thetacr < 0.01 )  thetacr = 0.01;                 // limit
  }

  return thetacr;
}


// ---------------------------------------------------------------------------------------
//  GetUtau:  compute friction velocity for grain (ks = d90)
// ---------------------------------------------------------------------------------------

double SED::GetUtau( double   Us,
                     double   H,
                     double   dhds,
                     PROJECT* project )
{
  if( H < project->hmin  ||  Us <= 1.0e-6 )  return 0.0;

  double grav  = project->g;
  double kappa = project->kappa;

  // compute sqrt(cf), cf = friction coefficient
  double sqrt_cf = kappa / log( 12.0 * H / d90 ); // roughness height ks = d90
  double Utau    = sqrt_cf * Us;                  // Utau = sqrt(cf) * Us

  if( isFS(slope,kSLOPE_Tau)  &&  dhds > 0.0 )    // NAKAGAWA et al. (1980)
  {
    double Fro = Us*Us/grav/H;          // Froude number (of undisturbed flow?)
    double df  = 1.0 - betaSlope*dhds/Fro;

    if( df > 2.0 )  df = 2.0;

    if( df > 0.0 )  Utau *= sqrt( df );
    else            Utau  = 0.0;
  }

  return Utau;
}


// ---------------------------------------------------------------------------------------
// compute mass balance Vtot for element
// ---------------------------------------------------------------------------------------

double SED::Balance( ELEM*   elem,           // input:  pointer to element
                     int     shape,          //         linear / qudratic interpolation
                     double* etaQb,          //         reduction factor on rigid bed
                     double  V[kMaxNodes2D], // output: change of volume (per node)
                     double* A,              //         area of element
                     double* Vmax )          //         maximum available sediment volume
{                                            //         integral of (Z - Zo) over element
  // element shape parameter -------------------------------------------------------------
  int ncn = elem->Getncn();

  SHAPE* lShape = elem->GetLShape();

  // initializations ---------------------------------------------------------------------
  *A = *Vmax = 0.0;

  // local element coordinates -----------------------------------------------------------
  double x[kMaxNodes2D], y[kMaxNodes2D];

  for( int i=0; i<ncn; i++ )
  {
    x[i] = elem->nd[i]->x - elem->nd[0]->x;
    y[i] = elem->nd[i]->y - elem->nd[0]->y;
  }

  // integrate over element area: A and Vmax ---------------------------------------------
  for( int g=0; g<lShape->ngp; g++ )
  {
    double* dfdxPtr = lShape->dfdx[g];
    double* dfdyPtr = lShape->dfdy[g];

    double  trafo[2][2];

    double  detj   = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
    double  weight = detj * lShape->weight[g];

    double* m = lShape->f[g];
    double  maxdz = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      NODE* node = elem->nd[i];
      if( node->z > node->zero )  maxdz += m[i] * (node->z - node->zero);
    }

    *A    += weight;
    *Vmax += weight * maxdz;
  }

  return Balance( elem, shape, etaQb, V );
}


double SED::Balance( ELEM*   elem,           // input:  pointer to element
                     int     shape,          //         linear / qudratic interpolation
                     double* etaQb,          //         reduction factor on rigid bed
                     double  V[kMaxNodes2D] )// output: change of volume (per node)
{
  // element shape parameter -------------------------------------------------------------
  int ncn = elem->Getncn();
  int nnd = elem->Getnnd();

  SHAPE* lineShape;
  if( shape == kLinear )  lineShape = SHAPE::get( kLine, 2 );
  else                    lineShape = SHAPE::get( kLine, 3 );

  // initializations ---------------------------------------------------------------------
  for( int i=0; i<nnd; i++ )  V[i] = 0.0;

  // local element coordinates -----------------------------------------------------------
  double x[kMaxNodes2D], y[kMaxNodes2D];

  for( int i=0; i<ncn; i++ )
  {
    x[i] = elem->nd[i]->x - elem->nd[0]->x;
    y[i] = elem->nd[i]->y - elem->nd[0]->y;
  }

  // compute transport through element boundary ------------------------------------------
  for( int i=0; i<ncn; i++ )
  {
    int j = (i + 1) % ncn;

    NODE* ndi  = elem->Getnode(i);
    NODE* ndj  = elem->Getnode(j);

    double ui  = sx[ndi->Getno()];
    double vi  = sy[ndi->Getno()];
    double qbi = etaQb[ndi->Getno()] * ndi->v.Qb;

    double uj  = sx[ndj->Getno()];
    double vj  = sy[ndj->Getno()];
    double qbj = etaQb[ndj->Getno()] * ndj->v.Qb;

    if( shape == kLinear )
    {
      for( int g=0; g<lineShape->ngp; g++ )
      {
        double* n      = lineShape->f[g];
        double* dn     = lineShape->dfdx[g];
        double  weight = lineShape->weight[g];

        // note: length of element side is included by (nx,ny) ---------------------------
        double nx = -dn[0]*y[i] - dn[1]*y[j];
        double ny =  dn[0]*x[i] + dn[1]*x[j];

        double ugi  = n[0] * ui;
        double vgi  = n[0] * vi;
        double qbgi = n[0] * qbi;

        double ugj  = n[1] * uj;
        double vgj  = n[1] * vj;
        double qbgj = n[1] * qbj;

        V[i] += weight * qbgi * (ugi*nx + vgi*ny);
        V[j] += weight * qbgj * (ugj*nx + vgj*ny);
      }
    }
    else
    {
      int m = i + ncn;

      NODE* ndm = elem->Getnode(m);

      double um  = sx[ndm->Getno()];
      double vm  = sy[ndm->Getno()];
      double qbm = etaQb[ndm->Getno()] * ndm->v.Qb;

      for( int g=0; g<lineShape->ngp; g++ )
      {
        double* n      = lineShape->f[g];
        double* dn     = lineShape->dfdx[g];
        double  weight = lineShape->weight[g];

        // note: length of element side is included by (nx,ny) ---------------------------
        double nx = -dn[0]*y[i] - dn[1]*y[j] - dn[2]*y[m];
        double ny =  dn[0]*x[i] + dn[1]*x[j] + dn[2]*x[m];

        double ugi  = n[0] * ui;
        double vgi  = n[0] * vi;
        double qbgi = n[0] * qbi;

        double ugj  = n[1] * uj;
        double vgj  = n[1] * vj;
        double qbgj = n[1] * qbj;

        double ugm  = n[2] * um;
        double vgm  = n[2] * vm;
        double qbgm = n[2] * qbm;

        V[i] += weight * qbgi * (ugi*nx + vgi*ny);
        V[j] += weight * qbgj * (ugj*nx + vgj*ny);
        V[m] += weight * qbgm * (ugm*nx + vgm*ny);
      }
    }
  }

  // compute Vtot as sum of V[i] ---------------------------------------------------------
  double Vtot = 0.0;
  for( int i=0; i<nnd; i++ )  Vtot += V[i];

  return Vtot;
}
