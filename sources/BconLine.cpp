// =======================================================================================
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
// =======================================================================================

#include "Defs.h"
#include "Report.h"
#include "Shape.h"
#include "Node.h"
#include "Elem.h"
#include "Grid.h"
#include "Model.h"
#include "Subdom.h"
#include "Time.h"

#include "Bcon.h"


BCONLINE::BCONLINE() : BCON()
{
  ndat = 0;
}


BCONLINE::~BCONLINE()
{
}


// ---------------------------------------------------------------------------------------
// generate boundary conditions for lines

int BCONLINE::GenBcon( MODEL  *model,
                       TIME   *actualTime,
                       int     b,
                       BCON   *bcon,
                       double *preQ )
{
  char text[120];

  GRID* rg = model->region;

  sprintf( text, "\n (BCONLINE::GenBcon)     %s %d\n",
                 "generating boundary conditions for line", no );
  REPORT::rpt.Output( text, 3 );


  // initialize node flags ---------------------------------------------------------------

  for( int i=0; i<rg->Getnp(); i++ )
  {
    rg->Getnode(i)->mark = false;
  }


  // -------------------------------------------------------------------------------------
  // complement inlet and outlet boundary conditions

  if( isFS(kind, BCON::kQTInlet) )  SF( kind, BCON::kQInlet );
  if( isFS(kind, BCON::kQInlet) )   SF( kind, BCON::kInlet );

  if( isFS(kind, BCON::kSQOutlet) ) SF( kind, BCON::kOutlet );
  if( isFS(kind, BCON::kSTOutlet) ) SF( kind, BCON::kOutlet );


  // -------------------------------------------------------------------------------------
  // Work on user specified boundary conditions along control lines. Boundary conditions
  // have been read from the timestep file (TIMEINT::Input_XXXXXX()). The user might
  // have specified boundary conditions over "count" lines per timestep.
  //
  // The following boundary conditions are possible:
  //
  // 1. inlet boundary conditions
  //
  //    a) BCON::kInlet,    count =1: One specific discharge q=(U*H,V*H), which
  //                                  is constant along the control line.
  //                        count >1: A specific discharge profile q(x). The value
  //                                  of q at nodes must be interpolated.
  //    b) BCON::kQInlet,   count =1: The integral discharge Q over the control line.
  //                                  A specific discharge profile q must be computed
  //                                  from flow depth and bottom roughness.
  //    c) BCON::kQTInlet,  count>=1: The discharge hydrograph Q(T). Very much like
  //                                  boundary condtion BCON::kQInlet.
  //
  // 2. outlet boundary conditions
  //
  //    a) BCON::kOutlet,   count =1: One water surface elevation S, which
  //                                  is constant along the control line.
  //                        count >1: A water surface elevation profile S(x). The value
  //                                  of S at nodes will be interpolated.
  //    b) BCON::kSQOutlet, count>=1: Water surface elevation depends on discharge S(Q).
  //                                  The Dicharge Q is computed for the control line.
  //    c) BCON::kSTOutlet, count>=1: Water surface elevation depends on time S(T).
  //
  // 3. further boundary conditions
  //    BCON::kOpenBnd,     count =1:
  //    BCON::kSlip,        count>=1: Compute values at nodes from specified value(x,y).
  //    BCON::kSetUV,       count>=1:
  //    BCON::kSetS,        count>=1:
  //    BCON::kSetKD,       count>=1:
  //    BCON::kSetC,        count>=1:
  //    BCON::kRateC,       count>=1:
  // -------------------------------------------------------------------------------------

  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kInlet
  ////////////////////////////////////////////////////////////////////////////////////////

  if( isFS(kind, BCON::kInlet) )
  {
    b = Inlet( model, actualTime, b, bcon );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kOutlet
  ////////////////////////////////////////////////////////////////////////////////////////

  else if( isFS(kind, BCON::kOutlet) )
  {
    b = Outlet( model, actualTime, b, bcon, preQ );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // further boundary conditions ...
  ////////////////////////////////////////////////////////////////////////////////////////

  long ki = kind;

  while( ki )
  {
    b = Further( model, actualTime, b, bcon, &ki );
  }


  ////////////////////////////////////////////////////////////////////////////////////////

  return b;
}


int BCONLINE::FindInterval( double  x,
                            double  y,
                            int     n,
                            double* xarray,
                            double* yarray,
                            double* l )
{
  double dx, dy, lx, ly;

  for( int i=0; i<n-1; i++ )
  {
    dx = xarray[i+1] - xarray[i];
    lx = x           - xarray[i];

    if( fabs(dx) > kEpsilon ) lx /= dx;
    else                      lx  = dx = 0.0;

    if( yarray )
    {
      dy = yarray[i+1] - yarray[i];
      ly = y           - yarray[i];

      if( fabs(dy) > kEpsilon ) ly /= dy;
      else                      ly  = dy = 0.0;
    }
    else
    {
      dy = 0.0;
      ly = 0.0;
    }

    if( fabs(dx) < kEpsilon  &&  fabs(lx) < kEpsilon )
    {
      if( ly > -0.001  &&  ly < 1.001 )
      {
        *l = ly;
        return i;
      }
    }

    else if( fabs(dy) < kEpsilon  &&  fabs(ly) < kEpsilon )
    {
      if( lx > -0.001  &&  lx < 1.001 )
      {
        *l = lx;
        return i;
      }
    }

    else
    {
      if( lx > -0.001  &&  lx < 1.001  &&  ly > -0.001  &&  ly < 1.001)
      {
        *l = ( lx + ly ) / 2.0;
        return i;
      }
    }
  }

  return -1;
}


BCON* BCONLINE::FindBcon( int no, int b, BCON* bcon )
{
  BCON* bc = NULL;

  for( int i=0; i<b; i++ )
  {
    if( bcon[i].no == no )
    {
      bc = &bcon[i];
      break;
    }
  }

  return bc;
}


//////////////////////////////////////////////////////////////////////////////////////////

int BCONLINE::Inlet( MODEL* model, TIME* actualTime, int b, BCON* bcon )
{
  int  firstMiss = true;
  char text[120];

  GRID* ct = model->control;
  GRID* rg = model->region;

  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kQInlet
  ////////////////////////////////////////////////////////////////////////////////////////

  if( isFS(kind, BCON::kQInlet) )
  {
    double Qspec  = U[0];
    double Sspec  = S[0];

    // determine discharge for actual time -----------------------------------------------
    if( isFS(kind, BCON::kQTInlet) )
    {
      double L;

      int j = FindInterval( actualTime->Getsec(), 0.0, ndat, x, NULL, &L );

      if( j >= 0 )
      {
        Qspec = U[j]  +  L * (U[j+1] - U[j]);
        Sspec = S[j]  +  L * (S[j+1] - S[j]);
      }
    }

    // compute a specific discharge profile ----------------------------------------------
    double q[3];
    NODE*  node[3];

    double Q = 0.0;

    for( int c=0; c<ct->Getne(); c++ )
    {
      ELEM* el = ct->Getelem(c);

      if( el->type != no )  continue;

      SHAPE* qShape = el->GetQShape();

      int nnd = qShape->nnd;

      node[0] = el->nd[0];                 // corner nodes
      node[1] = el->nd[1];
      node[2] = el->nd[2];                 // midside node

      // initialize the specific discharge profile q = pow(h,3/2)
      for( int i=0; i<nnd; i++ )
      {
        if( isFS(node[i]->flag, NODE::kDry) || isFS(node[i]->flag, NODE::kMarsh) )
        {
          q[i] = 0.0;
        }
        else
        {
          double h = Sspec - node[i]->z;
          if( h > 0.0 )  q[i] = sqrt( h*h*h );
          else           q[i] = 0.0;

          double cf = node[i]->cf;
          if( cf > 0.0 )  q[i] /= sqrt(cf);
        }

        node[i]->mark = true;
      }

      // integrate the specific discharge profile
      for( int g=0; g<qShape->ngp; g++ )       // loop on GAUSS points
      {
        double* N  = qShape->f[g];             // quadratic shape
        double* dN = qShape->dfdx[g];

        // compute normal vector at Gauss point g
        double dx = dN[0]*node[0]->x + dN[1]*node[1]->x + dN[2]*node[2]->x;
        double dy = dN[0]*node[0]->y + dN[1]*node[1]->y + dN[2]*node[2]->y;

        double qg = N[0]*q[0] + N[1]*q[1] + N[2]*q[2];

        double weight = qShape->weight[g] * sqrt(dx*dx + dy*dy);

        Q += weight * qg;
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: sum of discharge Q from all subdomains
#   ifdef _MPI_
    Q = model->subdom->Mpi_sum( Q );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

//  double Qratio = (Q > 0.0)? (-fabs(Qspec/Q)) : (0.0);
    double Qratio = -Qspec/Q;

    // set boundary conditions at nodes --------------------------------------------------

    for( int i=0; i<rg->Getnp(); i++ )
    {
      NODE* nd = rg->Getnode(i);

      // regard only marked nodes on the control line
      if( nd->mark )
      {
        double q = 0.0;

        if( isFS(nd->flag, NODE::kDry) || isFS(nd->flag, NODE::kMarsh) )
        {
          q = 0.0;
        }
        else
        {
          double h = Sspec - nd->z;
          if( h > 0.0 )  q = sqrt( h*h*h );
          else           q = 0.0;

          double cf = nd->cf;
          if( cf > 0.0 )  q /= sqrt(cf);
        }

        BCON* bc = FindBcon( nd->Getno(), b, bcon );

        if( !bc )
        {
          bcon[b].no     = nd->Getno();
          bcon[b].kind   = kind;

          bcon[b].val->U = q * Qratio;
          bcon[b].val->S = Sspec;

          b++;
        }
        else if( bc->Consistent(kind) )
        {
          bc->kind  |= kind;

          bc->val->U = q * Qratio;
          bc->val->S = Sspec;
        }
      }
    }

    if( Q == 0.0 )
    {
      sprintf( text, "\n                         %s %d\n",
                     "no discharge computed for line", no );
      REPORT::rpt.Output( text, 3 );
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kInlet
  ////////////////////////////////////////////////////////////////////////////////////////

  else if( isFS(kind, BCON::kInlet) )
  {
    if( ndat == 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            BCON* bc = FindBcon( nd->Getno(), b, bcon );

            if( !bc )
            {
              bcon[b].no   = nd->Getno();
              bcon[b].kind = kind;

              bcon[b].val->U = -U[0];        // Change the sign of flow on inlet, so
              bcon[b].val->V = -V[0];        // that the user may specify positive values
                                             // for q = H*U. SC, 17.01.2006.
              b++;
            }
            else if( bc->Consistent(kind) )
            {
              bc->kind  |= kind;

              bc->val->U = -U[0];
              bc->val->V = -V[0];
            }
          }
        }
      }
    }

    // -----------------------------------------------------------------------------------

    else // if( count > 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            double L;
            int j = FindInterval( el->nd[i]->x, el->nd[i]->y, ndat, x, y, &L );

            if( j >= 0 )
            {
              double Uspec = U[j]  +  L * (U[j+1] - U[j]);
              double Vspec = V[j]  +  L * (V[j+1] - V[j]);

              BCON* bc = FindBcon( nd->Getno(), b, bcon );

              if( !bc )
              {
                bcon[b].no   = nd->Getno();
                bcon[b].kind = kind;

                bcon[b].val->U = -Uspec;     // Change the sign of flow on inlet, so
                bcon[b].val->V = -Vspec;     // that the user may specify positive values
                                             // for q = H*U. SC, 17.01.2006.
                b++;
              }
              else if( bc->Consistent(kind) )
              {
                bc->kind  |= kind;

                bc->val->U = -Uspec;
                bc->val->V = -Vspec;
              }
            }

            else
            {
              if( firstMiss )
              {
                sprintf( text, "                 %s %d\n",
                               "no boundary condition at node",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );

                firstMiss = false;
              }

              else
              {
                sprintf( text, "                 %s %d\n",
                               "                             ",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );
              }
            }
          }
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////

  return b;
}


//////////////////////////////////////////////////////////////////////////////////////////

int BCONLINE::Outlet( MODEL*model, TIME *actualTime, int b, BCON *bcon, double *preQ )
{
  int  firstMiss = true;
  char text[120];

  GRID *ct = model->control;

  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kSQOutlet
  ////////////////////////////////////////////////////////////////////////////////////////

  if( isFS(kind, BCON::kSQOutlet) )
  {
    double Q = 0.0;

    // compute discharge through line
    for( int c=0; c<ct->Getne(); c++ )
    {
      ELEM* el = ct->Getelem(c);

      if( el->type != no )  continue;

      SHAPE* lShape = el->GetLShape();
      SHAPE* qShape = el->GetQShape();

      NODE* node[3];

      node[0] = el->nd[0];                 // corner nodes
      node[1] = el->nd[1];
      node[2] = el->nd[2];                 // midside node

      for( int g=0; g<qShape->ngp; g++ )   // loop on GAUSS points
      {
        double* M  = lShape->f[g];         // linear shape
        double* N  = qShape->f[g];         // quadratic shape
        double* dN = qShape->dfdx[g];

        // compute normal vector at Gauss point g
        // since the normal is not reduced to unit length it
        // implies the transformation of the integrand
        double nx =  dN[0]*node[0]->y + dN[1]*node[1]->y + dN[2]*node[2]->y;
        double ny = -dN[0]*node[0]->x - dN[1]*node[1]->x - dN[2]*node[2]->x;

        double H =    M[0] * (node[0]->v.S  -  node[0]->z)
                   +  M[1] * (node[1]->v.S  -  node[1]->z);

        if ( H <= 0.0 )  H = 0.0;

        double U = N[0]*node[0]->v.U + N[1]*node[1]->v.U + N[2]*node[2]->v.U;
        double V = N[0]*node[0]->v.V + N[1]*node[1]->v.V + N[2]*node[2]->v.V;

        double weight = qShape->weight[g];

        Q += weight * H * (U * nx  +  V * ny);
      }
    }

    Q = fabs( Q );   // not interested in flow direction !?

    //////////////////////////////////////////////////////////////////////////////////////
    // MPI: sum of discharge Q from all subdomains
#   ifdef _MPI_
    Q = model->subdom->Mpi_sum( Q );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    double  SQ = 0.0;
    double  SG = 0.0;

    double *LQ = U;
    double *LS = S;

    if( Q <= LQ[0]  ||  ndat <= 1 )
    {
      SQ = LS[0];
      SG = gct[0].So;
    }
    else
    {
      for( int j=1; j<ndat; j++ )
      {
        SQ = LS[j];
        SG = gct[j].So;

        if( Q <= LQ[j] )
        {
          SQ = LS[j-1]     + (LS[j]-LS[j-1])         * (Q-LQ[j-1]) / (LQ[j]-LQ[j-1]);
          SG = gct[j-1].So + (gct[j].So-gct[j-1].So) * (Q-LQ[j-1]) / (LQ[j]-LQ[j-1]);
          break;
        }
      }
    }

    // relaxed change of Q
    if( preQ  &&  *preQ > 0.0 )  SQ = (SQ + *preQ) / 2.0;

    for( int c=0; c<ct->Getne(); c++ )
    {
      ELEM* el = ct->Getelem(c);

      if( el->type != no ) continue;

      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        NODE* nd = el->nd[i];

        if( !nd->mark )
        {
          nd->mark = true;

          BCON* bc = FindBcon( nd->Getno(), b, bcon );

          if( !bc )
          {
            bcon[b].no   = nd->Getno();
            bcon[b].kind = kind;
            bcon[b].val->gct.node = nd->Getname();
            bcon[b].val->gct.nocg = gct[0].nocg;
            bcon[b].val->gct.So   = SG;

            if( preQ  &&  *preQ > 0.0 ) // relaxed change
            {
              bcon[b].val->S = (SQ + bcon[b].val->S) / 2.0;
            }
            else
            {
              bcon[b].val->S = SQ;
            }
            b++;
          }
          else if( bc->Consistent(kind) )
          {
            bc->kind |= kind;
            bc->val->gct.nocg = gct[0].nocg;
            bc->val->gct.nocg = SG;

            if( preQ  &&  *preQ > 0.0 ) // relaxed change
            {
              bc->val->S = (SQ + bcon[b].val->S) / 2.0;
            }
            else
            {
              bc->val->S = SQ;
            }
          }
        }
      }
    }

    if( preQ ) *preQ = SQ;
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kSTOutlet
  ////////////////////////////////////////////////////////////////////////////////////////

  else if( isFS(kind, BCON::kSTOutlet) )  // determine outlet condition S for actual time
  {
    double Sspec   = S[0];
    double gauge   = gct[0].nocg;
    double gaugeSo = gct[0].So;

    double L;
    int j = FindInterval( actualTime->Getsec(), 0.0, ndat, x, NULL, &L );


    if( j >= 0 )
    {
      Sspec   = S[j]       +  L * (S[j+1]      - S[j]);
      gauge   = gct[j].nocg;
      gaugeSo = gct[j].So  +  L * (gct[j+1].So - gct[j].So);
    }

    for( int c=0; c<ct->Getne(); c++ )
    {
      ELEM *el = ct->Getelem(c);

      if( el->type != no )  continue;

      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        NODE *nd = el->nd[i];

        if( !nd->mark )
        {
          nd->mark = true;

          BCON *bc = FindBcon( nd->Getno(), b, bcon );

          if( !bc )
          {
            bcon[b].no            = nd->Getno();
            bcon[b].kind          = kind;
            bcon[b].val->S        = Sspec;
            bcon[b].val->gct.node = nd->Getname();
            bcon[b].val->gct.nocg = gauge;
            bcon[b].val->gct.So   = gaugeSo;
            b++;
          }
          else if( bc->Consistent(kind) )
          {
            bc->kind         |= kind;
            bc->val->S        = Sspec;
            bc->val->gct.nocg = gauge;
            bc->val->gct.So   = gaugeSo;
          }
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kOutlet
  ////////////////////////////////////////////////////////////////////////////////////////

  else if( isFS(kind, BCON::kOutlet) )
  {
    if( ndat == 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM *el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE *nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            BCON *bc = FindBcon( nd->Getno(), b, bcon );

            if( !bc )
            {
              bcon[b].no            = nd->Getno();
              bcon[b].kind          = kind;
              bcon[b].val->S        = S[0];
              bcon[b].val->gct.node = nd->Getname();
              bcon[b].val->gct.nocg = gct[0].nocg;
              bcon[b].val->gct.So   = gct[0].So;
              b++;
            }
            else if( bc->Consistent(kind) )
            {
              bc->kind         |= kind;
              bc->val->S        = S[0];
              bc->val->gct.nocg = gct[0].nocg;
              bc->val->gct.So   = gct[0].So;
            }
          }
        }
      }
    }

    // -----------------------------------------------------------------------------------

    else // if( count > 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            double L;
            int j = FindInterval( el->nd[i]->x, el->nd[i]->y, ndat, x, y, &L );

            if( j >= 0 )
            {
              double Sspec = V[j]  +  L * (V[j+1] - V[j]);

              BCON* bc = FindBcon( nd->Getno(), b, bcon );

              if( !bc )
              {
                bcon[b].no            = nd->Getno();
                bcon[b].kind          = kind;
                bcon[b].val->S        = Sspec;
                bcon[b].val->gct.node = nd->Getname();
                bcon[b].val->gct.nocg = gct[0].nocg;
                bcon[b].val->gct.So   = gct[0].So;
                b++;
              }
              else if( bc->Consistent(kind) )
              {
                bc->kind         |= kind;
                bc->val->U        = Sspec;
                bc->val->gct.nocg = gct[0].nocg;
                bc->val->gct.So   = gct[0].So;
              }
            }

            else
            {
              if( firstMiss )
              {
                sprintf( text, "                 %s %d\n",
                               "no boundary condition at node",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );

                firstMiss = false;
              }

              else
              {
                sprintf( text, "                 %s %d\n",
                               "                             ",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );
              }
            }
          }
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////

  return b;
}


//////////////////////////////////////////////////////////////////////////////////////////

int BCONLINE::Further( MODEL* model, TIME* actualTime, int b, BCON* bcon, long* ki )
{
  int  firstMiss = true;
  char text[120];

  GRID* ct = model->control;
  GRID* rg = model->region;


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kOpenBnd
  ////////////////////////////////////////////////////////////////////////////////////////

  if( isFS(kind, BCON::kOpenBnd) )
  {
    for( int c=0; c<ct->Getne(); c++ )
    {
      ELEM* el = ct->Getelem(c);

      if( el->type != no )  continue;

      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        NODE* nd = el->nd[i];

        if( !nd->mark )
        {
          nd->mark = true;

          BCON* bc = FindBcon( nd->Getno(), b, bcon );

          if( !bc )
          {
            bcon[b].no   = nd->Getno();
            bcon[b].kind = kind;

            b++;
          }
          else if( bc->Consistent(kind) )
          {
            bc->kind  |= kind;
          }
        }
      }
    }

    CF( *ki, BCON::kOpenBnd );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kSlip || BCON::kSetUV
  ////////////////////////////////////////////////////////////////////////////////////////

  if( isFS(kind, BCON::kSlip)  ||  isFS(kind, BCON::kSetUV) )
  {
    if( ndat == 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            BCON* bc = FindBcon( nd->Getno(), b, bcon );

            if( !bc )
            {
              bcon[b].no   = nd->Getno();
              bcon[b].kind = kind;

              bcon[b].val->U = U[0];
              bcon[b].val->V = V[0];

              b++;
            }
            else if( bc->Consistent(kind) )
            {
              bc->kind  |= kind;

              bc->val->U = U[0];
              bc->val->V = V[0];
            }
          }
        }
      }
    }

    // ---------------------------------------------------------------------------------------------

    else // if( count > 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            double L;
            int j = FindInterval( el->nd[i]->x, el->nd[i]->y, ndat, x, y, &L );

            if( j >= 0 )
            {
              double Uspec = U[j]  +  L * (U[j+1] - U[j]);
              double Vspec = V[j]  +  L * (V[j+1] - V[j]);

              BCON* bc = FindBcon( nd->Getno(), b, bcon );

              if( !bc )
              {
                bcon[b].no   = nd->Getno();
                bcon[b].kind = kind;

                bcon[b].val->U = Uspec;
                bcon[b].val->V = Vspec;

                b++;
              }
              else if( bc->Consistent(kind) )
              {
                bc->kind  |= kind;

                bc->val->U = Uspec;
                bc->val->V = Vspec;
              }
            }

            else
            {
              if( firstMiss )
              {
                sprintf( text, "                 %s %d\n",
                               "no boundary condition at node",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );

                firstMiss = false;
              }

              else
              {
                sprintf( text, "                 %s %d\n",
                               "                             ",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );
              }
            }
          }
        }
      }
    }

    CF( *ki, BCON::kSlip | BCON::kSetUV );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kLoglaw
  ////////////////////////////////////////////////////////////////////////////////////////

  if( isFS(kind, BCON::kLoglaw) )
  {
    for( int c=0; c<ct->Getne(); c++ )
    {
      ELEM* el = ct->Getelem(c);

      if( el->type != no )  continue;

      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        NODE* nd = el->nd[i];

        if( !nd->mark )
        {
          nd->mark = true;

          BCON* bc = FindBcon( nd->Getno(), b, bcon );

          if( !bc )
          {
            bcon[b].no   = nd->Getno();
            bcon[b].kind = kind;

            bcon[b].val->dw = dw;
            bcon[b].val->kw = kw;

            b++;
          }
          else if( bc->Consistent(kind) )
          {
            bc->kind  |= kind;

            bc->val->dw = dw;
            bc->val->kw = kw;
          }
        }
      }
    }

    CF( *ki, BCON::kLoglaw );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kSetS
  ////////////////////////////////////////////////////////////////////////////////////////

  else if( isFS(kind, BCON::kSetS) )
  {
    if( ndat == 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            BCON* bc = FindBcon( nd->Getno(), b, bcon );

            if( !bc )
            {
              bcon[b].no     = nd->Getno();
              bcon[b].kind   = kind;
              bcon[b].val->S = S[0];
              b++;
            }
            else if( bc->Consistent(kind) )
            {
              bc->kind  |= kind;
              bc->val->S = S[0];
            }
          }
        }
      }
    }

    // ---------------------------------------------------------------------------------------------

    else // if( count > 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            double L;
            int j = FindInterval( nd->x, nd->y, ndat, x, y, &L );

            if( j >= 0 )
            {
              double Sspec = S[j]  +  L * (S[j+1] - S[j]);

              BCON* bc = FindBcon( nd->Getno(), b, bcon );

              if( !bc )
              {
                bcon[b].no     = nd->Getno();
                bcon[b].kind   = kind;
                bcon[b].val->S = Sspec;
                b++;
              }
              else if( bc->Consistent(kind) )
              {
                bc->kind  |= kind;
                bc->val->S = Sspec;
              }
            }

            else
            {
              if( firstMiss )
              {
                sprintf( text, "                 %s %d\n",
                               "no boundary condition at node",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );

                firstMiss = false;
              }

              else
              {
                sprintf( text, "                 %s %d\n",
                               "                             ",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );
              }
            }
          }
        }
      }
    }

    CF( *ki, BCON::kSetS );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kSetKD
  ////////////////////////////////////////////////////////////////////////////////////////

  else if( isFS(kind, BCON::kSetKD) )
  {
    if( ndat == 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            BCON* bc = FindBcon( nd->Getno(), b, bcon );

            if( !bc )
            {
              bcon[b].no     = nd->Getno();
              bcon[b].kind   = kind;
              bcon[b].val->K = K[0];
              bcon[b].val->D = D[0];
              b++;
            }
            else if( bc->Consistent(kind) )
            {
              bc->kind  |= kind;
              bc->val->K = K[0];
              bc->val->D = D[0];
            }
          }
        }
      }
    }

    // ---------------------------------------------------------------------------------------------

    else // if( count > 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            double L;
            int j = FindInterval( nd->x, nd->y, ndat, x, y, &L );

            if( j >= 0 )
            {
              double Kspec = K[j]  +  L * (K[j+1] - K[j]);
              double Dspec = D[j]  +  L * (D[j+1] - D[j]);

              BCON* bc = FindBcon( nd->Getno(), b, bcon );

              if( !bc )
              {
                bcon[b].no     = nd->Getno();
                bcon[b].kind   = kind;
                bcon[b].val->K = Kspec;
                bcon[b].val->D = Dspec;
                b++;
              }
              else if( bc->Consistent(kind) )
              {
                bc->kind  |= kind;
                bc->val->K = Kspec;
                bc->val->D = Dspec;
              }
            }

            else
            {
              if( firstMiss )
              {
                sprintf( text, "                 %s %d\n",
                               "no boundary condition at node",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );

                firstMiss = false;
              }

              else
              {
                sprintf( text, "                 %s %d\n",
                               "                             ",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );
              }
            }
          }
        }
      }
    }

    CF( *ki, BCON::kSetKD );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kSetC
  ////////////////////////////////////////////////////////////////////////////////////////

  else if( isFS(kind, BCON::kSetC) )
  {
    if( ndat == 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            BCON* bc = FindBcon( nd->Getno(), b, bcon );

            if( !bc )
            {
              bcon[b].no        = nd->Getno();
              bcon[b].kind      = kind;
              bcon[b].val->C[0] = C[0];
              b++;
            }
            else if( bc->Consistent(kind) )
            {
              bc->kind     |= kind;
              bc->val->C[0] = C[0];
            }
          }
        }
      }
    }

    // ---------------------------------------------------------------------------------------------

    else // if( count > 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            double L;
            int j = FindInterval( nd->x, nd->y, ndat, x, y, &L );

            if( j >= 0 )
            {
              double Cspec = C[j]  +  L * (C[j+1] - C[j]);

              BCON* bc = FindBcon( nd->Getno(), b, bcon );

              if( !bc )
              {
                bcon[b].no        = nd->Getno();
                bcon[b].kind      = kind;
                bcon[b].val->C[0] = Cspec;
                b++;
              }
              else if( bc->Consistent(kind) )
              {
                bc->kind     |= kind;
                bc->val->C[0] = Cspec;
              }
            }

            else
            {
              if( firstMiss )
              {
                sprintf( text, "                 %s %d\n",
                               "no boundary condition at node",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );

                firstMiss = false;
              }

              else
              {
                sprintf( text, "                 %s %d\n",
                               "                             ",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );
              }
            }
          }
        }
      }
    }

    CF( *ki, BCON::kSetC );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // BCON::kRateC
  ////////////////////////////////////////////////////////////////////////////////////////

  else if( isFS(kind, BCON::kRateC) )
  {
    if( ndat == 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            BCON* bc = FindBcon( nd->Getno(), b, bcon );

            if( !bc )
            {
              bcon[b].no         = nd->Getno();
              bcon[b].kind       = kind;
              bcon[b].val->Qb[0] = Qb[0];
              b++;
            }
            else if( bc->Consistent(kind) )
            {
              bc->kind      |= kind;
              bc->val->Qb[0] = Qb[0];
            }
          }
        }
      }
    }

    // ---------------------------------------------------------------------------------------------

    else // if( count > 1 )
    {
      for( int c=0; c<ct->Getne(); c++ )
      {
        ELEM* el = ct->Getelem(c);

        if( el->type != no )  continue;

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          NODE* nd = el->nd[i];

          if( !nd->mark )
          {
            nd->mark = true;

            double L;
            int j = FindInterval( nd->x, nd->y, ndat, x, y, &L );

            if( j >= 0 )
            {
              double Qbspec = Qb[j]  +  L * (Qb[j+1] - Qb[j]);

              BCON* bc = FindBcon( nd->Getno(), b, bcon );

              if( !bc )
              {
                bcon[b].no         = nd->Getno();
                bcon[b].kind       = kind;
                bcon[b].val->Qb[0] = Qbspec;
                b++;
              }
              else if( bc->Consistent(kind) )
              {
                bc->kind      |= kind;
                bc->val->Qb[0] = Qbspec;
              }
            }

            else
            {
              if( firstMiss )
              {
                sprintf( text, "                 %s %d\n",
                               "no boundary condition at node",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );

                firstMiss = false;
              }

              else
              {
                sprintf( text, "                 %s %d\n",
                               "                             ",
                               el->nd[i]->Getname() );
                REPORT::rpt.Output( text, 3 );
              }
            }
          }
        }
      }
    }

    CF( *ki, BCON::kRateC );
  }

  ////////////////////////////////////////////////////////////////////////////////////////

  else
  {
    *ki = 0l;
  }

  ////////////////////////////////////////////////////////////////////////////////////////

  return b;
}
