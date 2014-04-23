// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_SL2D
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C) 2011 - 2014  by  P.M. SCHROEDER  (sc)
//
// This program is free software; you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with this program; if
// not, write to the
//
// Free Software Foundation, Inc.
// 59 Temple Place
// Suite 330
// Boston
// MA 02111-1307 USA
//
// -------------------------------------------------------------------------------------------------
//
// P.M. Schroeder
// Walzbachtal / Germany
// michael.schroeder@hnware.de
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

// #define DEBUG 1

#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsSL2D.h"


int EQS_SL2D::Coefs( ELEM*    elem,
                     PROJECT* project,
                     double** estifm,
                     double*  force )
{
  if( isFS(elem->flag, ELEM::kDry) )  return 0;

  if( isFS(elem->flag,ELEM::kBound) )
    Bound( elem, project, estifm, force );
  else
    Region( elem, project, estifm, force );

  return 1;
}



void EQS_SL2D::Bound( ELEM*    elem,
                      PROJECT* project,
                      double** estifm,
                      double*  force )
{
  SED* sed = &project->sed;

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int ncn = lShape->nnd;
  int nnd = qShape->nnd;
  int ngp = qShape->ngp;

  if( force )
  {
    for( int i=0; i<maxEleq; i++ )  force[i] = 0.0;
  }

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      for( int j=0; j<maxEleq; j++ )
      {
        estifm[i][j] = 0.0;
      }
    }
  }

  NODE* node[3];

  node[0] = elem->nd[0];      // corner nodes
  node[1] = elem->nd[1];
  node[2] = elem->nd[2];      // midside node


  if(    !isFS(node[2]->bc.kind, BCON::kRateC) )  return;
      //&& !isFS(node[2]->bc.kind, BCON::kInlet)
      //&& !isFS(node[2]->bc.kind, BCON::kOutlet) )  return;


  for( int g=0; g<ngp; g++ )            // loop on GAUSS points
  {
    double* m  = lShape->f[g];          // linear shape
    double* n  = qShape->f[g];          // quadratic shape
    double* dn = qShape->dfdx[g];

    // weight of Gauss point g -----------------------------------------------------------

    double weight = qShape->weight[g];

    // compute normal vector -------------------------------------------------------------
    // since the normal is not reduced to unit length it
    // implies the transformation of the integrand

    double nx = 0.0;
    double ny = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      nx += dn[i]*node[i]->y;
      ny -= dn[i]*node[i]->x;
    }

    // compute some parameter at Gauss point ---------------------------------------------

    double H = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      H += m[i] * (node[i]->v.S - node[i]->z);
    }

    double U = 0.0;
    double V = 0.0;
    double C = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      U += n[i] * node[i]->v.U;
      V += n[i] * node[i]->v.V;
      C += n[i] * node[i]->v.C;
    }

    // compute sediment transport --------------------------------------------------------

    double qs = 0.0;

    if( isFS(node[2]->bc.kind, BCON::kRateC) )
    {
      // use experimental upstream sediment transport ... --------------------------------
      for( int i=0; i<nnd; i++ )  qs += n[i] * node[i]->bc.val->C[0];
    }
    else
    {
      // ... or compute sediment transport at gauss point --------------------------------
      qs = H * (U*nx + V*ny) * C;
    }


    // force vector ----------------------------------------------------------------------

    if( force )
    {
      double f = weight * qs;

      for( int i=0; i<nnd; i++ )  force[i] -= n[i] * f;
    }

    // compute estifm, if no transport rate is specified ---------------------------------

    if( estifm )
    {
      if( !isFS(node[2]->bc.kind, BCON::kRateC) )
      {
        double df = weight * H * (U*nx + V*ny);

        for( int i=0; i<nnd; i++ )
        {
          for( int j=0; j<nnd; j++ )
          {
            estifm[i][i] += n[i] * df * n[j];
          }
        }
      }
    }
  }
}


void EQS_SL2D::Region( ELEM*    elem,
                       PROJECT* project,
                       double** estifm,
                       double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  double dt  = project->timeint.incTime.Getsec();
  double adt = project->timeint.thetaSedi * dt;
  double pdt = (1.0 -  project->timeint.thetaSedi) * dt;

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int ngp = qShape->ngp;            /* number of GAUSS points */
  int nnd = qShape->nnd;            /* total number of nodes */
  int ncn = lShape->nnd;            /* number of corner nodes */


  if( force )  for( int i=0; i<maxEleq; i++ )  force[i] = 0.0;

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      double* estifmPtr = estifm[i];

      for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    x[i] = elem->nd[i]->x - *x;
    y[i] = elem->nd[i]->y - *y;
  }

  x[0] = y[0] = 0.0;

  TYPE* type = TYPE::Getid( elem->type );

  double rho  = project->rho;

  double taus = project->sed.taus;
  double tauc = project->sed.tauc;
  double us   = project->sed.us;
  double M    = project->sed.M;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve adv equation

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double  trafo[2][2];
    double* dfdxPtr = qShape->dfdx[g];
    double* dfdyPtr = qShape->dfdy[g];

    double detj = qShape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );

    double weight = detj * qShape->weight[g];


    // -----------------------------------------------------------------------------------
    // compute values of linear shape functions at GP i

    double* m = lShape->f[g];


    // -----------------------------------------------------------------------------------
    // compute values of quadratic shape functions at GP i

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = qShape->f[g];

    for( int j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j] + trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j] + trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP i

    // integrate with linear shape -------------------------------------------------------

    double h  = 0.0;
    double hl = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE *node = elem->nd[j];

      h  += m[j] * (node->v.S - node->z);
      hl += m[j] * (node->zor - node->zero);
    }

    if( h < 0.0 )  h = 0.0;


    // integrate with quadratic shape ----------------------------------------------------

    double U = 0.0;
    double V = 0.0;

    double C     = 0.0;
    double dCdx  = 0.0;
    double dCdy  = 0.0;
    double Co    = 0.0;
    double dCodx = 0.0;
    double dCody = 0.0;

    double Exx   = 0.0;
    double Exy   = 0.0;
    double Eyy   = 0.0;

    for( int j=0; j<nnd; j++ )
    {
      NODE *node = elem->nd[j];

      double ndC  = node->v.C;
      double ndCo = node->vo.C;

      Exx += n[j] * node->exx * node->vt;
      Exy += n[j] * node->exy * node->vt;
      Eyy += n[j] * node->eyy * node->vt;

      U     +=    n[j] * node->v.U;
      V     +=    n[j] * node->v.V;

      C     +=    n[j] * ndC;
      dCdx  += dndx[j] * ndC;
      dCdy  += dndy[j] * ndC;

      Co    +=    n[j] * ndCo;
      dCodx += dndx[j] * ndCo;
      dCody += dndy[j] * ndCo;
    }

    // compute shear stress --------------------------------------------------------------

    double Us   = sqrt( U*U + V*V );
    double Utau = project->sed.GetUtau( Us, h, 0.0, project );

    double tau  = project->rho * Utau * Utau;


    // -----------------------------------------------------------------------------------
    // compute right side

    if( force )
    {
      double f, fx, fy;

      f   = h * (C - Co);                                   // time

      //f  += h * adt * (U * dCdx   +  V * dCdy);             // advection
      //f  += h * pdt * (U * dCodx  +  V * dCody);

      fx  = h * adt * (Exx * dCdx   +  Exy * dCdy);         // diffusion
      fx += h * pdt * (Exx * dCodx  +  Exy * dCody);

      fy  = h * adt * (Exy * dCdx   +  Eyy * dCdy);
      fy += h * pdt * (Exy * dCodx  +  Eyy * dCody);

      fx -= h * U * (adt * C + pdt * Co);                   // advection
      fy -= h * V * (adt * C + pdt * Co);

      if( tau < taus )                                      // sedimentation
      {
        f += adt * (1.0 - tau/taus) * us * C;
        f += pdt * (1.0 - tau/taus) * us * Co;
      }

      if( tau > tauc  &&  hl > 0.0 )                        // erosion
      {
        f -= dt * (tau/tauc - 1.0) * M;
      }

      f  *= weight;
      fx *= weight;
      fy *= weight;

      for( int j=0; j<nnd; j++ )
      {
        force[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }
    }


    // -----------------------------------------------------------------------------------
    // compute components of NEWTON-RAPHSON Jacobi matrix

    if( estifm )
    {
      double t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double df__;//, df_x, df_y;
      double dfx_, dfxx, dfxy;
      double dfy_, dfyx, dfyy;

      df__  = weight * h;

      if( tau < taus )
      {
        df__ += weight * adt * ( 1.0 - tau/taus ) * us;
      }

      //df_x  = weight * h * adt * U;
      //df_y  = weight * h * adt * V;

      dfx_  = -weight * h * adt * U;
      dfy_  = -weight * h * adt * V;

      dfxx  =  weight * h * adt * Exx;
      dfxy  =
      dfyx  =  weight * h * adt * Exy;
      dfyy  =  weight * h * adt * Eyy;

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j];//  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j]  +  dfxx * dndx[j]  +  dfxy * dndy[j];
        ty[j] = dfy_ * n[j]  +  dfyx * dndx[j]  +  dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        double* estifmPtr = estifm[j];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }
    }
  } // end of loop on all Gauss points


  // set specified flow conditions (sink or source) --------------------------------------

  //if( force )
  //{
  //  for( int i=0; i<nnd; i++ )
  //  {
  //    BCON* bc = elem->nd[i]->bc;

  //    if( bc  &&  isFS(bc->kind, BCON::kFlowC) )
  //    {
  //      force[i] += bc->spec.C;
  //    }
  //  }
  //}
}
