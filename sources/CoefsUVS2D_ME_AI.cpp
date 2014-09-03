// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_ME_AI
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

#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsUVS2D_ME_AI.h"

//#define kDebug

////////////////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_ME_AI::Region( ELEM     *elem,
                              PROJECT  *project,
                              double  **estifm,
                              double   *force )
{
  // -----------------------------------------------------------------------------------------------
  // initializations

  SHAPE *lShape = elem->GetLShape();
  SHAPE *bShape = elem->GetBShape();

  TYPE *type = TYPE::Getid( elem->type );

  int ngp = lShape->ngp;                // number of GAUSS points

  int nnd = elem->GetQShape()->nnd;     // number of nodes
  int ncn = lShape->nnd;                // number of corner nodes
  int nbn = bShape->nnd;                // number of nodes for bubble shape function

  double gravity = project->g;

//int eqidU = 0;
  int eqidV = nnd;
  int eqidS = 2 * nnd;

  if( force )
  {
    for( int i=0; i<maxEleq; i++ )  force[i] = 0.0;
  }

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      for( int j=0; j<maxEleq; j++ )  estifm[i][j] = 0.0;
    }
  }


  // -----------------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<ncn; i++ )
  {
    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;


  // -----------------------------------------------------------------------------------------------
  // use GAUSS point integration to solve momentum equations for x- and
  // y-direction (U and V) and continuity equation (H)

  double area = 0.0;

  for( int g=0; g<ngp; g++ )
  {
    // ---------------------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double* dfdxPtr = lShape->dfdx[g];
    double* dfdyPtr = lShape->dfdy[g];

    double  trafo[2][2];

    double detj   = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
    double weight = detj * lShape->weight[g];

    area += weight;


    // ---------------------------------------------------------------------------------------------
    // compute values of linear shape functions at GP g

    double *m = lShape->f[g];

    double  dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    for( int i=0; i<ncn; i++ )
    {
      dmdx[i] = trafo[0][0] * dfdxPtr[i] + trafo[0][1] * dfdyPtr[i];
      dmdy[i] = trafo[1][0] * dfdxPtr[i] + trafo[1][1] * dfdyPtr[i];
    }


    // ---------------------------------------------------------------------------------------------
    // compute values of bubble shape functions at GP g

    double *b = bShape->f[g];

    double  dbdx[kMaxNodes2D], dbdy[kMaxNodes2D];

    for( int i=0; i<nbn; i++ )
    {
      double dfdx = bShape->dfdx[g][i];
      double dfdy = bShape->dfdy[g][i];

      dbdx[i] = trafo[0][0] * dfdx + trafo[0][1] * dfdy;
      dbdy[i] = trafo[1][0] * dfdx + trafo[1][1] * dfdy;
    }


    // ----------------------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP g
    //             horizontal velocities: U and V
    //             flow depth           : H
    //             Source or Sink       : SS
    //             bottom elevation     : a
    //             eddy viscosity       : vt, cf
    //             Reynolds stresses    : uu, uv, and vv

    // integrate H, a and cf with linear shape

    double H    = 0.0;
    double dHdt = 0.0;
    double dHdx = 0.0;
    double dHdy = 0.0;

    double dadx = 0.0;
    double dady = 0.0;

    double cf   = 0.0;

    double SS   = 0.0;

    double uu   = 0.0;
    double uv   = 0.0;
    double vv   = 0.0;

    double dUdt = 0.0;
    double dVdt = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      NODE *node = elem->nd[i];
      BCON *bcon = &node->bc;

      double ndZ = node->z;
      double ndH = node->v.S - ndZ;

      if( ndH <= 0.0 )  ndH = project->hmin;

      double ndSS;

      // -------------------------------------------------------------------------------------------
      // Source or Sink
      //
      // ndSS = -Q * ncn / A;
      //
      // bcon->val->Q  - Q
      // bcon->val->A  - area of connected elements (see BCONSET::InitBcon)
      // ncn           - number of corner nodes
      //
      if( isFS(bcon->kind, BCON::kSource) )
      {
        ndSS = -bcon->val->Q * ncn / bcon->val->A;
      }
      else
      {
        ndSS = 0.0;
      }
      // -------------------------------------------------------------------------------------------

      H    +=    m[i] * ndH;
      dHdx += dmdx[i] * ndH;
      dHdy += dmdy[i] * ndH;

      dadx += dmdx[i] * ndZ;
      dady += dmdy[i] * ndZ;

      dHdt +=    m[i] * node->v.dSdt;

      cf   +=    m[i] * node->cf;

      SS   +=    m[i] * ndSS;            // Source or Sink

      uu   +=    m[i] * node->uu;
      uv   +=    m[i] * node->uv;
      vv   +=    m[i] * node->vv;

      dUdt +=    m[i] * node->v.dUdt;
      dVdt +=    m[i] * node->v.dVdt;
    }

//  if( H <= 0.0 )  H = project->hmin;


    // ---------------------------------------------------------------------------------------------

    double U    = 0.0;
    double dUdx = 0.0;
    double dUdy = 0.0;

    double V    = 0.0;
    double dVdx = 0.0;
    double dVdy = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      NODE*  node = elem->nd[i];

      double ndU  = node->v.U;
      double ndV  = node->v.V;

      U    +=    b[i] * ndU;
      dUdx += dbdx[i] * ndU;
      dUdy += dbdy[i] * ndU;

      V    +=    b[i] * ndV;
      dVdx += dbdx[i] * ndV;
      dVdy += dbdy[i] * ndV;
    }

    U    +=    b[ncn] * elem->U;
    dUdx += dbdx[ncn] * elem->U;
    dUdy += dbdy[ncn] * elem->U;

    V    +=    b[ncn] * elem->V;
    dVdx += dbdx[ncn] * elem->V;
    dVdy += dbdy[ncn] * elem->V;

    // compute dispersion coefficients -------------------------------------------------------------

    double Dxx = 0.0;
    double Dxy = 0.0;
    double Dyy = 0.0;

    if( project->actualDisp > 0 )
    {
      for( int i=0; i<ncn; i++ )
      {
        NODE* node = elem->nd[i];

        Dxx += m[i] * node->Dxx;
        Dxy += m[i] * node->Dxy;
        Dyy += m[i] * node->Dyy;
      }
    }


    // compute eddy viscosity according to ELDER's assumption ----------------------------

    double vtxx = 0.0;
    double vtxy = 0.0;
    double vtyy = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      NODE* node = elem->nd[i];

      vtxx += m[i] * node->exx * node->vt;
      vtxy += m[i] * node->exy * node->vt;
      vtyy += m[i] * node->eyy * node->vt;
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vtxx < type->vt )  vtxx = type->vt;
      if( vtyy < type->vt )  vtyy = type->vt;
    }

    // add eddy viscosity to kinematic viscosity -----------------------------------------

    vtxx += project->vk;
    vtxy += project->vk;
    vtyy += project->vk;


    // ---------------------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double Ures = sqrt( U*U + V*V );           // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double *forcePtr;


      // -------------------------------------------------------------------------------------------
      // compute x-momentum

      f   = H * dUdt;                                  // time
      f  += H * (U * dUdx  +  V * dUdy);               // convection

      fx  = H * (vtxx * dUdx + vtxy * dUdy);           // eddy viscosity (approximation)
      fy  = H * (vtxy * dUdx + vtyy * dUdy);

      fx -= H * uu;                                    // turbulence
      fy -= H * uv;

      f  += H * gravity * dadx;                        // gravity
      fx -= H * H * gravity / 2.0;

      f  += cf * Ures * U;                             // bottom friction

      fx += H * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );   // dispersion
      fy += H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int i=0; i<nbn; i++ )
      {
        forcePtr[i] -= b[i] * f  +  dbdx[i] * fx  +  dbdy[i] * fy;
      }


      // -------------------------------------------------------------------------------------------
      // compute y-momentum

      f   = H * dVdt;                                  // time
      f  += H * (U * dVdx  +  V * dVdy);               // convection

      fx  = H * (vtxx * dVdx + vtxy * dVdy);           // eddy viscosity (approximation)
      fy  = H * (vtxy * dVdx + vtyy * dVdy);

      fx -= H * uv;                                    // turbulence
      fy -= H * vv;

      f  += H * gravity * dady;                        // gravity
      fy -= H * H * gravity / 2.0;

      f  += cf * Ures * V;                             // bottom friction

      fx += H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );     // dispersion
      fy += H * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + eqidV;

      for( int i=0; i<nbn; i++ )
      {
        forcePtr[i] -= b[i] * f  +  dbdx[i] * fx  +  dbdy[i] * fy;
      }


      // -------------------------------------------------------------------------------------------
      // compute continuity equation

      f  = dHdt  +  H * (dUdx + dVdy)  +  U * dHdx  +  V * dHdy;
      f += SS;      // source or sink

      f *= weight;

      forcePtr = force + eqidS;

      for( int i=0; i<ncn; i++ )
      {
        forcePtr[i] -= m[i] * f;
      }
    }


    if( estifm )
    {
      double t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double df__, df_x, df_y, dfx_, dfy_, dfxx, dfxy, dfyx, dfyy;


      // -------------------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      double iUres;

      if( Ures > 1.0e-8 ) iUres = 1.0 / Ures;
      else                iUres = 0.0;


      // --- U-derivative of x-momentum ------------------------------------------------------------

      df__  =  weight * H * relaxThdt_UV;

      df__ +=  weight * H * dUdx;
      df_x  =  weight * H * U;
      df_y  =  weight * H * V;

      dfxx  =  weight * H * vtxx;
      dfxy  =  weight * H * vtxy;

      dfyx  =  weight * H * vtxy;
      dfyy  =  weight * H * vtyy;

      df__ +=  weight * cf * (iUres * U*U  +  Ures);

      dfx_  =  weight * H * ( 2.0*U*Dxx - 2.0*V*Dxy );
      dfy_  =  weight * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        t[i]  = df__ * b[i]  +  df_x * dbdx[i]  +  df_y * dbdy[i];
        tx[i] = dfx_ * b[i]  +  dfxx * dbdx[i]  +  dfxy * dbdy[i];
        ty[i] = dfy_ * b[i]  +  dfyx * dbdx[i]  +  dfyy * dbdy[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        double *estifmPtr = estifm[i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // --- V-derivative of x-momentum ------------------------------------------------------------

      df__  =  weight * H * dUdy;

      df__ +=  weight * cf * iUres * U * V;

      dfx_  =  weight * H * ( 2.0*V*Dyy - 2.0*V*Dxy );
      dfy_  =  weight * H * ( V*(Dxx-Dyy) - 2.0*V*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        t[i]  = df__ * b[i];
        tx[i] = dfx_ * b[i];
        ty[i] = dfy_ * b[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        double *estifmPtr = estifm[i] + eqidV;

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // --- H-derivative of x-momentum ------------------------------------------------------------

      df__  =  weight * dUdt;

      df__ +=  weight * (U * dUdx + V * dUdy);

      dfx_  =  weight * (vtxx * dUdx  +  vtxy * dUdy);
      dfy_  =  weight * (vtxy * dUdx  +  vtyy * dUdy);

      dfx_ -=  weight * uu;
      dfy_ -=  weight * uv;

      df__ +=  weight * gravity * dadx;
      dfx_ -=  weight * gravity * H;

      dfx_ +=  weight * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );
      dfy_ +=  weight * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( int i=0; i<nbn; i++ )
      {
        double *estifmPtr = estifm[i] + eqidS;

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // --- U-derivative of y-momentum ------------------------------------------------------------

      df__  =  weight * H * dVdx;

      df__ +=  weight * cf * iUres * U * V;

      dfx_  =  weight * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );
      dfy_  =  weight * H * ( 2.0*V*Dxy + 2.0*U*Dyy );

      for( int i=0; i<nbn; i++ )
      {
        t[i]  = df__ * b[i];
        tx[i] = dfx_ * b[i];
        ty[i] = dfy_ * b[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        double *estifmPtr = estifm[eqidV + i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // --- V-derivative of y-momentum ------------------------------------------------------------

      df__  =  weight * H * relaxThdt_UV;

      df__ +=  weight * H * dVdy;
      df_x  =  weight * H * U;
      df_y  =  weight * H * V;

      dfxx  =  weight * H * vtxx;
      dfxy  =  weight * H * vtxy;

      dfyx  =  weight * H * vtxy;
      dfyy  =  weight * H * vtyy;

      df__ +=  weight * cf * (iUres * V*V  +  Ures);

      dfx_  =  weight * H * ( U*(Dxx-Dyy) - 2.0*V*Dxy );
      dfy_  =  weight * H * ( 2.0*V*Dxx + 2.0*U*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        t[i]  = df__ * b[i]  +  df_x * dbdx[i]  +  df_y * dbdy[i];
        tx[i] = dfx_ * b[i]  +  dfxx * dbdx[i]  +  dfxy * dbdy[i];
        ty[i] = dfy_ * b[i]  +  dfyx * dbdx[i]  +  dfyy * dbdy[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        double *estifmPtr = estifm[eqidV + i] + eqidV;

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // H-derivative of y-momentum ----------------------------------------------------------------

      df__  =  weight * dVdt;

      df__ +=  weight * (U * dVdx + V * dVdy);

      dfx_  =  weight * (vtxx * dVdx + vtxy * dVdy);
      dfy_  =  weight * (vtxy * dVdx + vtyy * dVdy);

      dfx_ -=  weight * uv;
      dfy_ -=  weight * vv;

      df__ +=  weight * gravity * dady;
      dfy_ -=  weight * gravity * H;

      dfx_ +=  weight * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );
      dfy_ +=  weight * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      for( int i=0; i<ncn; i++ )
      {
        t[i]  = df__ * m[i];
        tx[i] = dfx_ * m[i];
        ty[i] = dfy_ * m[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        double *estifmPtr = estifm[eqidV + i] + eqidS;

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // U-derivative of continuity ----------------------------------------------------------------

      df__ = weight * dHdx;
      df_x = weight * H;

      for( int i=0; i<nbn; i++ )
      {
        t[i] = df__ * b[i]  +  df_x * dbdx[i];
      }

      for( int i=0; i<ncn; i++ )
      {
        double *estifmPtr = estifm[eqidS + i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += m[i]*t[j];
        }
      }


      // V-derivative of continuity ----------------------------------------------------------------

      df__ = weight * dHdy;
      df_y = weight * H;

      for( int i=0; i<nbn; i++ )
      {
        t[i] = df__ * b[i] + df_y * dbdy[i];
      }

      for( int i=0; i<ncn; i++ )
      {
        double *estifmPtr = estifm[eqidS + i] + eqidV;

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += m[i]*t[j];
        }
      }


      // H-derivative of continuity ----------------------------------------------------------------

      df__  = weight * relaxThdt_H;
      df__ += weight * (dUdx + dVdy);
      df_x  = weight * U;
      df_y  = weight * V;

      for( int j=0; j<ncn; j++ )
      {
        t[j] = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
      }

      for( int i=0; i<ncn; i++ )
      {
        double *estifmPtr = estifm[eqidS + i] + eqidS;

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[j] += m[i]*t[j];
        }
      }
    }
  } // END of loop over all GAUSS points


  // -----------------------------------------------------------------------------------------------
  // apply transformation

  Rotate2D( ncn, eqidV, elem->nd, estifm, force );


  // -----------------------------------------------------------------------------------------------
  // eliminate the two equation at the bubble node
  Eliminate( elem, estifm, force, ncn, nbn, eqidV, eqidS );


  // -----------------------------------------------------------------------------------------------
  // insert experimental upstream boundary forces at nodes with inflow
  // boundary condition  (qfix = constant):

  for( int i=0; i<ncn; i++ )
  {
    BCON* bcon = &elem->nd[i]->bc;

    if( isFS(bcon->kind,BCON::kInlet) )
    {
      NODE* node = elem->nd[i];

      // set equation row i to zero
      if( estifm )
      {
        for( int j=0; j<maxEleq; j++ )  estifm[i][j] = 0.0;
      }

      // compute flow at node i in normal direction
      double Un    = node->v.U * bcon->niox  +  node->v.V * bcon->nioy;
      double H     = node->v.S - node->z;
      double specQ = bcon->val->U;

      if( H <= 0.0 )
      {
        H = project->hmin;
        specQ = 0.0;
      }

      // force vector
      if( force )  force[i] = area * (specQ - Un*H);

      // stiffness matrix
      if( estifm )
      {
        estifm[i][i]         = area * H;
        estifm[i][eqidS + i] = area * Un;
      }
    }
  }
}
