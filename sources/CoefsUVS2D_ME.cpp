// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_ME
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

#include "EqsUVS2D_ME.h"


void EQS_UVS2D_ME::Bound( ELEM*    elem,
                          PROJECT* project,
                          double** estifm,
                          double*  force )
{
  SHAPE* lShape  = elem->GetLShape();

  int    ncn     = lShape->nnd;
  int    nnd     = elem->GetQShape()->nnd;

  int    startV  = nnd;
  int    startS  = 2 * nnd;

  double gravity = project->g;

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


  int solidFlag;

  if(    !isFS(elem->nd[2]->bc.kind, BCON::kLoglaw)
      ||  isFS(elem->flag, ELEM::kInlet)
      ||  isFS(elem->flag, ELEM::kOutlet)
      ||  isFS(elem->flag, ELEM::kOpenBnd) )
  {
    solidFlag = false;
  }

  else
  {
    solidFlag = true;
  }


  NODE* node[3];

  node[0] = elem->nd[0];      // corner nodes
  node[1] = elem->nd[1];


  // row index for force vector and element stiffness matrix -----------------------------

  int r0U = 0;        int r1U = 1;
  int r0V = startV;   int r1V = startV + 1;


  // column index for element stiffness matrix -------------------------------------------

  int c0H = startS;   int c1H = startS + 1;


  for( int g=0; g<lShape->ngp; g++ )   // loop on GAUSS points
  {
    double *m  = lShape->f[g];         // linear shape
    double *dm = lShape->dfdx[g];


    // test midside node for boundary condition outFlow ----------------------------------

    double H;

    if( isFS(node[2]->bc.kind, BCON::kOutlet) )
    {
      // use experimental downstream flow depth ... --------------------------------------
      double H0 = node[0]->bc.val->S - node[0]->z;
      double H1 = node[1]->bc.val->S - node[1]->z;

      if( H0 <= 0.0 )  H0 = project->hmin;
      if( H1 <= 0.0 )  H1 = project->hmin;

      H = m[0] * H0  +  m[1] * H1;
    }
    else
    {
     // ... or compute flow depth at gauss point -----------------------------------------
      double H0 = node[0]->v.S - node[0]->z;
      double H1 = node[1]->v.S - node[1]->z;

      if( H0 <= 0.0 )  H0 = project->hmin;
      if( H1 <= 0.0 )  H1 = project->hmin;

      H = m[0] * H0  +  m[1] * H1;
    }

//  if( H <= 0.0 )  H = project->hmin;


    // -----------------------------------------------------------------------------------
    // compute normal vector
    // since the normal is not reduced to unit length it
    // implies the transformation of the integrand

    double nx =  dm[0]*node[0]->y + dm[1]*node[1]->y;
    double ny = -dm[0]*node[0]->x - dm[1]*node[1]->x;

    double len = sqrt( nx*nx + ny*ny );


    // weight of Gauss point j -----------------------------------------------------------

    double weight = lShape->weight[g];

    double U    = 0.0;
    double V    = 0.0;
    double Ures = 0.0;
    double cf   = 0.0;

    if( solidFlag )
    {
      U = m[0]*node[0]->v.U + m[1]*node[1]->v.U;
      V = m[0]*node[0]->v.V + m[1]*node[1]->v.V;

      Ures = sqrt( U * U  +  V * V );

      cf = m[0]*node[0]->cfw + m[1]*node[1]->cfw;
    }


    // -----------------------------------------------------------------------------------

    if( force )
    {
      double fU = weight * H * ( gravity*H/2.0 * nx );
      double fV = weight * H * ( gravity*H/2.0 * ny );


      // wall roughness ------------------------------------------------------------------

      if( solidFlag )
      {
        fU += weight * len * H * cf * U * Ures;
        fV += weight * len * H * cf * V * Ures;
      }


      // add fU, fV to force vector ------------------------------------------------------

      force[r0U] -= m[0] * fU;        force[r0V] -= m[0] * fV;
      force[r1U] -= m[1] * fU;        force[r1V] -= m[1] * fV;
    }


    // compute estifm, if no flow depth specified ----------------------------------------

    if( estifm )
    {
      double dfUdU = 0.0;
      double dfUdV = 0.0;
      double dfUdH = 0.0;
      double dfVdU = 0.0;
      double dfVdV = 0.0;
      double dfVdH = 0.0;

      if( !isFS(node[2]->bc.kind, BCON::kOutlet)  )
      {
        dfUdH += weight * ( gravity*H * nx );
        dfVdH += weight * ( gravity*H * ny );
      }

      if( solidFlag )
      {
        dfUdH += weight * len * cf * U * Ures;
        dfVdH += weight * len * cf * V * Ures;

        double iUres;
        if ( Ures > 1.0e-8 )  iUres = 1.0 / Ures;
        else                  iUres = 0.0;

        dfUdU += weight * len * cf * H * (U * U * iUres  +  Ures);
        dfUdV += weight * len * cf * H * (U * V * iUres);

        dfVdU += weight * len * cf * H * (U * V * iUres);
        dfVdV += weight * len * cf * H * (V * V * iUres  +  Ures);
      }

      estifm[r0U][r0U] += m[0] * dfUdU * m[0];
      estifm[r0U][r1U] += m[0] * dfUdU * m[1];

      estifm[r1U][r0U] += m[1] * dfUdU * m[0];
      estifm[r1U][r1U] += m[1] * dfUdU * m[1];

      estifm[r0U][r0V] += m[0] * dfUdV * m[0];
      estifm[r0U][r1V] += m[0] * dfUdV * m[1];

      estifm[r1U][r0V] += m[1] * dfUdV * m[0];
      estifm[r1U][r1V] += m[1] * dfUdV * m[1];

      estifm[r0V][r0U] += m[0] * dfVdU * m[0];
      estifm[r0V][r1U] += m[0] * dfVdU * m[1];

      estifm[r1V][r0U] += m[1] * dfVdU * m[0];
      estifm[r1V][r1U] += m[1] * dfVdU * m[1];

      estifm[r0V][r0V] += m[0] * dfVdV * m[0];
      estifm[r0V][r1V] += m[0] * dfVdV * m[1];

      estifm[r1V][r0V] += m[1] * dfVdV * m[0];
      estifm[r1V][r1V] += m[1] * dfVdV * m[1];

      estifm[r0U][c0H] += m[0] * dfUdH * m[0];
      estifm[r0V][c0H] += m[0] * dfVdH * m[0];

      estifm[r1U][c0H] += m[1] * dfUdH * m[0];
      estifm[r1V][c0H] += m[1] * dfVdH * m[0];

      estifm[r0U][c1H] += m[0] * dfUdH * m[1];
      estifm[r0V][c1H] += m[0] * dfVdH * m[1];

      estifm[r1U][c1H] += m[1] * dfUdH * m[1];
      estifm[r1V][c1H] += m[1] * dfVdH * m[1];
    }
  }


  // apply transformation ----------------------------------------------------------------

  Rotate2D( nnd, elem->nd, 3, estifm, force );


  // -------------------------------------------------------------------------------------
  // In case of experimental boundary forces on inlet set equation row to zero.
  // Accounting for boundary forces is done in method EQS_UVS2D::Region().

  for( int i=0; i<ncn; i++ )
  {
    BCON* bcon = &elem->nd[i]->bc;

    if( isFS(bcon->kind, BCON::kInlet) )
    {
      // remove  equation row dfU and force vector ---------------------------------------

      if( estifm )
      {
        for( int j=0; j<maxEleq; j++ )
        {
          estifm[i][j] = 0.0;
        }
      }

      if( force )  force[i] = 0.0;
    }
  }
}


void EQS_UVS2D_ME::Region( ELEM     *elem,
                           PROJECT  *project,
                           double  **estifm,
                           double   *force )
{
  // -------------------------------------------------------------------------------------
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

  // index to bubble nodes (bid) and equations (uid, vis, sid)
  int *bid;

  if( elem->shape == kTri ) bid = tri_id;
  else                      bid = quad_id;

  if( isFS(elem->nd[ncn]->flag, NODE::kRotat) )  CF( elem->nd[ncn]->flag, NODE::kRotat );

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


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // use GAUSS point integration to solve momentum equations for x- and
  // y-direction (U and V) and continuity equation (H)

  double area = 0.0;

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with quadratic shape functions

    double* dfdxPtr = lShape->dfdx[g];
    double* dfdyPtr = lShape->dfdy[g];

    double  trafo[2][2];

    double detj   = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
    double weight = detj * lShape->weight[g];

    area += weight;


    // -----------------------------------------------------------------------------------
    // compute values of linear shape functions at GP g

    double* m = lShape->f[g];

    double dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    for( int i=0; i<ncn; i++ )
    {
      dmdx[i] = trafo[0][0] * dfdxPtr[i] + trafo[0][1] * dfdyPtr[i];
      dmdy[i] = trafo[1][0] * dfdxPtr[i] + trafo[1][1] * dfdyPtr[i];
    }


    // -----------------------------------------------------------------------------------
    // compute values of bubble shape functions at GP g

    double *b = bShape->f[g];

    double  dbdx[kMaxNodes2D], dbdy[kMaxNodes2D];

    for( int i=0; i<nbn; i++ )
    {
      int k = bid[i];

      double dfdx = bShape->dfdx[g][k];
      double dfdy = bShape->dfdy[g][k];

      dbdx[k] = trafo[0][0] * dfdx + trafo[0][1] * dfdy;
      dbdy[k] = trafo[1][0] * dfdx + trafo[1][1] * dfdy;
    }


    // ------------------------------------------------------------------------------------
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

    for( int i=0; i<ncn; i++ )
    {
      NODE *node = elem->nd[i];
      BCON *bcon = &node->bc;

      double ndZ = node->z;
      double ndH = node->v.S - ndZ;

      if( ndH <= 0.0 )  ndH = project->hmin;

      double ndSS;

      // -----------------------------------------------------------------------------------
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
      // ---------------------------------------------------------------------------------

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
    }

//  if( H <= 0.0 )  H = project->hmin;


    // -----------------------------------------------------------------------------------

    double U    = 0.0;
    double dUdt = 0.0;
    double dUdx = 0.0;
    double dUdy = 0.0;

    double V    = 0.0;
    double dVdt = 0.0;
    double dVdx = 0.0;
    double dVdy = 0.0;

    for( int i=0; i<nbn; i++ )
    {
      int k = bid[i];

      NODE*  node = elem->nd[k];

      double ndU  = node->v.U;
      double ndV  = node->v.V;

      U    +=    b[k] * ndU;
      dUdt +=    b[k] * node->v.dUdt;
      dUdx += dbdx[k] * ndU;
      dUdy += dbdy[k] * ndU;

      V    +=    b[k] * ndV;
      dVdt +=    b[k] * node->v.dVdt;
      dVdx += dbdx[k] * ndV;
      dVdy += dbdy[k] * ndV;
    }

    // compute dispersion coefficients ---------------------------------------------------

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


    // compute eddy viscosity ------------------------------------------------------------

    double vt = 0.0;

    if( isFS(project->actualTurb, BCONSET::kVtConstant) )
    {
      vt = type->vt;
    }
    else
    {
      for( int i=0; i<ncn; i++ )  vt += m[i] * elem->nd[i]->vt;
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < type->vt )  vt = type->vt;
    }


    // add kinematic viscosity to eddy viscosity  ----------------------------------------

    vt += project->vk;


    // -----------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double Ures = sqrt( U*U + V*V );           // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double *forcePtr;


      // ---------------------------------------------------------------------------------
      // compute x-momentum

      f   = H * dUdt;                                  // time
      f  += H * (U * dUdx  +  V * dUdy);               // convection

      fx  = H * vt * (dUdx + dUdx);                    // eddy viscosity
      fy  = H * vt * (dUdy + dVdx);

      fx -= H * uu;                                    // turbulence
      fy -= H * uv;

      f  -= H * gravity * dadx;                        // gravity
      fx -= H * H * gravity / 2.0;

      f  += cf * Ures * U;                             // bottom friction

      fx +=  H * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );  // dispersion
      fy +=  H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];
        forcePtr[i] -= b[k] * f  +  dbdx[k] * fx  +  dbdy[k] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute y-momentum

      f   = H * dVdt;                                  // time
      f  += H * (U * dVdx  +  V * dVdy);               // convection

      fx  = H * vt * (dVdx + dUdy);                    // eddy viscosity
      fy  = H * vt * (dVdy + dVdy);

      fx -= H * uv;                                    // turbulence
      fy -= H * vv;

      f  +=  H * gravity * dady;                       // gravity
      fy -=  H * H * gravity / 2.0;

      f  += cf * Ures * V;                             // bottom friction

      fx +=  H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );    // dispersion
      fy +=  H * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + eqidV;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];
        forcePtr[i] -= b[k] * f  +  dbdx[k] * fx  +  dbdy[k] * fy;
      }


      // ---------------------------------------------------------------------------------
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


      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      double iUres;

      if( Ures > 1.0e-8 ) iUres = 1.0 / Ures;
      else                iUres = 0.0;


      // --- U-derivative of x-momentum --------------------------------------------------

      df__  =  weight * H * relaxThdt_UV;

      df__ +=  weight * H * dUdx;
      df_x  =  weight * H * U;
      df_y  =  weight * H * V;

      dfxx  =  weight * H * vt * 2.0;
      dfyy  =  weight * H * vt;

      df__ +=  weight * cf * (iUres * U*U  +  Ures);

      dfx_  =  weight * H * ( 2.0*U*Dxx - 2.0*V*Dxy );
      dfy_  =  weight * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        t[i]  = df__ * b[k]  +  df_x * dbdx[k]  +  df_y * dbdy[k];
        tx[i] = dfx_ * b[k]  +  dfxx * dbdx[k];
        ty[i] = dfy_ * b[k]                     +  dfyy * dbdy[k];
      }

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        double *estifmPtr = estifm[i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[k]*t[j] + dbdx[k]*tx[j] + dbdy[k]*ty[j];
        }
      }


      // --- V-derivative of x-momentum --------------------------------------------------

      df__  =  weight * H * dUdy;

      dfyx  =  weight * H * vt;

      df__ +=  weight * cf * iUres * U * V;

      dfx_  =  weight * H * ( 2.0*V*Dyy - 2.0*V*Dxy );
      dfy_  =  weight * H * ( V*(Dxx-Dyy) - 2.0*V*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        t[i]  = df__ * b[k];
        tx[i] = dfx_ * b[k];
        ty[i] = dfy_ * b[k]  +  dfyx * dbdx[k];
      }

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        double *estifmPtr = estifm[i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[eqidV + j] += b[k]*t[j] + dbdx[k]*tx[j] + dbdy[k]*ty[j];
        }
      }


      // --- H-derivative of x-momentum --------------------------------------------------

      df__  =  weight * dUdt;
      df__ +=  weight * (U * dUdx + V * dUdy);

      dfx_  =  weight * vt * (dUdx + dUdx);
      dfy_  =  weight * vt * (dUdy + dVdx);

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
        int k = bid[i];

        double *estifmPtr = estifm[i];

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[eqidS + j] += b[k]*t[j] + dbdx[k]*tx[j] + dbdy[k]*ty[j];
        }
      }


      // --- U-derivative of y-momentum --------------------------------------------------

      df__  =  weight * H * dVdx;

      dfxy  =  weight * H * vt;

      df__ +=  weight * cf * iUres * U * V;

      dfx_  =  weight * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );
      dfy_  =  weight * H * ( 2.0*V*Dxy + 2.0*U*Dyy );

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        t[i]  = df__ * b[k];
        tx[i] = dfx_ * b[k]  +  dfxy * dbdy[k];
        ty[i] = dfy_ * b[k];
      }

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        double *estifmPtr = estifm[eqidV + i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[k]*t[j] + dbdx[k]*tx[j] + dbdy[k]*ty[j];
        }
      }


      // --- V-derivative of y-momentum --------------------------------------------------

      df__  =  weight * H * relaxThdt_UV;

      df__ +=  weight * H * dVdy;
      df_x  =  weight * H * U;
      df_y  =  weight * H * V;

      dfxx  =  weight * H * vt;
      dfyy  =  weight * H * vt * 2.0;

      df__ +=  weight * cf * (iUres * V*V  +  Ures);

      dfx_  =  weight * H * ( U*(Dxx-Dyy) - 2.0*V*Dxy );
      dfy_  =  weight * H * ( 2.0*V*Dxx + 2.0*U*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        t[i]  = df__ * b[k]  +  df_x * dbdx[k]  +  df_y * dbdy[k];
        tx[i] = dfx_ * b[k]  +  dfxx * dbdx[k];
        ty[i] = dfy_ * b[k]  +  dfyy * dbdy[k];
      }

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        double *estifmPtr = estifm[eqidV + i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[eqidV + j] += b[k]*t[j] + dbdx[k]*tx[j] + dbdy[k]*ty[j];
        }
      }


      // H-derivative of y-momentum ------------------------------------------------------

      df__  =  weight * dVdt;
      df__ +=  weight * (U * dVdx + V * dVdy);

      dfx_  =  weight * vt * (dVdx + dUdy);
      dfy_  =  weight * vt * (dVdy + dVdy);

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
        int k = bid[i];

        double *estifmPtr = estifm[eqidV + i];

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[eqidS + j] += b[k]*t[j] + dbdx[k]*tx[j] + dbdy[k]*ty[j];
        }
      }


      // U-derivative of continuity ------------------------------------------------------

      df__ = weight * dHdx;
      df_x = weight * H;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        t[i] = df__ * b[k]  +  df_x * dbdx[k];
      }

      for( int i=0; i<ncn; i++ )
      {
        double *estifmPtr = estifm[eqidS + i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += m[i]*t[j];
        }
      }


      // V-derivative of continuity ------------------------------------------------------

      df__ = weight * dHdy;
      df_y = weight * H;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        t[i] = df__ * b[k] + df_y * dbdy[k];
      }

      for( int i=0; i<ncn; i++ )
      {
        double *estifmPtr = estifm[eqidS + i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[eqidV + j] += m[i]*t[j];
        }
      }


      // H-derivative of continuity ------------------------------------------------------

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
        double *estifmPtr = estifm[eqidS + i];

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[eqidS + j] += m[i]*t[j];
        }
      }
    }
  } // END of loop over all GAUSS points


  // -------------------------------------------------------------------------------------
  // apply transformation

  Rotate2D( nnd, elem->nd, 3, estifm, force );


  // -------------------------------------------------------------------------------------
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
        H     = project->hmin;
        specQ = 0.0;
      }

      // force vector
      if( force )  force[i] = area * (specQ - Un*H);

      // stiffness matrix
      if( estifm )
      {
        estifm[i][i]      = area * H;
        estifm[i][sid[i]] = area * Un;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // eliminate momentum equation V for element centers (partial Gauss elimination)

  double* VPtr = estifm[vid[ncn]];

  for( int i=0; i<nbn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in x-momentum equation

    MPtr = estifm[uid[i]];

    factor = MPtr[vid[ncn]] / VPtr[vid[ncn]];

    force[uid[i]] -= factor * force[vid[ncn]];

    for( int j=0; j<nbn; j++ )
    {
      MPtr[uid[j]] -= factor * VPtr[uid[j]];
      MPtr[vid[j]] -= factor * VPtr[vid[j]];
    }

    MPtr[vid[ncn]] = 0.0;
  }

  for( int i=0; i<nbn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in y-momentum equation

    MPtr = estifm[vid[i]];

    factor = MPtr[vid[ncn]] / VPtr[vid[ncn]];

    force[vid[i]] -= factor * force[vid[ncn]];

    for( int j=0; j<nbn; j++ )
    {
      MPtr[uid[j]] -= factor * VPtr[uid[j]];
      MPtr[vid[j]] -= factor * VPtr[vid[j]];
    }

    MPtr[vid[ncn]] = 0.0;
  }


  // save elimination equation -----------------------------------------------------------

  double* VElimEqPtr = VElimEq[elem->Getno()];

  for( int i=0; i<nbn; i++ )
  {
    VElimEqPtr[i]     = VPtr[uid[i]];
    VElimEqPtr[i+nbn] = VPtr[vid[i]];
  }

  VElimEqPtr[2*nbn]   = VPtr[vid[ncn]];   // diagonal entry
  VElimEqPtr[2*nbn+1] = force[vid[ncn]];  // RHS


  // -------------------------------------------------------------------------------------
  // eliminate momentum equation U for element centers (partial Gauss elimination)

  double* UPtr = estifm[uid[ncn]];

  for( int i=0; i<nbn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in x-momentum equation

    MPtr = estifm[uid[i]];

    factor = MPtr[uid[ncn]] / UPtr[uid[ncn]];

    force[uid[i]] -= factor * force[uid[ncn]];

    for( int j=0; j<ncn; j++ )
    {
      MPtr[uid[j]] -= factor * UPtr[uid[j]];
      MPtr[vid[j]] -= factor * UPtr[vid[j]];
    }

    MPtr[uid[3]] = 0.0;


    // eliminate in y-momentum equation

    MPtr = estifm[vid[i]];

    factor = MPtr[uid[ncn]] / UPtr[uid[ncn]];

    force[vid[i]] -= factor * force[uid[ncn]];

    for( int j=0; j<ncn; j++ )
    {
      MPtr[uid[j]] -= factor * UPtr[uid[j]];
      MPtr[vid[j]] -= factor * UPtr[vid[j]];
    }

    MPtr[uid[ncn]] = 0.0;
  }


  // save elimination equation -----------------------------------------------------------

  double* UElimEqPtr = UElimEq[elem->Getno()];

  for( int i=0; i<nbn; i++ )
  {
    UElimEqPtr[i]     = UPtr[uid[i]];
    UElimEqPtr[i+nbn] = UPtr[vid[i]];
  }

  UElimEqPtr[2*nbn]   = UPtr[uid[ncn]];   // diagonal entry
  UElimEqPtr[2*nbn+1] = force[uid[ncn]];  // RHS
}
