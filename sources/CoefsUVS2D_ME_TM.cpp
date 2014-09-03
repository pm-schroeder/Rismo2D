// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_ME_TM
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

#include "EqsUVS2D_ME_TM.h"


////////////////////////////////////////////////////////////////////////////////////////////////////

int EQS_UVS2D_ME_TM::Coefs( ELEM*     elem,
                            PROJECT  *project,
                            double  **estifm,
                            double   *force )
{
  if( isFS(elem->flag, ELEM::kDry) ) return 0;

  if( isFS(elem->flag, ELEM::kBound) )
  {
    Bound( elem, project, estifm, force );
  }
  else
  {
    Region( elem, project, estifm, force );
  }

  return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_ME_TM::Bound( ELEM*    elem,
                             PROJECT* project,
                             double** estifm,
                             double*  force )
{
  SHAPE* lShape  = elem->GetLShape();

  int    ncn     = lShape->nnd;
  int    nnd     = elem->GetQShape()->nnd;

//int eqidU = 0;
  int eqidV = nnd;
  int eqidS = 2 * nnd;

  double gravity = project->g;

  // time weighting of pressure/gravity "thp" and forces "thf"
  double thp = project->timeint.thetaFlow;
  double thf = project->timeint.thetaFlow;

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
  node[2] = elem->nd[2];      // midside node


  // row index for force vector and element stiffness matrix ---------------------------------------

  int r0U = 0;       int r1U = 1;
  int r0V = eqidV;   int r1V = eqidV + 1;


  // column index for element stiffness matrix -----------------------------------------------------

  int c0H = eqidS;   int c1H = eqidS + 1;


  for( int g=0; g<lShape->ngp; g++ )   // loop on GAUSS points
  {
    double *m  = lShape->f[g];         // linear shape
    double *dm = lShape->dfdx[g];


    // test midside node for boundary condition outFlow --------------------------------------------

    double H;

    if( isFS(node[2]->bc.kind, BCON::kOutlet) )
    {
      // use experimental downstream flow depth ... ------------------------------------------------
      double H0 = node[0]->bc.val->S - node[0]->z;
      double H1 = node[1]->bc.val->S - node[1]->z;

      if( H0 <= 0.0 ) H0 = project->hmin;
      if( H1 <= 0.0 ) H1 = project->hmin;

      H = m[0] * H0  +  m[1] * H1;
    }
    else
    {
     // ... or compute flow depth at gauss point ---------------------------------------------------
      double H0 = node[0]->v.S - node[0]->z;
      double H1 = node[1]->v.S - node[1]->z;

      if( H0 <= 0.0 ) H0 = project->hmin;
      if( H1 <= 0.0 ) H1 = project->hmin;

      H = m[0] * H0  +  m[1] * H1;
    }

//  if( H <= 0.0 )  H = project->hmin;


    // ---------------------------------------------------------------------------------------------
    // compute normal vector
    // since the normal is not reduced to unit length it
    // implies the transformation of the integrand

    double nx =  dm[0]*node[0]->y + dm[1]*node[1]->y;
    double ny = -dm[0]*node[0]->x - dm[1]*node[1]->x;

    double len = sqrt( nx*nx + ny*ny );


    // weight of Gauss point j ---------------------------------------------------------------------

    double weight = lShape->weight[g];

    double U    = 0.0;
    double V    = 0.0;
    double Ures = 0.0;
    double cfw  = 0.0;

    if( solidFlag )
    {
      U = m[0]*node[0]->v.U + m[1]*node[1]->v.U;
      V = m[0]*node[0]->v.V + m[1]*node[1]->v.V;

      Ures = sqrt( U * U  +  V * V );

      cfw = m[0]*node[0]->cfw + m[1]*node[1]->cfw;
    }


    // ---------------------------------------------------------------------------------------------

    if( force )
    {
      double fU = weight * thp * H * ( gravity*H/2.0 * nx );
      double fV = weight * thp * H * ( gravity*H/2.0 * ny );


      // wall roughness ----------------------------------------------------------------------------

      if( solidFlag )
      {
        fU += weight * thf * len * H * cfw * U * Ures;
        fV += weight * thf * len * H * cfw * V * Ures;
      }


      // add fU, fV to force vector ----------------------------------------------------------------

      force[r0U] -= m[0] * fU;        force[r0V] -= m[0] * fV;
      force[r1U] -= m[1] * fU;        force[r1V] -= m[1] * fV;
    }


    // compute estifm, if no flow depth specified --------------------------------------------------

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
        dfUdH += weight * thp * ( gravity*H * nx );
        dfVdH += weight * thp * ( gravity*H * ny );
      }

      if( solidFlag )
      {
        dfUdH += weight * thf * len * cfw * U * Ures;
        dfVdH += weight * thf * len * cfw * V * Ures;

        double iUres;
        if ( Ures > 1.0e-8 )  iUres = 1.0 / Ures;
        else                  iUres = 0.0;

        dfUdU += weight * thf * len * cfw * H * (U * U * iUres  +  Ures);
        dfUdV += weight * thf * len * cfw * H * (U * V * iUres);

        dfVdU += weight * thf * len * cfw * H * (U * V * iUres);
        dfVdV += weight * thf * len * cfw * H * (V * V * iUres  +  Ures);
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


  // apply transformation --------------------------------------------------------------------------

  Rotate2D( ncn, eqidV, elem->nd, estifm, force );


  // -----------------------------------------------------------------------------------------------
  // In case of experimental boundary forces on inlet set equation row to zero.
  // Accounting for boundary forces is done in method EQS_UVS2D_TM::Region().

  for( int i=0; i<nnd; i++ )
  {
    BCON* bcon = &elem->nd[i]->bc;

    if( isFS(bcon->kind, BCON::kInlet) )
    {
      // remove  equation row dfU and force vector -------------------------------------------------

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

////////////////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_ME_TM::Region( ELEM*    elem,
                              PROJECT* project,
                              double** estifm,
                              double*  force )
{
  // -----------------------------------------------------------------------------------------------
  // initializations

  SHAPE* lShape = elem->GetLShape();
  SHAPE* bShape = elem->GetBShape();

  TYPE* type = TYPE::Getid( elem->type );

  int ngp = lShape->ngp;                // number of GAUSS points

  int nnd = elem->GetQShape()->nnd;     // number of nodes
  int ncn = lShape->nnd;                // number of corner nodes
  int nbn = bShape->nnd;                // number of nodes for bubble shape function

  double gravity = project->g;

//int eqidU = 0;
  int eqidV  = nnd;
  int eqidS  = 2 * nnd;

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

  double dt  = project->timeint.incTime.Getsec();
  double th  = project->timeint.thetaFlow;
  double tha = project->timeint.thetaFlow;
  double thd = project->timeint.thetaFlow;
  double thp = project->timeint.thetaFlow;
  double thf = project->timeint.thetaFlow;
  double thc = project->timeint.thetaFlow;

  // -----------------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<ncn; i++ )
  {
    NODE *nd = elem->nd[i];

    x[i] = nd->x - x[0];
    y[i] = nd->y - y[0];
  }

  x[0] = y[0] = 0.0;


  // -----------------------------------------------------------------------------------------------
  // use GAUSS point integration to solve momentum equations for x- and
  // y-direction (U and V) and continuity equation (H)

  double area = 0.0;

  for( int g=0; g<ngp; g++ ) // START of loop over all GAUSS point
  {
    // ---------------------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with quadratic shape functions

    double* dfdxPtr = lShape->dfdx[g];
    double* dfdyPtr = lShape->dfdy[g];

    double trafo[2][2];

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


    // ---------------------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP g
    //             horizontal velocities: U and V
    //             water surface        : S
    //             flow depth           : H
    //             Source or Sink       : SS
    //             bottom elevation     : a
    //             eddy viscosity       : vt, cf
    //             Reynolds stresses    : uu, uv, and vv

    // integrate H, a and cf with linear shape

    double S     = 0.0;
    double So    = 0.0;
    double dSodt = 0.0;

    double H     = 0.0;
    double dHdx  = 0.0;
    double dHdy  = 0.0;

    double dadx  = 0.0;
    double dady  = 0.0;

    double cf    = 0.0;

    double SS    = 0.0;

    double uu    = 0.0;
    double uv    = 0.0;
    double vv    = 0.0;

    double dUodt = 0.0;
    double dVodt = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      NODE* node = elem->nd[i];
      BCON* bcon = &node->bc;

      double ndS = node->v.S;
      double ndZ = node->z;
      double ndH = ndS - ndZ;

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

      S     +=    m[i] * ndS;
      So    +=    m[i] * node->vo.S;
      dSodt +=    m[i] * node->vo.dSdt;

      H     +=    m[i] * ndH;
      dHdx  += dmdx[i] * ndH;
      dHdy  += dmdy[i] * ndH;

      dadx  += dmdx[i] * ndZ;
      dady  += dmdy[i] * ndZ;

      cf    +=    m[i] * node->cf;

      SS    +=    m[i] * ndSS;            // Source or Sink

      uu    +=    m[i] * node->uu;
      uv    +=    m[i] * node->uv;
      vv    +=    m[i] * node->vv;

      dUodt +=    m[i] * node->vo.dUdt;
      dVodt +=    m[i] * node->vo.dVdt;
    }


    // integrate U, V, uu, uv and vv with quadratic shape

    double U     = 0.0;
    double Uo    = 0.0;
    double dUdx  = 0.0;
    double dUdy  = 0.0;

    double V     = 0.0;
    double Vo    = 0.0;
    double dVdx  = 0.0;
    double dVdy  = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      NODE* node = elem->nd[i];

      double ndU = node->v.U;
      double ndV = node->v.V;

      U     +=    b[i] * ndU;
      Uo    +=    b[i] * node->vo.U;
      dUdx  += dbdx[i] * ndU;
      dUdy  += dbdy[i] * ndU;

      V     +=    b[i] * ndV;
      Vo    +=    b[i] * node->vo.V;
      dVdx  += dbdx[i] * ndV;
      dVdy  += dbdy[i] * ndV;
    }

    U    +=    b[ncn] * elem->U;
    Uo   +=    b[ncn] * elem->Uo;
    dUdx += dbdx[ncn] * elem->U;
    dUdy += dbdy[ncn] * elem->U;

    V    +=    b[ncn] * elem->V;
    Vo   +=    b[ncn] * elem->Vo;
    dVdx += dbdx[ncn] * elem->V;
    dVdy += dbdy[ncn] * elem->V;


    // compute dispersion coefficients -------------------------------------------------------------

    double Dxx = 0.0;
    double Dxy = 0.0;
    double Dyy = 0.0;

    if( project->actualDisp > 0 )
    {
      for( int j=0; j<ncn; j++ )
      {
        NODE *node = elem->nd[j];
        Dxx += m[j] * node->Dxx;
        Dxy += m[j] * node->Dxy;
        Dyy += m[j] * node->Dyy;
      }
    }


    // compute eddy viscosity ----------------------------------------------------------------------

    double vt = 0.0;

    if( isFS(project->actualTurb, BCONSET::kVtConstant) )
    {
      vt = type->vt;
    }
    else
    {
      for( int j=0; j<ncn; j++ )  vt += m[j] * elem->nd[j]->vt;
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < type->vt )  vt = type->vt;
    }


    // add kinematic viscosity to eddy viscosity  ----------------------------------------

    vt += project->vk;


    // ---------------------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double  Ures = sqrt( U*U + V*V );          // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // -------------------------------------------------------------------------------------------
      // compute x-momentum equation

      f   =  H * (U - Uo) / dt;                             // time
      f  -=  H * (1.0 - th) * dUodt;

      f  +=  tha * H * (U * dUdx  +  V * dUdy);             // convection

      fx  =  thd * H * vt * (dUdx + dUdx);                  // eddy viscosity
      fy  =  thd * H * vt * (dUdy + dVdx);

      fx -=  thd * H * uu;                                  // turbulence
      fy -=  thd * H * uv;

      f  +=  thp * H * gravity * dadx;                      // gravity
      fx -=  thp * H * H * gravity / 2.0;

      f  +=  thf * cf * Ures * U;                           // bottom friction

      fx +=  thd * H * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy ); // dispersion
      fy +=  thd * H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int i=0; i<nbn; i++ )
      {
        forcePtr[i] -= b[i] * f  +  dbdx[i] * fx  +  dbdy[i] * fy;
      }


      // -------------------------------------------------------------------------------------------
      // compute y-momentum equation

      f   =  H * (V - Vo) / dt;                             // time
      f  -=  H * (1.0 - th) * dVodt;

      f  +=  tha * H * (U * dVdx  +  V * dVdy);             // convection

      fx  =  thd * H * vt * (dVdx + dUdy);                  // eddy viscosity
      fy  =  thd * H * vt * (dVdy + dVdy);

      fx -=  thd * H * uv;                                  // turbulence
      fy -=  thd * H * vv;

      f  +=  thp * H * gravity * dady;                      // gravity
      fy -=  thp * H * H * gravity / 2.0;

      f  +=  thf * cf * Ures * V;                           // bottom friction

      fx +=  thd * H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );   // dispersion
      fy +=  thd * H * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + eqidV;

      for( int i=0; i<nbn; i++ )
      {
        forcePtr[i] -= b[i] * f  +  dbdx[i] * fx  +  dbdy[i] * fy;
      }


      // -------------------------------------------------------------------------------------------
      // compute continuity equation (solve only for corner nodes)

      f  = (S - So) / dt;
      f -= (1.0 - th) * dSodt;
      f += thc * (H * (dUdx + dVdy) + U * dHdx + V * dHdy);
      f += SS;                                              // Source or Sink

      f *= weight;

      forcePtr = force + eqidS;

      for( int i=0; i<ncn; i++ )
      {
        forcePtr[i] -= m[i] * f;
      }
    }


    if( estifm )
    {
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfx_, dfxx, dfxy, dfy_, dfyx, dfyy;
      double* estifmPtr;

      // -------------------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      double iUres;

      if( Ures > 1.0e-9 ) iUres = 1.0 / Ures;
      else                iUres = 0.0;


      // U-derivative of x-momentum ----------------------------------------------------------------

      df__  =  weight * H / dt;
      df__ +=  weight * tha * H * dUdx;
      df_x  =  weight * tha * H * U;
      df_y  =  weight * tha * H * V;

      dfxx  =  weight * thd * H * vt * 2.0;
      dfyy  =  weight * thd * H * vt;

      df__ +=  weight * thf * cf * (iUres * U*U  +  Ures);

      dfx_  =  weight * thd * H * ( 2.0*U*Dxx - 2.0*V*Dxy );
      dfy_  =  weight * thd * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        t[i]  = df__ * b[i]  +  df_x * dbdx[i]  +  df_y * dbdy[i];
        tx[i] = dfx_ * b[i]  +  dfxx * dbdx[i];
        ty[i] = dfy_ * b[i]                     +  dfyy * dbdy[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        estifmPtr = estifm[i];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // V-derivative of x-momentum ----------------------------------------------------------------

      df__  =  weight * tha * H * dUdy;

      dfyx  =  weight * thd * H * vt;

      df__ +=  weight * thf * cf * iUres * U * V;

      dfx_  =  weight * thd * H * ( 2.0*V*Dyy - 2.0*V*Dxy );
      dfy_  =  weight * thd * H * ( V*(Dxx-Dyy) - 2.0*V*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        t[i]  = df__ * b[i];
        tx[i] = dfx_ * b[i];
        ty[i] = dfy_ * b[i]  +  dfyx * dbdx[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        estifmPtr = estifm[i] + eqidV;

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // H-derivative of x-momentum ----------------------------------------------------------------

      df__  =  weight * (U - Uo) / dt;

      df__ +=  weight * tha * (U * dUdx + V * dUdy);

      dfx_  =  weight * thd * vt * (dUdx + dUdx);
      dfy_  =  weight * thd * vt * (dUdy + dVdx);

      dfx_ -=  weight * thd * uu;
      dfy_ -=  weight * thd * uv;

      df__ +=  weight * thp * gravity * dadx;
      dfx_ -=  weight * thp * gravity * H;

      dfx_ +=  weight * thd * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );
      dfy_ +=  weight * thd * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      for( int i=0; i<ncn; i++ )
      {
        t[i]  = df__ * m[i];
        tx[i] = dfx_ * m[i];
        ty[i] = dfy_ * m[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        estifmPtr = estifm[i] + eqidS;

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // U-derivative of y-momentum ----------------------------------------------------------------

      df__  =  weight * tha * H * dVdx;

      dfxy  =  weight * thd * H * vt;

      df__ +=  weight * thf * cf * iUres * U * V;

      dfx_  =  weight * thd * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );
      dfy_  =  weight * thd * H * ( 2.0*V*Dxy + 2.0*U*Dyy );

      for( int i=0; i<nbn; i++ )
      {
        t[i]  = df__ * b[i];
        tx[i] = dfx_ * b[i]  +  dfxy * dbdy[i];
        ty[i] = dfy_ * b[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        estifmPtr = estifm[i + eqidV];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // V-derivative of y-momentum ----------------------------------------------------------------

      df__  =  weight * H / dt;
      df__ +=  weight * tha * H * dVdy;
      df_x  =  weight * tha * H * U;
      df_y  =  weight * tha * H * V;

      dfxx  =  weight * thd * H * vt;
      dfyy  =  weight * thd * H * vt * 2.0;

      df__ +=  weight * thf * cf * (iUres * V*V  +  Ures);

      dfx_  =  weight * thd * H * ( U*(Dxx-Dyy) - 2.0*V*Dxy );
      dfy_  =  weight * thd * H * ( 2.0*V*Dxx + 2.0*U*Dxy );

      for( int i=0; i<nbn; i++ )
      {
        t[i]  = df__ * b[i]  +  df_x * dbdx[i]  +  df_y * dbdy[i];
        tx[i] = dfx_ * b[i]  +  dfxx * dbdx[i];
        ty[i] = dfy_ * b[i]                     +  dfyy * dbdy[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        estifmPtr = estifm[i + eqidV] + eqidV;

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // H-derivative of y-momentum ----------------------------------------------------------------

      df__  =  weight * (V - Vo) / dt;

      df__ +=  weight * tha * (U * dVdx + V * dVdy);

      dfx_  =  weight * thd * vt * (dVdx + dUdy);
      dfy_  =  weight * thd * vt * (dVdy + dVdy);

      dfx_ -=  weight * thd * uv;
      dfy_ -=  weight * thd * vv;

      df__ +=  weight * thp * gravity * dady;
      dfy_ -=  weight * thp * gravity * H;

      dfx_ +=  weight * thd * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );
      dfy_ +=  weight * thd * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      for( int i=0; i<ncn; i++ )
      {
        t[i]  = df__ * m[i];
        tx[i] = dfx_ * m[i];
        ty[i] = dfy_ * m[i];
      }

      for( int i=0; i<nbn; i++ )
      {
        estifmPtr = estifm[i + eqidV] + eqidS;

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[j] += b[i]*t[j] + dbdx[i]*tx[j] + dbdy[i]*ty[j];
        }
      }


      // U-derivative of continuity ----------------------------------------------------------------

      df__ = weight * thc * dHdx;
      df_x = weight * thc * H;

      for( int i=0; i<nbn; i++ )
      {
        t[i] = df__ * b[i]  +  df_x * dbdx[i];
      }

      for( int i=0; i<ncn; i++ )
      {
        estifmPtr = estifm[i + eqidS];

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += m[i]*t[j];
        }
      }


      // V-derivative of continuity ----------------------------------------------------------------

      df__ = weight * thc * dHdy;
      df_y = weight * thc * H;

      for( int i=0; i<nbn; i++ )
      {
        t[i] = df__ * b[i] + df_y * dbdy[i];
      }

      for( int i=0; i<ncn; i++ )
      {
        estifmPtr = estifm[i + eqidS] + eqidV;

        for( int j=0; j<nbn; j++ )
        {
          estifmPtr[j] += m[i]*t[j];
        }
      }


      // H-derivative of continuity ----------------------------------------------------------------

      df__  = weight / dt;
      df__ += weight * thc * (dUdx + dVdy);
      df_x  = weight * thc * U;
      df_y  = weight * thc * V;

      for( int i=0; i<ncn; i++ )
      {
        t[i] = df__ * m[i]  +  df_x * dmdx[i]  +  df_y * dmdy[i];
      }

      for( int i=0; i<ncn; i++ )
      {
        estifmPtr = estifm[i + eqidS] + eqidS;

        for( int j=0; j<ncn; j++ )
        {
          estifmPtr[j] += m[i]*t[j];
        }
      }
    }
  } // END of loop over all GAUSS points



  // apply transformation --------------------------------------------------------------------------

  Rotate2D( ncn, eqidV, elem->nd, estifm, force );


  // -----------------------------------------------------------------------------------------------
  // eliminate the two equation at the bubble node
  Eliminate( elem, estifm, force, ncn, nbn, eqidV, eqidS );


  // -----------------------------------------------------------------------------------------------
  // insert experimental upstream boundary forces at nodes with inflow
  // boundary condition  (qfix = constant):
  // flow depth at midside nodes is computed from values at corner nodes
  // the x-momentum equation is replaced by:
  //                     f  =  area * (qfix - Un * H)
  //                     where: Un, H are computed values

  for( int i=0; i<ncn; i++ )      // corner nodes
  {
    NODE* node = elem->nd[i];
    BCON* bcon = &node->bc;

    if( isFS(bcon->kind, BCON::kInlet) )
    {
      // set equation row dfU to zero --------------------------------------------------------------
      if( estifm )
      {
        for( int j=0; j<maxEleq; j++ )
        {
          estifm[i][j] = 0.0;
        }
      }

      // compute flow in normal direction at node i ------------------------------------------------
      double Un    = node->v.U * bcon->niox  +  node->v.V * bcon->nioy;
      double H     = node->v.S - node->z;
      double specQ = bcon->val->U;

      if( H <= 0.0 )
      {
        H = project->hmin;
        specQ = 0.0;
      }

      // force vector ------------------------------------------------------------------------------
      if( force )  force[i] = area * (specQ - Un*H);

      // jacobi matrix -----------------------------------------------------------------------------
      if( estifm )
      {
        estifm[i][i]         = area * H;
        estifm[i][i + eqidS] = area * Un;
      }
    }
  }
}
