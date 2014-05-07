// ======================================================================================
//
// Copyright (C) 1992-2013  by  P.M. SCHROEDER
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
// compute element stiffness matrix with anisotrop turbulence (Elder model)
// solve for momentum equations (U,V) and continuity equation (S)
//
// methods                     description
// -------------------------   -------------------------------------------------------
// EQS_UVS2D::Coefs()
// EQS_UVS2D::Region()         two-dimensional region elements
//
// EQS_UVS2D::Region_40306()   two-dimensional region elements (Rismo version 4.03.06)
// ======================================================================================
//
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 18.08.2012     sc     Anisotrop formulation of the method EQS_UVS2D_AI::Region().
//                       Class EqsUVS2D_AI derived from class EqsUVS2D.
// ======================================================================================

#include "Defs.h"
#include "Report.h"
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsUVS2D_AI.h"


// ======================================================================================

int EQS_UVS2D_AI::Coefs( ELEM*    elem,
                         PROJECT* project,
                         double** estifm,
                         double*  force )
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


// ======================================================================================
/*
void EQS_UVS2D_AI::Region( ELEM*    elem,
                           PROJECT* project,
                           double** estifm,
                           double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int ngp = qShape->ngp;         // number of GAUSS points
  int nnd = qShape->nnd;         // number of nodes in all
  int ncn = lShape->nnd;         // number of corner nodes

  TYPE* type = TYPE::Getid( elem->type );

  int startV = nnd;
  int startS = 2 * nnd;

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

  double gravity = project->g;

  double dt  = project->timeint.incTime.Getsec();
  double th  = project->timeint.thetaFlow;
  double tha = project->timeint.thetaFlow;
  double thd = project->timeint.thetaFlow;
  double thp = project->timeint.thetaFlow;
  double thf = project->timeint.thetaFlow;
  double thc = project->timeint.thetaFlow;

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

  for( int g=0; g<ngp; g++ ) // START of loop over all GAUSS point
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with quadratic shape functions

    double* dfdxPtr = qShape->dfdx[g];
    double* dfdyPtr = qShape->dfdy[g];

    double trafo[2][2];

    double detj = qShape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );

    double weight = detj * qShape->weight[g];

    area += weight;


    // -----------------------------------------------------------------------------------
    // compute values of quadratic shape functions at GP g

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = qShape->f[g];

    for( int j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // ----------------------------------------------------------------------------------
    // compute values of linear shape functions at GP g

    dfdxPtr = lShape->dfdx[g];
    dfdyPtr = lShape->dfdy[g];

    double dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    double* m = lShape->f[g];

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // ------------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP g
    //             horizontal velocities: U and V
    //             water surface        : S
    //             flow depth           : H
    //             Source or Sink       : SS
    //             bottom elevation     : a
    //             eddy viscosity       : vt, cf
    //             Reynolds stresses    : uu, uv, and vv

    // integrate S, SoH, a and cf with linear shape

    double S    = 0.0;
    double So   = 0.0;

    double H    = 0.0;
    double Ho   = 0.0;
    double dHdt = 0.0;
    double dHdx = 0.0;
    double dHdy = 0.0;

    double dadx = 0.0;
    double dady = 0.0;

    double cf   = 0.0;

    double SS   = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE* node = elem->nd[j];
      BCON* bcon = &node->bc;

      double ndZ = node->z;

      double ndS = node->v.S;
      double ndH = ndS - ndZ;
      if( ndH <= 0.0 )  ndH = project->hmin;

      double ndSo = node->vo.S;
      double ndHo = ndSo - ndZ;
      if( ndHo <= 0.0 )  ndHo = project->hmin;

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
      // -----------------------------------------------------------------------------------

      S     +=    m[j] * ndS;
      So    +=    m[j] * ndSo;

      H     +=    m[j] * ndH;
      Ho    +=    m[j] * ndHo;
      dHodt +=    m[j] * node->vo.dSdt;
      dHdx  += dmdx[j] * ndH;
      dHdy  += dmdy[j] * ndH;

      dadx  += dmdx[j] * ndZ;
      dady  += dmdy[j] * ndZ;

      cf    +=    m[j] * node->cf;

      SS    +=    m[j] * ndSS;            // Source or Sink
    }

//  if( H <= 0.0 )  H = project->hmin;


    // integrate U, V, uu, uv and vv with quadratic shape

    double U     = 0.0;
    double Uo    = 0.0;
    double dUodt = 0.0;
    double dUdx  = 0.0;
    double dUdy  = 0.0;

    double V     = 0.0;
    double Vo    = 0.0;
    double dVodt = 0.0;
    double dVdx  = 0.0;
    double dVdy  = 0.0;

    double uu    = 0.0;
    double uv    = 0.0;
    double vv    = 0.0;

    for( int j=0; j<nnd; j++ )
    {
      NODE* node = elem->nd[j];

      double ndU = node->v.U;
      double ndV = node->v.V;

      U     +=    n[j] * ndU;
      Uo    +=    n[j] * node->vo.U;
      dUodt +=    n[j] * node->vo.dUdt;
      dUdx  += dndx[j] * ndU;
      dUdy  += dndy[j] * ndU;

      V     +=    n[j] * ndV;
      Vo    +=    n[j] * node->vo.V;
      dVodt +=    n[j] * node->vo.dVdt;
      dVdx  += dndx[j] * ndV;
      dVdy  += dndy[j] * ndV;

      uu    +=    n[j] * node->uu;
      uv    +=    n[j] * node->uv;
      vv    +=    n[j] * node->vv;
    }


    // compute dispersion coefficients ---------------------------------------------------

    double Dxx = 0.0;
    double Dxy = 0.0;
    double Dyy = 0.0;

    if( project->actualDisp > 0 )
    {
      for( int j=0; j<nnd; j++ )
      {
        NODE* node = elem->nd[j];
        Dxx += n[j] * node->Dxx;
        Dxy += n[j] * node->Dxy;
        Dyy += n[j] * node->Dyy;
      }
    }


    // compute eddy viscosity according to ELDER's assumption ----------------------------

    double vtxx = 0.0;
    double vtxy = 0.0;
    double vtyy = 0.0;

    for( int j=0; j<nnd; j++ )
    {
      NODE* node = elem->nd[j];

      vtxx += n[j] * node->exx * node->vt;
      vtxy += n[j] * node->exy * node->vt;
      vtyy += n[j] * node->eyy * node->vt;

// -- to force compatibility with former versions < 4.01, use following implementation
//      double vt = node->vt;
//      if( vt < project->minVt )  vt = project->minVt;
//      vtxx += n[j] * node->exx * vt;
//      vtxy += n[j] * node->exy * vt;
//      vtyy += n[j] * node->eyy * vt;
    }

// -- following implementation has been valid in versions from 4.01.00 to 4.01.06
//    if( isFS(project->actualTurb, BCONSET::kVtMin) )
//    {
//      if( vtxx < project->minVt )  vtxx = project->minVt;
//      if( vtyy < project->minVt )  vtyy = project->minVt;
//    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vtxx < project->minVtxx )  vtxx = project->minVtxx;
      if( vtyy < project->minVtyy )  vtyy = project->minVtyy;
    }

    // add eddy viscosity to kinematic viscosity -----------------------------------------

    vtxx += project->vk;
    vtxy += project->vk;
    vtyy += project->vk;


    // -----------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double  Ures = sqrt( U*U + V*V );          // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // ---------------------------------------------------------------------------------
      // compute x-momentum equation

      f   = H * (U - Uo) / dt;                         // time
      f  -= Ho * (1.0 - th) * dUodt;

      f  += tha * H * (U * dUdx  +  V * dUdy);         // convection

      fx  = thd * H * (vtxx * dUdx + vtxy * dUdy);     // eddy viscosity (approximation)
      fy  = thd * H * (vtxy * dUdx + vtyy * dUdy);

      fx -= thd * H * uu;                              // turbulence
      fy -= thd * H * uv;

      f  += thp * H * gravity * dadx;                  // gravity
      fx -= thp * H * H * gravity / 2.0;

      f  += thf * cf * Ures * U;                       // bottom friction

      // ------------------------------------------------ dispersion
      // fx += thd * H * Duu;
      // fy += thd * H * Duv;

      fx += thd * H * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );
      fy += thd * H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute y-momentum equation

      f   = H * (V - Vo) / dt;                         // time
      f  -= Ho * (1.0 - th) * dVodt;

      f  += tha * H * (U * dVdx  +  V * dVdy);         // convection

      fx  = thd * H * (vtxx * dVdx + vtxy * dVdy);     // eddy viscosity (approximation)
      fy  = thd * H * (vtxy * dVdx + vtyy * dVdy);

      fx -= thd * H * uv;                              // turbulence
      fy -= thd * H * vv;

      f  += thp * H * gravity * dady;                  // gravity
      fy -= thp * H * H * gravity / 2.0;

      f  += thf * cf * Ures * V;                       // bottom friction

      // ------------------------------------------------ dispersion
      // fx += thd * H * Duv;
      // fy += thd * H * Dvv;

      fx += thd * H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );
      fy += thd * H * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + startV;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute continuity equation (solve only for corner nodes)

      f  = (S - So) / dt;
      f -= (1.0 - th) * dHodt;
      f += thc * (H * (dUdx + dVdy) + U * dHdx + V * dHdy);
      f += SS;                                         // Source or Sink
      f *= weight;

      forcePtr = force + startS;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f;
      }
    }


    if( estifm )
    {
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfx_, dfxx, dfxy, dfy_, dfyx, dfyy;
      double* estifmPtr;

      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      double iUres;

      if( Ures > 1.0e-8 ) iUres = 1.0 / Ures;
      else                iUres = 0.0;


      // U-derivative of x-momentum ------------------------------------------------------

      df__  =  weight * H / dt;
      df__ +=  weight * tha * H * dUdx;
      df_x  =  weight * tha * H * U;
      df_y  =  weight * tha * H * V;

      dfxx  =  weight * thd * H * vtxx;
      dfxy  =  weight * thd * H * vtxy;

      dfyx  =  weight * thd * H * vtxy;
      dfyy  =  weight * thd * H * vtyy;

      df__ +=  weight * thf * cf * (iUres * U*U  +  Ures);

      dfx_  =  weight * thd * H * ( 2.0*U*Dxx - 2.0*V*Dxy );
      dfy_  =  weight * thd * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j]  +  dfxx * dndx[j]  +  dfxy * dndy[j];
        ty[j] = dfy_ * n[j]  +  dfyx * dndx[j]  +  dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // V-derivative of x-momentum ------------------------------------------------------

      df__  =  weight * tha * H * dUdy;

      df__ +=  weight * thf * cf * iUres * U * V;

      dfx_  =  weight * thd * H * ( 2.0*V*Dyy - 2.0*V*Dxy );
      dfy_  =  weight * thd * H * ( V*(Dxx-Dyy) - 2.0*V*Dxy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j] + startV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // H-derivative of x-momentum ------------------------------------------------------

      df__  =  weight * (U - Uo) / dt;
      df__ +=  weight * tha * (U * dUdx + V * dUdy);

      dfx_  =  weight * thd * (vtxx * dUdx  +  vtxy * dUdy);
      dfy_  =  weight * thd * (vtxy * dUdx  +  vtyy * dUdy);

      dfx_ -=  weight * thd * uu;
      dfy_ -=  weight * thd * uv;

      df__ +=  weight * thp * gravity * dadx;
      dfx_ -=  weight * thp * gravity * H;

      dfx_ +=  weight * thd * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );
      dfy_ +=  weight * thd * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j] + startS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // U-derivative of y-momentum ------------------------------------------------------

      df__  =  weight * tha * H * dVdx;

      df__ +=  weight * thf * cf * iUres * U * V;

      dfx_  =  weight * thd * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );
      dfy_  =  weight * thd * H * ( 2.0*V*Dxy + 2.0*U*Dyy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + startV];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // V-derivative of y-momentum ------------------------------------------------------

      df__  =  weight * H / dt;
      df__ +=  weight * tha * H * dVdy;
      df_x  =  weight * tha * H * U;
      df_y  =  weight * tha * H * V;

      dfxx  =  weight * thd * H * vtxx;
      dfxy  =  weight * thd * H * vtxy;

      dfyx  =  weight * thd * H * vtxy;
      dfyy  =  weight * thd * H * vtyy;

      df__ +=  weight * thf * cf * (iUres * V*V  +  Ures);

      dfx_  =  weight * thd * H * ( U*(Dxx-Dyy) - 2.0*V*Dxy );
      dfy_  =  weight * thd * H * ( 2.0*V*Dxx + 2.0*U*Dxy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j]  +  dfxx * dndx[j]  +  dfxy * dndy[j];
        ty[j] = dfy_ * n[j]  +  dfyx * dndx[j]  +  dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + startV] + startV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // H-derivative of y-momentum ------------------------------------------------------

      df__  =  weight * (V - Vo) / dt;
      df__ +=  weight * tha * (U * dVdx + V * dVdy);

      dfx_  =  weight * thd * (vtxx * dVdx + vtxy * dVdy);
      dfy_  =  weight * thd * (vtxy * dVdx + vtyy * dVdy);

      dfx_ -=  weight * thd * uv;
      dfy_ -=  weight * thd * vv;

      df__ +=  weight * thp * gravity * dady;
      dfy_ -=  weight * thp * gravity * H;

      dfx_ +=  weight * thd * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );
      dfy_ +=  weight * thd * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + startV] + startS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // U-derivative of continuity ------------------------------------------------------

      df__ = weight * thc * dHdx;
      df_x = weight * thc * H;

      for( int j=0; j<nnd; j++ )
      {
        t[j] = df__ * n[j]  +  df_x * dndx[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startS];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }


      // V-derivative of continuity ------------------------------------------------------

      df__ = weight * thc * dHdy;
      df_y = weight * thc * H;

      for( int j=0; j<nnd; j++ )
      {
        t[j] = df__ * n[j] + df_y * dndy[j];
      }


      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startS] + startV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }


      // H-derivative of continuity ------------------------------------------------------

      df__  = weight / dt;
      df__ += weight * thc * (dUdx + dVdy);
      df_x  = weight * thc * U;
      df_y  = weight * thc * V;


      for( int j=0; j<ncn; j++ )
      {
        t[j] = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
      }


      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startS] + startS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }
    }
  } // END of loop over all GAUSS point



  // apply transformation ----------------------------------------------------------------

  Rotate2D( nnd, elem->nd, 3, estifm, force );


  // -------------------------------------------------------------------------------------
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

    if(isFS(bcon->kind,BCON::kInlet) )
    {
      // set equation row dfU to zero ----------------------------------------------------

      if( estifm )
      {
        for( int j=0; j<maxEleq; j++ )
        {
          estifm[i][j] = 0.0;
        }
      }


      // compute flow in normal direction at node i --------------------------------------

      double Un    = node->v.U * bcon->niox  +  node->v.V * bcon->nioy;
      double H     = node->v.S - node->z;
      double specQ = bcon->val->U;

      if( H <= 0.0 )
      {
        H = project->hmin;
        specQ = 0.0;
      }


      // force vector --------------------------------------------------------------------

      if( force )  force[i] = area * (specQ - Un*H);


      // jacobi matrix -------------------------------------------------------------------

      if( estifm )
      {
        estifm[i][i]          = area * H;
        estifm[i][i + startS] = area * Un;
      }
    }
  }


  for( int i=ncn; i<nnd; i++ )    // midside nodes
  {
    NODE* node = elem->nd[i];
    BCON* bcon = &node->bc;

    if( isFS(bcon->kind, BCON::kInlet) )
    {
      // get left and right corner node to midside node i --------------------------------

      int l, r;

      qShape->getCornerNodes( i, &l, &r );

      NODE* lnode = elem->nd[l];
      NODE* rnode = elem->nd[r];

      int nofeqHl = startS + l;    // dfH at left corner node
      int nofeqHr = startS + r;    // dfH at right corner node


      // set equation row dfU to zero ----------------------------------------------------

      if( estifm )
      {
        for( int j=0; j<maxEleq; j++ )
        {
          estifm[i][j] = 0.0;
        }
      }


      // compute flow in normal direction at node i --------------------------------------

      double Un    = node->v.U * bcon->niox  +  node->v.V * bcon->nioy;
      double Hl    = lnode->v.S - lnode->z;
      double Hr    = rnode->v.S - rnode->z;
      double H     = ( Hl + Hr ) / 2.0;
      double specQ = bcon->val->U;

      if( H <= 0.0 )
      {
        H     = project->hmin;
        specQ = 0.0;
      }


      // force vector --------------------------------------------------------------------

      if( force )  force[i] = area * (specQ - Un*H);


      // jacobi matrix -------------------------------------------------------------------

      if( estifm )
      {
        estifm[i][i]       = area * H;
        estifm[i][nofeqHl] = area * Un / 2.0;
        estifm[i][nofeqHr] = area * Un / 2.0;
      }
    }
  }
}
*/

// ======================================================================================

void EQS_UVS2D_AI::Region( ELEM*    elem,
                            PROJECT* project,
                            double** estifm,
                            double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  int ngp = qShape->ngp;         // number of GAUSS points
  int nnd = qShape->nnd;         // number of nodes in all
  int ncn = lShape->nnd;         // number of corner nodes

  TYPE* type = TYPE::Getid( elem->type );

  int startV = nnd;
  int startS = 2 * nnd;

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

  double gravity = project->g;
  double area    = 0.0;


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

  for( int g=0; g<ngp; g++ ) // START of loop over all GAUSS point
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with quadratic shape functions

    double* dfdxPtr = qShape->dfdx[g];
    double* dfdyPtr = qShape->dfdy[g];

    double trafo[2][2];

    double detj = qShape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );

    double weight = detj * qShape->weight[g];

    area += weight;


    // -----------------------------------------------------------------------------------
    // compute values of quadratic shape functions at GP g

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = qShape->f[g];

    for( int j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // ----------------------------------------------------------------------------------
    // compute values of linear shape functions at GP g

    dfdxPtr = lShape->dfdx[g];
    dfdyPtr = lShape->dfdy[g];

    double dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    double* m = lShape->f[g];

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
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

    for( int j=0; j<ncn; j++ )
    {
      NODE* node = elem->nd[j];
      BCON* bcon = &node->bc;

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
      // -----------------------------------------------------------------------------------

      H    +=    m[j] * ndH;
      dHdx += dmdx[j] * ndH;
      dHdy += dmdy[j] * ndH;

      dadx += dmdx[j] * ndZ;
      dady += dmdy[j] * ndZ;

      dHdt +=    m[j] * node->v.dSdt;

      cf   +=    m[j] * node->cf;

      SS   +=    m[j] * ndSS;            // Source or Sink
    }

//  if( H <= 0.0 )  H = project->hmin;


    // integrate U, V, uu, uv and vv with quadratic shape

    double U    = 0.0;
    double dUdt = 0.0;
    double dUdx = 0.0;
    double dUdy = 0.0;

    double V    = 0.0;
    double dVdt = 0.0;
    double dVdx = 0.0;
    double dVdy = 0.0;

    double uu   = 0.0;
    double uv   = 0.0;
    double vv   = 0.0;

    for( int j=0; j<nnd; j++ )
    {
      NODE* node = elem->nd[j];

      double ndU = node->v.U;
      double ndV = node->v.V;

      U     +=    n[j] * ndU;
      dUdt  +=    n[j] * node->v.dUdt;
      dUdx  += dndx[j] * ndU;
      dUdy  += dndy[j] * ndU;

      V     +=    n[j] * ndV;
      dVdt  +=    n[j] * node->v.dVdt;
      dVdx  += dndx[j] * ndV;
      dVdy  += dndy[j] * ndV;

      uu    +=    n[j] * node->uu;
      uv    +=    n[j] * node->uv;
      vv    +=    n[j] * node->vv;
    }


    // compute dispersion coefficients ---------------------------------------------------

    double Dxx = 0.0;
    double Dxy = 0.0;
    double Dyy = 0.0;

    if( project->actualDisp > 0 )
    {
      for( int j=0; j<nnd; j++ )
      {
        NODE* node = elem->nd[j];
        Dxx += n[j] * node->Dxx;
        Dxy += n[j] * node->Dxy;
        Dyy += n[j] * node->Dyy;
      }
    }


    // compute eddy viscosity according to ELDER's assumption ----------------------------

    double vtxx = 0.0;
    double vtxy = 0.0;
    double vtyy = 0.0;

    for( int j=0; j<nnd; j++ )
    {
      NODE* node = elem->nd[j];

      vtxx += n[j] * node->exx * node->vt;
      vtxy += n[j] * node->exy * node->vt;
      vtyy += n[j] * node->eyy * node->vt;
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


    // -----------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double  Ures = sqrt( U*U + V*V );          // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // ---------------------------------------------------------------------------------
      // compute x-momentum equation

      f   = H * dUdt;                                  // time

      f  += H * (U * dUdx  +  V * dUdy);               // convection

      fx  = H * (vtxx * dUdx + vtxy * dUdy);           // eddy viscosity (approximation)
      fy  = H * (vtxy * dUdx + vtyy * dUdy);

      fx -= H * uu;                                    // turbulence
      fy -= H * uv;

      f  += H * gravity * dadx;                        // gravity
      fx -= H * H * gravity / 2.0;

      f  += cf * Ures * U;                             // bottom friction

      // ------------------------------------------------ dispersion
      // fx += H * Duu;
      // fy += H * Duv;

      fx += H * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );
      fy += H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute y-momentum equation

      f   = H * dVdt;                                  // time

      f  += H * (U * dVdx  +  V * dVdy);               // convection

      fx  = H * (vtxx * dVdx + vtxy * dVdy);           // eddy viscosity (approximation)
      fy  = H * (vtxy * dVdx + vtyy * dVdy);

      fx -= H * uv;                                    // turbulence
      fy -= H * vv;

      f  += H * gravity * dady;                        // gravity
      fy -= H * H * gravity / 2.0;

      f  += cf * Ures * V;                             // bottom friction

      // ------------------------------------------------ dispersion
      // fx += H * Duv;
      // fy += H * Dvv;

      fx += H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );
      fy += H * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + startV;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute continuity equation (solve only for corner nodes)

      f  = dHdt  +  H * (dUdx + dVdy)  +  U * dHdx  +  V * dHdy;

      f *= weight;
      f += SS * weight;            // Source or Sink

      forcePtr = force + startS;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f;
      }
    }


    if( estifm )
    {
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfx_, dfxx, dfxy, dfy_, dfyx, dfyy;
      double* estifmPtr;

      // ---------------------------------------------------------------------------------
      // compute components of NEWTON-RAPHSON Jacobi matrix

      double iUres;

      if( Ures > 1.0e-8 ) iUres = 1.0 / Ures;
      else                iUres = 0.0;


      // U-derivative of x-momentum ------------------------------------------------------

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

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j]  +  dfxx * dndx[j]  +  dfxy * dndy[j];
        ty[j] = dfy_ * n[j]  +  dfyx * dndx[j]  +  dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // V-derivative of x-momentum ------------------------------------------------------

      df__  =  weight * H * dUdy;

      df__ +=  weight * cf * iUres * U * V;

      dfx_  =  weight * H * ( 2.0*V*Dyy - 2.0*V*Dxy );
      dfy_  =  weight * H * ( V*(Dxx-Dyy) - 2.0*V*Dxy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j] + startV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // H-derivative of x-momentum ------------------------------------------------------

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

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j] + startS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // U-derivative of y-momentum ------------------------------------------------------

      df__  =  weight * H * dVdx;

      df__ +=  weight * cf * iUres * U * V;

      dfx_  =  weight * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );
      dfy_  =  weight * H * ( 2.0*V*Dxy + 2.0*U*Dyy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + startV];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // V-derivative of y-momentum ------------------------------------------------------

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

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j]  +  dfxx * dndx[j]  +  dfxy * dndy[j];
        ty[j] = dfy_ * n[j]  +  dfyx * dndx[j]  +  dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + startV] + startV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // H-derivative of y-momentum ------------------------------------------------------

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

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        tx[j] = dfx_ * m[j];
        ty[j] = dfy_ * m[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + startV] + startS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // U-derivative of continuity ------------------------------------------------------

      df__ = weight * dHdx;
      df_x = weight * H;

      for( int j=0; j<nnd; j++ )
      {
        t[j] = df__ * n[j]  +  df_x * dndx[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startS];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }


      // V-derivative of continuity ------------------------------------------------------

      df__ = weight * dHdy;
      df_y = weight * H;

      for( int j=0; j<nnd; j++ )
      {
        t[j] = df__ * n[j] + df_y * dndy[j];
      }


      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startS] + startV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
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


      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startS] + startS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }
    }
  } // END of loop over all GAUSS point



  // apply transformation ----------------------------------------------------------------

  Rotate2D( nnd, elem->nd, 3, estifm, force );


  // -------------------------------------------------------------------------------------
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

    if(isFS(bcon->kind,BCON::kInlet) )
    {
      // set equation row dfU to zero ----------------------------------------------------

      if( estifm )
      {
        for( int j=0; j<maxEleq; j++ )
        {
          estifm[i][j] = 0.0;
        }
      }


      // compute flow in normal direction at node i --------------------------------------

      double Un    = node->v.U * bcon->niox  +  node->v.V * bcon->nioy;
      double H     = node->v.S - node->z;
      double specQ = bcon->val->U;

      if( H <= 0.0 )
      {
        H = project->hmin;
        specQ = 0.0;
      }


      // force vector --------------------------------------------------------------------

      if( force )  force[i] = area * (specQ - Un*H);


      // jacobi matrix -------------------------------------------------------------------

      if( estifm )
      {
        estifm[i][i]          = area * H;
        estifm[i][i + startS] = area * Un;
      }
    }
  }


  for( int i=ncn; i<nnd; i++ )    // midside nodes
  {
    NODE* node = elem->nd[i];
    BCON* bcon = &node->bc;

    if( isFS(bcon->kind, BCON::kInlet) )
    {
      // get left and right corner node to midside node i --------------------------------

      int l, r;

      qShape->getCornerNodes( i, &l, &r );

      NODE* lnode = elem->nd[l];
      NODE* rnode = elem->nd[r];

      int nofeqHl = startS + l;    // dfH at left corner node
      int nofeqHr = startS + r;    // dfH at right corner node


      // set equation row dfU to zero ----------------------------------------------------

      if( estifm )
      {
        for( int j=0; j<maxEleq; j++ )
        {
          estifm[i][j] = 0.0;
        }
      }


      // compute flow in normal direction at node i --------------------------------------

      double Un    = node->v.U * bcon->niox  +  node->v.V * bcon->nioy;
      double Hl    = lnode->v.S - lnode->z;
      double Hr    = rnode->v.S - rnode->z;
      double H     = ( Hl + Hr ) / 2.0;
      double specQ = bcon->val->U;

      if( H <= 0.0 )
      {
        H     = project->hmin;
        specQ = 0.0;
      }


      // force vector --------------------------------------------------------------------

      if( force )  force[i] = area * (specQ - Un*H);


      // jacobi matrix -------------------------------------------------------------------

      if( estifm )
      {
        estifm[i][i]       = area * H;
        estifm[i][nofeqHl] = area * Un / 2.0;
        estifm[i][nofeqHr] = area * Un / 2.0;
      }
    }
  }
}
