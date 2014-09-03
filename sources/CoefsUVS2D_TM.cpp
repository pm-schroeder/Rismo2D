// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_TM
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

#include "EqsUVS2D_TM.h"



// ======================================================================================

int EQS_UVS2D_TM::Coefs( ELEM*    elem,
                         PROJECT* project,
                         double** estifm,
                         double*  force )
{
  if( isFS(elem->flag, ELEM::kDry) ) return 0;

  if( this->timegrad )
  {
    if( isFS(elem->flag, ELEM::kBound) )
    {
      Bound_dt( elem, project, estifm, force );
    }
    else
    {
      Region_dt( elem, project, estifm, force );
    }
  }

  else
  {
    if( isFS(elem->flag, ELEM::kBound) )
    {
      Bound( elem, project, estifm, force );
    }
    else
    {
      Region( elem, project, estifm, force );
    }
  }

  return 1;
}



// ======================================================================================

void EQS_UVS2D_TM::Bound( ELEM*    elem,
                          PROJECT* project,
                          double** estifm,
                          double*  force )
{
  SHAPE* lShape  = elem->GetLShape();
  SHAPE* qShape  = elem->GetQShape();

  int    nnd     = qShape->nnd;

  int    eqidV   = nnd;
  int    eqidS   = 2 * nnd;

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


  // row index for force vector and element stiffness matrix -----------------------------

  int r0U = 0;       int r1U = 1;           int r2U = 2;
  int r0V = eqidV;   int r1V = eqidV + 1;   int r2V = eqidV + 2;


  // column index for element stiffness matrix -------------------------------------------

  int c0H = eqidS;   int c1H = eqidS + 1;


  for( int g=0; g<qShape->ngp; g++ )   // loop on GAUSS points
  {
    double* m  = lShape->f[g];         // linear shape
    double* n  = qShape->f[g];         // quadratic shape
    double* dn = qShape->dfdx[g];


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

    double nx =  dn[0]*node[0]->y + dn[1]*node[1]->y + dn[2]*node[2]->y;
    double ny = -dn[0]*node[0]->x - dn[1]*node[1]->x - dn[2]*node[2]->x;

    double len = sqrt( nx*nx + ny*ny );


    // weight of Gauss point j -----------------------------------------------------------

    double weight = qShape->weight[g];

    double U    = 0.0;
    double V    = 0.0;
    double Ures = 0.0;
    double cfw  = 0.0;

    if( solidFlag )
    {
      U = n[0]*node[0]->v.U + n[1]*node[1]->v.U + n[2]*node[2]->v.U;
      V = n[0]*node[0]->v.V + n[1]*node[1]->v.V + n[2]*node[2]->v.V;

      Ures = sqrt( U * U  +  V * V );

      cfw = m[0]*node[0]->cfw + m[1]*node[1]->cfw;
    }


    // -----------------------------------------------------------------------------------

    if( force )
    {
      double fU = weight * thp * H * ( gravity*H/2.0 * nx );
      double fV = weight * thp * H * ( gravity*H/2.0 * ny );


      // wall roughness ------------------------------------------------------------------

      if( solidFlag )
      {
        fU += weight * thf * len * H * cfw * U * Ures;
        fV += weight * thf * len * H * cfw * V * Ures;
      }


      // add fU, fV to force vector ------------------------------------------------------

      force[r0U] -= n[0] * fU;        force[r0V] -= n[0] * fV;
      force[r1U] -= n[1] * fU;        force[r1V] -= n[1] * fV;
      force[r2U] -= n[2] * fU;        force[r2V] -= n[2] * fV;
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

      estifm[r0U][r0U] += n[0] * dfUdU * n[0];
      estifm[r0U][r1U] += n[0] * dfUdU * n[1];
      estifm[r0U][r2U] += n[0] * dfUdU * n[2];

      estifm[r1U][r0U] += n[1] * dfUdU * n[0];
      estifm[r1U][r1U] += n[1] * dfUdU * n[1];
      estifm[r1U][r2U] += n[1] * dfUdU * n[2];

      estifm[r2U][r0U] += n[2] * dfUdU * n[0];
      estifm[r2U][r1U] += n[2] * dfUdU * n[1];
      estifm[r2U][r2U] += n[2] * dfUdU * n[2];

      estifm[r0U][r0V] += n[0] * dfUdV * n[0];
      estifm[r0U][r1V] += n[0] * dfUdV * n[1];
      estifm[r0U][r2V] += n[0] * dfUdV * n[2];

      estifm[r1U][r0V] += n[1] * dfUdV * n[0];
      estifm[r1U][r1V] += n[1] * dfUdV * n[1];
      estifm[r1U][r2V] += n[1] * dfUdV * n[2];

      estifm[r2U][r0V] += n[2] * dfUdV * n[0];
      estifm[r2U][r1V] += n[2] * dfUdV * n[1];
      estifm[r2U][r2V] += n[2] * dfUdV * n[2];

      estifm[r0V][r0U] += n[0] * dfVdU * n[0];
      estifm[r0V][r1U] += n[0] * dfVdU * n[1];
      estifm[r0V][r2U] += n[0] * dfVdU * n[2];

      estifm[r1V][r0U] += n[1] * dfVdU * n[0];
      estifm[r1V][r1U] += n[1] * dfVdU * n[1];
      estifm[r1V][r2U] += n[1] * dfVdU * n[2];

      estifm[r2V][r0U] += n[2] * dfVdU * n[0];
      estifm[r2V][r1U] += n[2] * dfVdU * n[1];
      estifm[r2V][r2U] += n[2] * dfVdU * n[2];

      estifm[r0V][r0V] += n[0] * dfVdV * n[0];
      estifm[r0V][r1V] += n[0] * dfVdV * n[1];
      estifm[r0V][r2V] += n[0] * dfVdV * n[2];

      estifm[r1V][r0V] += n[1] * dfVdV * n[0];
      estifm[r1V][r1V] += n[1] * dfVdV * n[1];
      estifm[r1V][r2V] += n[1] * dfVdV * n[2];

      estifm[r2V][r0V] += n[2] * dfVdV * n[0];
      estifm[r2V][r1V] += n[2] * dfVdV * n[1];
      estifm[r2V][r2V] += n[2] * dfVdV * n[2];

      estifm[r0U][c0H] += n[0] * dfUdH * m[0];
      estifm[r0V][c0H] += n[0] * dfVdH * m[0];

      estifm[r1U][c0H] += n[1] * dfUdH * m[0];
      estifm[r1V][c0H] += n[1] * dfVdH * m[0];

      estifm[r2U][c0H] += n[2] * dfUdH * m[0];
      estifm[r2V][c0H] += n[2] * dfVdH * m[0];

      estifm[r0U][c1H] += n[0] * dfUdH * m[1];
      estifm[r0V][c1H] += n[0] * dfVdH * m[1];

      estifm[r1U][c1H] += n[1] * dfUdH * m[1];
      estifm[r1V][c1H] += n[1] * dfVdH * m[1];

      estifm[r2U][c1H] += n[2] * dfUdH * m[1];
      estifm[r2V][c1H] += n[2] * dfVdH * m[1];
    }
  }


  // apply transformation ----------------------------------------------------------------

  Rotate2D( nnd, eqidV, elem->nd, estifm, force );


  // -------------------------------------------------------------------------------------
  // In case of experimental boundary forces on inlet set equation row to zero.
  // Accounting for boundary forces is done in method EQS_UVS2D_TM::Region().

  for( int i=0; i<nnd; i++ )
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


// ======================================================================================

void EQS_UVS2D_TM::Region( ELEM*    elem,
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
    NODE *nd = elem->nd[i];

    x[i] = nd->x - x[0];
    y[i] = nd->y - y[0];
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


    // -------------------------------------------------------------------------------------
    // compute values of quadratic shape functions at GP g

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = qShape->f[g];

    for( int j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // ------------------------------------------------------------------------------------
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

    for( int j=0; j<ncn; j++ )
    {
      NODE* node = elem->nd[j];
      BCON* bcon = &node->bc;

      double ndS = node->v.S;
      double ndZ = node->z;
      double ndH = ndS - ndZ;

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

      S    +=    m[j] * ndS;
      So   +=    m[j] * node->vo.S;

      H    +=    m[j] * ndH;
      dSodt +=   m[j] * node->vo.dSdt;
      dHdx += dmdx[j] * ndH;
      dHdy += dmdy[j] * ndH;

      dadx += dmdx[j] * ndZ;
      dady += dmdy[j] * ndZ;

      cf   +=    m[j] * node->cf;

      SS   +=    m[j] * ndSS;            // Source or Sink
    }


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


    // compute eddy viscosity ------------------------------------------------------------

    double vt = 0.0;

    if( isFS(project->actualTurb, BCONSET::kVtConstant) )
    {
      vt = type->vt;
    }
    else
    {
      for( int j=0; j<nnd; j++ )  vt += n[j] * elem->nd[j]->vt;
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < type->vt )  vt = type->vt;
    }


    // add kinematic viscosity to eddy viscosity  ----------------------------------------

    vt += project->vk;


    // -----------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double  Ures = sqrt( U*U + V*V );          // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // ---------------------------------------------------------------------------------
      // compute x-momentum equation

      f   =  H * (U - Uo) / dt;                             // time
      f  -=  H * (1.0 - th) * dUodt;

      f  +=  tha * H * (U * dUdx  +  V * dUdy);             // convection
      fx  =  0.0;
      fy  =  0.0;

      fx +=  thd * H * vt * (dUdx + dUdx);                  // eddy viscosity
      fy +=  thd * H * vt * (dUdy + dVdx);

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // fx +=  thd * vt * (U * dHdx + U * dHdx);
      // fy +=  thd * vt * (U * dHdy + V * dHdx);
      // -------------------------------------------------------------------------------------------

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

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute y-momentum equation

      f   =  H * (V - Vo) / dt;                             // time
      f  -=  H * (1.0 - th) * dVodt;

      f  +=  tha * H * (U * dVdx  +  V * dVdy);             // convection
      fx  =  0.0;
      fy  =  0.0;

      fx +=  thd * H * vt * (dVdx + dUdy);                  // eddy viscosity
      fy +=  thd * H * vt * (dVdy + dVdy);

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // fx +=  thd * vt * (V * dHdx + U * dHdy);
      // fy +=  thd * vt * (V * dHdy + V * dHdy);
      // -------------------------------------------------------------------------------------------

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

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute continuity equation (solve only for corner nodes)

      f  = (S - So) / dt;
      f -= (1.0 - th) * dSodt;
      f += thc * (H * (dUdx + dVdy) + U * dHdx + V * dHdy);
      f += SS;                                              // Source or Sink
      f *= weight;

      forcePtr = force + eqidS;

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

      if( Ures > 1.0e-9 ) iUres = 1.0 / Ures;
      else                iUres = 0.0;


      // U-derivative of x-momentum ------------------------------------------------------

      df__  =  weight * H / dt;
      df__ +=  weight * tha * H * dUdx;
      df_x  =  weight * tha * H * U;
      df_y  =  weight * tha * H * V;
      dfx_  =  0.0;
      dfy_  =  0.0;

      dfxx  =  weight * thd * H * vt * 2.0;
      dfyy  =  weight * thd * H * vt;

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // dfx_ +=  weight * thd * vt * dHdx * 2.0;
      // dfy_ +=  weight * thd * vt * dHdy;
      // -------------------------------------------------------------------------------------------

      df__ +=  weight * thf * cf * (iUres * U*U  +  Ures);

      dfx_ +=  weight * thd * H * ( 2.0*U*Dxx - 2.0*V*Dxy );
      dfy_ +=  weight * thd * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j]  +  dfxx * dndx[j];
        ty[j] = dfy_ * n[j]                     +  dfyy * dndy[j];
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
      dfx_  =  0.0;
      dfy_  =  0.0;

      dfyx  =  weight * thd * H * vt;

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // dfy_ +=  weight * thd * vt * dHdx;
      // -------------------------------------------------------------------------------------------

      df__ +=  weight * thf * cf * iUres * U * V;

      dfx_ +=  weight * thd * H * ( 2.0*V*Dyy - 2.0*V*Dxy );
      dfy_ +=  weight * thd * H * ( V*(Dxx-Dyy) - 2.0*V*Dxy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j];
        tx[j] = dfx_ * n[j];
        ty[j] = dfy_ * n[j]  +  dfyx * dndx[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j] + eqidV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // H-derivative of x-momentum ------------------------------------------------------

      df__  =  weight * (U - Uo) / dt;

      df__ +=  weight * tha * (U * dUdx + V * dUdy);
      dfx_  =  0.0;
      dfy_  =  0.0;

      dfx_ +=  weight * thd * vt * (dUdx + dUdx);
      dfy_ +=  weight * thd * vt * (dUdy + dVdx);

      dfx_ -=  weight * thd * uu;
      dfy_ -=  weight * thd * uv;

      df__ +=  weight * thp * gravity * dadx;
      dfx_ -=  weight * thp * gravity * H;

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // dfxx  =  weight * thd * vt * U * 2.0;
      // dfyx  =  weight * thd * vt * V;
      // dfyy  =  weight * thd * vt * U;
      // -------------------------------------------------------------------------------------------

      dfxx  =  0.0;
      dfyx  =  0.0;
      dfyy  =  0.0;

      dfx_ +=  weight * thd * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );
      dfy_ +=  weight * thd * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        tx[j] = dfx_ * m[j]  +  dfxx * dmdx[j];
        ty[j] = dfy_ * m[j]  +  dfyx * dmdx[j]  +  dfyy * dmdy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j] + eqidS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // U-derivative of y-momentum ------------------------------------------------------

      df__  =  weight * tha * H * dVdx;
      dfx_  =  0.0;
      dfy_  =  0.0;

      dfxy  =  weight * thd * H * vt;

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // dfx_ +=  weight * thd * vt * dHdy;
      // -------------------------------------------------------------------------------------------

      df__ +=  weight * thf * cf * iUres * U * V;

      dfx_ +=  weight * thd * H * ( V*(Dxx-Dyy) + 2.0*U*Dxy );
      dfy_ +=  weight * thd * H * ( 2.0*V*Dxy + 2.0*U*Dyy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j];
        tx[j] = dfx_ * n[j]  +  dfxy * dndy[j];
        ty[j] = dfy_ * n[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + eqidV];

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
      dfx_  =  0.0;
      dfy_  =  0.0;

      dfxx  =  weight * thd * H * vt;
      dfyy  =  weight * thd * H * vt * 2.0;

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // dfx_ +=  weight * thd * vt * dHdx;
      // dfy_ +=  weight * thd * vt * dHdy * 2.0;
      // -------------------------------------------------------------------------------------------

      df__ +=  weight * thf * cf * (iUres * V*V  +  Ures);

      dfx_ +=  weight * thd * H * ( U*(Dxx-Dyy) - 2.0*V*Dxy );
      dfy_ +=  weight * thd * H * ( 2.0*V*Dxx + 2.0*U*Dxy );

      for( int j=0; j<nnd; j++ )
      {
        t[j]  = df__ * n[j]  +  df_x * dndx[j]  +  df_y * dndy[j];
        tx[j] = dfx_ * n[j]  +  dfxx * dndx[j];
        ty[j] = dfy_ * n[j]  +  dfyy * dndy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + eqidV] + eqidV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j]*t[k] + dndx[j]*tx[k] + dndy[j]*ty[k];
        }
      }


      // H-derivative of y-momentum ------------------------------------------------------

      df__  =  weight * (V - Vo) / dt;

      df__ +=  weight * tha * (U * dVdx + V * dVdy);
      dfx_  =  0.0;
      dfy_  =  0.0;

      dfx_ +=  weight * thd * vt * (dVdx + dUdy);
      dfy_ +=  weight * thd * vt * (dVdy + dVdy);

      dfx_ -=  weight * thd * uv;
      dfy_ -=  weight * thd * vv;

      df__ +=  weight * thp * gravity * dady;
      dfy_ -=  weight * thp * gravity * H;

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // dfxx  =  weight * thd * vt * V;
      // dfxy  =  weight * thd * vt * U;
      // dfyy  =  weight * thd * vt * V * 2.0;
      // -------------------------------------------------------------------------------------------

      dfxx  =  0.0;
      dfxy  =  0.0;
      dfyy  =  0.0;

      dfx_ +=  weight * thd * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );
      dfy_ +=  weight * thd * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        tx[j] = dfx_ * m[j]  +  dfxx * dmdx[j]  +  dfxy * dmdy[j];
        ty[j] = dfy_ * m[j]                     +  dfyy * dmdy[j];
      }

      for( int j=0; j<nnd; j++ )
      {
        estifmPtr = estifm[j + eqidV] + eqidS;

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
        estifmPtr = estifm[j + eqidS];

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
        estifmPtr = estifm[j + eqidS] + eqidV;

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
        estifmPtr = estifm[j + eqidS] + eqidS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }
    }
  } // END of loop over all GAUSS points



  // apply transformation --------------------------------------------------------------------------

  Rotate2D( nnd, eqidV, elem->nd, estifm, force );


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

  for( int i=ncn; i<nnd; i++ )    // midside nodes
  {
    NODE* node = elem->nd[i];
    BCON* bcon = &node->bc;

    if( isFS(bcon->kind, BCON::kInlet) )
    {
      // get left and right corner node to midside node i ------------------------------------------

      int l, r;

      qShape->getCornerNodes( i, &l, &r );

      NODE* lnode = elem->nd[l];
      NODE* rnode = elem->nd[r];

      int nofeqHl = eqidS + l;    // dfH at left corner node
      int nofeqHr = eqidS + r;    // dfH at right corner node


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
      double Hl    = lnode->v.S - lnode->z;
      double Hr    = rnode->v.S - rnode->z;
      double H     = ( Hl + Hr ) / 2.0;
      double specQ = bcon->val->U;

      if( H <= 0.0 )
      {
        H     = project->hmin;
        specQ = 0.0;
      }


      // force vector ------------------------------------------------------------------------------

      if( force )  force[i] = area * (specQ - Un*H);


      // jacobi matrix -----------------------------------------------------------------------------

      if( estifm )
      {
        estifm[i][i]       = area * H;
        estifm[i][nofeqHl] = area * Un / 2.0;
        estifm[i][nofeqHr] = area * Un / 2.0;
      }
    }
  }
}


// ======================================================================================

void EQS_UVS2D_TM::Bound_dt( ELEM*    elem,
                             PROJECT* project,
                             double** estifm,
                             double*  force )
{
  SHAPE* lShape  = elem->GetLShape();
  SHAPE* qShape  = elem->GetQShape();

  int    nnd     = qShape->nnd;

  int    eqidV   = nnd;

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
  node[2] = elem->nd[2];      // midside node


  // row index for force vector ----------------------------------------------------------

  int r0U = 0;       int r1U = 1;           int r2U = 2;
  int r0V = eqidV;   int r1V = eqidV + 1;   int r2V = eqidV + 2;

  for( int g=0; g<qShape->ngp; g++ )   // loop on GAUSS points
  {
    double* m  = lShape->f[g];         // linear shape
    double* n  = qShape->f[g];         // quadratic shape
    double* dn = qShape->dfdx[g];


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


    // -----------------------------------------------------------------------------------
    // compute normal vector
    // since the normal is not reduced to unit length it
    // implies the transformation of the integrand

    double nx =  dn[0]*node[0]->y + dn[1]*node[1]->y + dn[2]*node[2]->y;
    double ny = -dn[0]*node[0]->x - dn[1]*node[1]->x - dn[2]*node[2]->x;

    double len = sqrt( nx*nx + ny*ny );


    // weight of Gauss point j -----------------------------------------------------------

    double weight = qShape->weight[g];

    double U    = 0.0;
    double V    = 0.0;
    double Ures = 0.0;
    double cfw  = 0.0;

    if( solidFlag )
    {
      U = n[0]*node[0]->v.U + n[1]*node[1]->v.U + n[2]*node[2]->v.U;
      V = n[0]*node[0]->v.V + n[1]*node[1]->v.V + n[2]*node[2]->v.V;

      Ures = sqrt( U * U  +  V * V );

      cfw = m[0]*node[0]->cfw + m[1]*node[1]->cfw;
    }


    // -----------------------------------------------------------------------------------

    double fU = weight * H * ( gravity*H/2.0 * nx );
    double fV = weight * H * ( gravity*H/2.0 * ny );


    // wall roughness ------------------------------------------------------------------

    if( solidFlag )
    {
      fU += weight * len * H * cfw * U * Ures;
      fV += weight * len * H * cfw * V * Ures;
    }


    // add fU, fV to force vector ------------------------------------------------------

    force[r0U] -= n[0] * fU;        force[r0V] -= n[0] * fV;
    force[r1U] -= n[1] * fU;        force[r1V] -= n[1] * fV;
    force[r2U] -= n[2] * fU;        force[r2V] -= n[2] * fV;
  }


  // apply transformation ----------------------------------------------------------------

  Rotate2D( nnd, eqidV, elem->nd, estifm, force );


  // -------------------------------------------------------------------------------------
  // In case of experimental boundary forces on inlet set equation row to zero.
  // Accounting for boundary forces is done in method EQS_UVS2D_TM::Region().

  for( int i=0; i<nnd; i++ )
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

// ======================================================================================

void EQS_UVS2D_TM::Region_dt( ELEM*    elem,
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

  double gravity = project->g;


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    NODE *nd = elem->nd[i];

    x[i] = nd->x - x[0];
    y[i] = nd->y - y[0];
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


    // -------------------------------------------------------------------------------------
    // compute values of quadratic shape functions at GP g

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = qShape->f[g];

    for( int j=0; j<nnd; j++ )
    {
      dndx[j] = trafo[0][0] * dfdxPtr[j]  +  trafo[0][1] * dfdyPtr[j];
      dndy[j] = trafo[1][0] * dfdxPtr[j]  +  trafo[1][1] * dfdyPtr[j];
    }


    // ------------------------------------------------------------------------------------
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

    // integrate H, a and cf with linear shape

    double S    = 0.0;

    double H    = 0.0;
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

      double ndS = node->v.S;
      double ndZ = node->z;
      double ndH = ndS - ndZ;

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

      S    +=    m[j] * ndS;

      H    +=    m[j] * ndH;
      dHdx += dmdx[j] * ndH;
      dHdy += dmdy[j] * ndH;

      dadx += dmdx[j] * ndZ;
      dady += dmdy[j] * ndZ;

      cf   +=    m[j] * node->cf;

      SS   +=    m[j] * ndSS;            // Source or Sink
    }


    // integrate U, V, uu, uv and vv with quadratic shape

    double U    = 0.0;
    double dUdx = 0.0;
    double dUdy = 0.0;

    double V    = 0.0;
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
      dUdx  += dndx[j] * ndU;
      dUdy  += dndy[j] * ndU;

      V     +=    n[j] * ndV;
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


    // compute eddy viscosity ------------------------------------------------------------

    double vt = 0.0;

    if( isFS(project->actualTurb, BCONSET::kVtConstant) )
    {
      vt = type->vt;
    }
    else
    {
      for( int j=0; j<nnd; j++ )  vt += n[j] * elem->nd[j]->vt;
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < type->vt )  vt = type->vt;
    }


    // add kinematic viscosity to eddy viscosity  ----------------------------------------

    vt += project->vk;


    // -----------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double  Ures = sqrt( U*U + V*V );          // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;

      // ---------------------------------------------------------------------------------
      // compute x-momentum equation

      f   =  H * (U * dUdx  +  V * dUdy);                   // convection
      fx  =  0.0;
      fy  =  0.0;

      fx +=  H * vt * (dUdx + dUdx);                        // eddy viscosity
      fy +=  H * vt * (dUdy + dVdx);

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // fx +=  vt * (U * dHdx + U * dHdx);
      // fy +=  vt * (U * dHdy + V * dHdx);
      // -------------------------------------------------------------------------------------------

      fx -=  H * uu;                                        // turbulence
      fy -=  H * uv;

      f  +=  H * gravity * dadx;                            // gravity
      fx -=  H * H * gravity / 2.0;

      f  +=  cf * Ures * U;                                 // bottom friction

      fx +=  H * ( U*U*Dxx - 2.0*U*V*Dxy + V*V*Dyy );       // dispersion
      fy +=  H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );

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

      f   =  H * (U * dVdx  +  V * dVdy);                   // convection
      fx  =  0.0;
      fy  =  0.0;

      fx +=  H * vt * (dVdx + dUdy);                        // eddy viscosity
      fy +=  H * vt * (dVdy + dVdy);

      // -------------------------------------------------------------------------------------------
      // The following depth-averaged Boussinesq approach leads to instabilities in LES.
      // fx +=  vt * (V * dHdx + U * dHdy);
      // fy +=  vt * (V * dHdy + V * dHdy);
      // -------------------------------------------------------------------------------------------

      fx -=  H * uv;                                        // turbulence
      fy -=  H * vv;

      f  +=  H * gravity * dady;                            // gravity
      fy -=  H * H * gravity / 2.0;

      f  +=  cf * Ures * V;                                 // bottom friction

      fx +=  H * ( U*V*(Dxx-Dyy) + (U*U-V*V)*Dxy );         // dispersion
      fy +=  H * ( V*V*Dxx + 2.0*U*V*Dxy + U*U*Dyy );

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + eqidV;

      for( int j=0; j<nnd; j++ )
      {
        forcePtr[j] -= n[j] * f  +  dndx[j] * fx  +  dndy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute continuity equation (solve only for corner nodes)

      f  = H * (dUdx + dVdy) + U * dHdx + V * dHdy;
      f += SS;                                              // Source or Sink
      f *= weight;

      forcePtr = force + eqidS;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f;
      }
    }


    if( estifm )
    {
      // U-derivative of x-momentum ------------------------------------------------------

      for( int j=0; j<nnd; j++ )
      {
        double* estifmPtr = estifm[j];

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j] * weight * n[k];
        }
      }

      // V-derivative of y-momentum ------------------------------------------------------

      for( int j=0; j<nnd; j++ )
      {
        double* estifmPtr = estifm[j + eqidV] + eqidV;

        for( int k=0; k<nnd; k++ )
        {
          estifmPtr[k] += n[j] * weight * n[k];
        }
      }


      // H-derivative of continuity ------------------------------------------------------

      for( int j=0; j<ncn; j++ )
      {
        double * estifmPtr = estifm[j + eqidS] + eqidS;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j] * weight * m[k];
        }
      }
    }
  } // END of loop over all GAUSS points



  // apply transformation --------------------------------------------------------------------------

  Rotate2D( nnd, eqidV, elem->nd, estifm, force );


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

  for( int i=ncn; i<nnd; i++ )    // midside nodes
  {
    NODE* node = elem->nd[i];
    BCON* bcon = &node->bc;

    if( isFS(bcon->kind, BCON::kInlet) )
    {
      // get left and right corner node to midside node i ------------------------------------------

      int l, r;

      qShape->getCornerNodes( i, &l, &r );

      NODE* lnode = elem->nd[l];
      NODE* rnode = elem->nd[r];

      int nofeqHl = eqidS + l;    // dfH at left corner node
      int nofeqHr = eqidS + r;    // dfH at right corner node


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
      double Hl    = lnode->v.S - lnode->z;
      double Hr    = rnode->v.S - rnode->z;
      double H     = ( Hl + Hr ) / 2.0;
      double specQ = bcon->val->U;

      if( H <= 0.0 )
      {
        H     = project->hmin;
        specQ = 0.0;
      }


      // force vector ------------------------------------------------------------------------------

      if( force )  force[i] = area * (specQ - Un*H);


      // jacobi matrix -----------------------------------------------------------------------------

      if( estifm )
      {
        estifm[i][i]       = area * H;
        estifm[i][nofeqHl] = area * Un / 2.0;
        estifm[i][nofeqHr] = area * Un / 2.0;
      }
    }
  }
}

