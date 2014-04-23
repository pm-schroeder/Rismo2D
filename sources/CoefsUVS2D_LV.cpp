// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_LV
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

#include "EqsUVS2D_LV.h"


#define kFiniteVolume  1


int EQS_UVS2D_LV::Coefs( ELEM*    elem,
                         PROJECT* project,
                         double** estifm,
                         double*  force )
{
  if( isFS(elem->flag, ELEM::kDry) )    return 0;
  if( isFS(elem->flag, ELEM::kBound) )  return 0;

  if( isFS(project->actualTurb, BCONSET::kVtAnisotrop) )
  {
    RegionAI( elem, project, estifm, force );
  }
  else
  {
    if( elem->shape == kTri )  RegionBT( elem, project, estifm, force );
    else                       Region( elem, project, estifm, force );
  }

  return 1;
}


void EQS_UVS2D_LV::Bound( ELEM*    elem,
                          PROJECT* project,
                          double** estifm,
                          double*  force,
                          int      r0U,
                          int      r0V,
                          int      r1U,
                          int      r1V,
                          int      cH )
{
  SHAPE*  lShape  = elem->GetLShape();

  TYPE* type = TYPE::Getid( elem->type );

  double  gravity = project->g;


  int solidFlag;

  if(    !isFS(elem->nd[2]->bc.kind, BCON::kSlip)
      ||  isFS(elem->flag, ELEM::kInlet)  ||  isFS(elem->flag, ELEM::kOutlet) )
  {
    solidFlag = false;
  }

  else
  {
    solidFlag = true;
  }


  NODE* node[3];

  node[0] = elem->nd[0];
  node[1] = elem->nd[1];
  node[2] = elem->nd[2];


  // -------------------------------------------------------------------------------------
  // determine water surface elevation

  double S = node[2]->el[0]->P;            // in adjacent region element

  if( isFS(node[2]->bc.kind, BCON::kOutlet) )
  {
    // use experimental downstream flow depth ...
    S = (node[0]->bc.val->S + node[1]->bc.val->S)/2.0;
  }


  for( int g=0; g<lShape->ngp; g++ )       // loop on GAUSS points
  {
    double* m  = lShape->f[g];             // linear shape
    double* dm = lShape->dfdx[g];


    // -----------------------------------------------------------------------------------
    // compute flow depth H

    double H =    m[0] * (node[0]->v.S - node[0]->z)
               +  m[1] * (node[1]->v.S - node[1]->z);

    if( H <= 0.0 )  H = project->hmin;


    // -----------------------------------------------------------------------------------
    // compute normal vector
    // since the normal is not reduced to unit length it
    // implies the transformation of the integrand

    double nx =  dm[0]*node[0]->y + dm[1]*node[1]->y;
    double ny = -dm[0]*node[0]->x - dm[1]*node[1]->x;

    double len = sqrt( nx*nx + ny*ny );


    // -----------------------------------------------------------------------------------
    // weight of Gauss point g

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

      cf = 0.0;//m[0]*node[0]->cfw + m[1]*node[1]->cfw;
    }


    if( force )
    {
      double fU = weight * gravity * H * S * nx;
      double fV = weight * gravity * H * S * ny;


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
      double dfUdU, dfUdV, dfUdH;
      double dfVdU, dfVdV, dfVdH;


      if( !isFS(node[2]->bc.kind,BCON::kOutlet) )
      {
        dfUdH = weight * gravity * H * nx;
        dfVdH = weight * gravity * H * ny;
      }

      else
      {
        dfUdH = 0.0;
        dfVdH = 0.0;
      }


      if( solidFlag )
      {
        double iUres;

        if ( Ures > 1.0e-8 )  iUres = 1.0 / Ures;
        else                  iUres = 0.0;

        dfUdU  = weight * len * cf * H * (U * U * iUres  +  Ures);
        dfUdV  = weight * len * cf * H * (U * V * iUres);

        dfVdU  = weight * len * cf * H * (U * V * iUres);
        dfVdV  = weight * len * cf * H * (V * V * iUres  +  Ures);

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
      }

      estifm[r0U][cH] += m[0] * dfUdH;
      estifm[r0V][cH] += m[0] * dfVdH;

      estifm[r1U][cH] += m[1] * dfUdH;
      estifm[r1V][cH] += m[1] * dfVdH;
    }
  }
}


void EQS_UVS2D_LV::Region( ELEM*    elem,
                           PROJECT* project,
                           double** estifm,
                           double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  TYPE* type = TYPE::Getid( elem->type );

  int ngp = qShape->ngp;         // number of GAUSS points
  int nnd = qShape->nnd;         // number of nodes in all
  int ncn = lShape->nnd;         // number of corner nodes

  double gravity = project->g;
  double S       = elem->P;

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
      for( int j=0; j<maxEleq; j++ )  estifm[i][j] = 0.0;
    }
  }


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double xmin, xmax, ymin, ymax;
  double x[kMaxNodes2D], y[kMaxNodes2D];

  xmin = xmax = x[0] = elem->nd[0]->x;
  ymin = ymax = y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    if( elem->nd[i]->x < xmin )  xmin = elem->nd[i]->x;
    if( elem->nd[i]->x > xmax )  xmax = elem->nd[i]->x;

    if( elem->nd[i]->y < ymin )  ymin = elem->nd[i]->y;
    if( elem->nd[i]->y > ymax )  ymax = elem->nd[i]->y;

    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;

  double dl = sqrt( (xmax-xmin)*(ymax-ymin) );
  double ls = type->lm * dl;


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

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j] + trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j] + trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at the GP i

    double dHdt = 0.0;

    double H    = 0.0;
    double dHdx = 0.0;
    double dHdy = 0.0;

    double U    = 0.0;
    double dUdt = 0.0;
    double dUdx = 0.0;
    double dUdy = 0.0;

    double V    = 0.0;
    double dVdt = 0.0;
    double dVdx = 0.0;
    double dVdy = 0.0;

    double cf   = 0.0;

    double uu   = 0.0;
    double uv   = 0.0;
    double vv   = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE*  node = elem->nd[j];

      double ndZ  = node->z;
      double ndH  = node->v.S - ndZ;
      double ndU  = node->v.U;
      double ndV  = node->v.V;

      dHdt +=    m[j] * node->v.dSdt;

      H    +=    m[j] * ndH;
      dHdx += dmdx[j] * ndH;
      dHdy += dmdy[j] * ndH;

      U    +=    m[j] * ndU;
      dUdx += dmdx[j] * ndU;
      dUdy += dmdy[j] * ndU;
      dUdt +=    m[j] * node->v.dUdt;

      V    +=    m[j] * ndV;
      dVdx += dmdx[j] * ndV;
      dVdy += dmdy[j] * ndV;
      dVdt +=    m[j] * node->v.dVdt;

      cf   +=    m[j] * node->cf;

      uu   +=    m[j] * node->uu;
      uv   +=    m[j] * node->uv;
      vv   +=    m[j] * node->vv;
    }


    if( H <= 0.0 )  H = project->hmin;


    // -----------------------------------------------------------------------------------

    double vt = 0.0;

    if( isFS(project->actualTurb, BCONSET::kVtConstant) )
    {
      vt = type->vt;
    }
    else
    {
      for( int j=0; j<ncn; j++ )
      {
        vt += m[j] * elem->nd[j]->vt;
      }
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < type->vt )  vt = type->vt;
    }


    // add eddy viscosity to kinematic viscosity -----------------------------------------

    vt += project->vk;


    // -----------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double Ures = sqrt( U*U + V*V );           // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;


      // ---------------------------------------------------------------------------------
      // compute x-momentum

      f   = H * dUdt;                          // time

      f  += H * (U * dUdx  +  V * dUdy);       // convection

      fx  = H * vt * (dUdx + dUdx);            // eddy viscosity
      fy  = H * vt * (dUdy + dVdx);

      fx -= H * uu;                            // turbulence
      fy -= H * uv;

      f  -= gravity * dHdx * S;                // gravity
      fx -= gravity * H * S;

      f  += cf * Ures * U;                     // bottom friction

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f  +  dmdx[j] * fx  +  dmdy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute y-momentum

      f   = H * dVdt;                          // time

      f  += H * (U * dVdx  +  V * dVdy);       // convection

      fx  = H * vt * (dVdx + dUdy);            // eddy viscosity
      fy  = H * vt * (dVdy + dVdy);

      fx -= H * uv;                            // turbulence
      fy -= H * vv;

      f  -= gravity * dHdy * S;                // gravity
      fy -= gravity * H * S;

      f  += cf * Ures * V;                     // bottom friction

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + startV;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f  +  dmdx[j] * fx  +  dmdy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute continuity

      f  = dHdt; //dSdt;

#     ifndef kFiniteVolume

      f += H * (dUdx + dVdy)  +  U * dHdx  +  V * dHdy;

#     endif

      f *= weight;

      force[startS] -= f;
    }


    if( estifm )
    {
      double* estifmPtr;
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfx_, dfy_, dfxx, dfxy, dfyx, dfyy;


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

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] = dfxx * dmdx[j];
        ty[j] = dfyy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j];

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // --- V-derivative of x-momentum --------------------------------------------------

      df__  =  weight * H * dUdy;
      dfyx  =  weight * H * vt;
      df__ +=  weight * cf * iUres * U * V;

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        ty[j] = dfyx * dmdx[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j] + startV;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdy[j]*ty[k];
        }
      }


      // --- S-derivative of x-momentum --------------------------------------------------

      df__  = -weight * gravity * dHdx;
      dfx_  = -weight * gravity * H;

      for( int j=0; j<ncn; j++ )
      {
        estifm[j][startS] += m[j]*df__ + dmdx[j]*dfx_;
      }


      // --- U-derivative of y-momentum --------------------------------------------------

      df__  =  weight * H * dVdx;
      dfxy  =  weight * H * vt;
      df__ +=  weight * cf * iUres * U * V;

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
        tx[j] = dfxy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startV];

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k];
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

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] = dfxx * dmdx[j];
        ty[j] = dfyy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startV] + startV;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // --- S-derivative of y-momentum --------------------------------------------------

      df__  = -weight * gravity * dHdy;
      dfy_  = -weight * gravity * H;

      for( int j=0; j<ncn; j++ )
      {
        estifm[j + startV][startS] += m[j]*df__ + dmdy[j]*dfy_;
      }


#     ifndef kFiniteVolume

      // --- U-derivative of continuity --------------------------------------------------

      df__ = weight * dHdx;
      df_x = weight * H;

      estifmPtr = estifm[startS];

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr[j] += df__*m[j] + df_x*dmdx[j];
      }


      // --- V-derivative of continuity --------------------------------------------------

      df__ = weight * dHdy;
      df_y = weight * H;


      estifmPtr = estifm[startS] + startV;

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr[j] += df__*m[j] + df_y*dmdy[j];
      }

#     endif


      // --- S-derivative of continuity --------------------------------------------------

      df__  = weight * relaxThdt_H;

      estifm[startS][startS] += df__;
    }
  }


# ifdef kFiniteVolume

  // -------------------------------------------------------------------------------------
  // compute flow through element boundary

  SHAPE* lineShape = SHAPE::get( kLine, 2 );

  for( int i=0; i<ncn; i++ )
  {
    int j = (i+1) % ncn;

    NODE* ndi = elem->nd[i];
    NODE* ndj = elem->nd[j];

    double Hi = ndi->v.S - ndi->z;
    double Ui = ndi->v.U;
    double Vi = ndi->v.V;

    double Hj = ndj->v.S - ndj->z;
    double Uj = ndj->v.U;
    double Vj = ndj->v.V;

    for( int g=0; g<lineShape->ngp; g++ )
    {
      double* m      = lineShape->f[g];
      double* dm     = lineShape->dfdx[g];
      double  weight = lineShape->weight[g];

      double nx =  dm[0]*y[i] + dm[1]*y[j];
      double ny = -dm[0]*x[i] - dm[1]*x[j];

      double Hg = m[0] * Hi  +  m[1] * Hj;
      double Ug = m[0] * Ui  +  m[1] * Uj;
      double Vg = m[0] * Vi  +  m[1] * Vj;

      if( force )
      {
        // -------------------------------------------------------------------------------
        // continuity

        force[startS] -= weight * Hg * (Ug*nx + Vg*ny);
      }


      if( estifm )
      {
        // -------------------------------------------------------------------------------
        // continuity

        // U-derivative of continuity
        estifm[startS][i] += weight * Hg * nx * m[0];
        estifm[startS][j] += weight * Hg * nx * m[1];

        // V-derivative of continuity
        estifm[startS][i + startV] += weight * Hg * ny * m[0];
        estifm[startS][j + startV] += weight * Hg * ny * m[1];
      }
    }
  }

# endif


  // -------------------------------------------------------------------------------------
  // add estifm and force for adjacent boundary elements

  for( int i=ncn; i<nnd; i++ )
  {
    NODE* node = elem->nd[i];

    if( isFS(node->flag, NODE::kBound) )
    {
      ELEM* bdel = project->M2D->Getbound(node->Getno());

      // row index for force vector and element stiffness matrix
      // and  column index for element stiffness matrix

      int r0U =          i-ncn;   int r1U =          (i-ncn + 1)%ncn;
      int r0V = startV + i-ncn;   int r1V = startV + (i-ncn + 1)%ncn;

      int cH  = startS;

      Bound( bdel, project, estifm, force, r0U, r0V, r1U, r1V, cH );
    }
  }


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
        estifm[i][startS] = area * Un;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // eliminate equation S (partial Gauss elimination)

  double* SPtr = estifm[startS];

  for( int i=0; i<ncn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in x-momentum equation --------------------------------------------------

    MPtr = estifm[i];

    factor = MPtr[startS] / SPtr[startS];

    force[i] -= factor * force[startS];

    for( int j=0; j<ncn; j++ )
    {
      MPtr[j]          -= factor * SPtr[j];
      MPtr[j + startV] -= factor * SPtr[j + startV];
    }

    MPtr[startS] = 0.0;


    // eliminate in y-momentum equation --------------------------------------------------

    MPtr = estifm[i + startV];

    factor = MPtr[startS] / SPtr[startS];

    force[i+startV] -= factor * force[startS];

    for( int j=0; j<ncn; j++ )
    {
      MPtr[j]          -= factor * SPtr[j];
      MPtr[j + startV] -= factor * SPtr[j + startV];
    }

    MPtr[startS] = 0.0;
  }


  // -------------------------------------------------------------------------------------
  // save elimination equation

  double* PElimEqPtr = PElimEq[elem->Getno()];

  for( int i=0; i<ncn; i++ )
  {
    PElimEqPtr[i]     = SPtr[i];
    PElimEqPtr[i+ncn] = SPtr[i + startV];
  }

  PElimEqPtr[2*ncn]   = SPtr[startS];
  PElimEqPtr[2*ncn+1] = force[startS];
}


//////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_LV::RegionBT( ELEM*    elem,
                             PROJECT* project,
                             double** estifm,
                             double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  SHAPE* lShape = elem->GetLShape();
  SHAPE* bShape = elem->GetBShape();
  SHAPE* qShape = elem->GetQShape();

  TYPE* type = TYPE::Getid( elem->type );

  int ngp = lShape->ngp;         // number of GAUSS points

  int nnd = qShape->nnd;         // number of nodes
  int ncn = lShape->nnd;         // number of corner nodes

  int nbn = bShape->nnd;         // number of bubble nodes (for bubble shape function)

  double gravity = project->g;
  double P       = elem->P;

  int bid[4] = { 0, 1, 2, 6 };

  int uid[4] = {   0,     1,     2,     3 }; //2*nnd   };
  int vid[4] = { nnd, nnd+1, nnd+2, nnd+3 }; //2*nnd+1 };
  int sid    = {                    2*nnd }; //2*nnd+2 };

  if( isFS(elem->nd[3]->flag, NODE::kRotat) )  CF( elem->nd[3]->flag, NODE::kRotat );

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

  double xmin, xmax, ymin, ymax;
  double x[kMaxNodes2D], y[kMaxNodes2D];

  xmin = xmax = x[0] = elem->nd[0]->x;
  ymin = ymax = y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    if( elem->nd[i]->x < xmin )  xmin = elem->nd[i]->x;
    if( elem->nd[i]->x > xmax )  xmax = elem->nd[i]->x;

    if( elem->nd[i]->y < ymin )  ymin = elem->nd[i]->y;
    if( elem->nd[i]->y > ymax )  ymax = elem->nd[i]->y;

    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;

  double dl = sqrt( (xmax-xmin)*(ymax-ymin) );
  double ls = type->lm * dl;


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

    double* b = bShape->f[g];

    double  dbdx[kMaxNodes2D], dbdy[kMaxNodes2D];

    for( int i=0; i<nbn; i++ )
    {
      int k = bid[i];

      double dfdx = bShape->dfdx[g][k];
      double dfdy = bShape->dfdy[g][k];

      dbdx[k] = trafo[0][0] * dfdx + trafo[0][1] * dfdy;
      dbdy[k] = trafo[1][0] * dfdx + trafo[1][1] * dfdy;
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP i

    double dHdt = 0.0;

    double H    = 0.0;
    double dHdx = 0.0;
    double dHdy = 0.0;

    double cf   = 0.0;

    double uu   = 0.0;
    double uv   = 0.0;
    double vv   = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE*  node = elem->nd[j];

      double ndZ  = node->z;
      double ndH  = node->v.S - ndZ;

      dHdt +=    m[j] * node->v.dSdt;

      H    +=    m[j] * ndH;
      dHdx += dmdx[j] * ndH;
      dHdy += dmdy[j] * ndH;

      cf   +=    m[j] * node->cf;

      uu   +=    m[j] * node->uu;
      uv   +=    m[j] * node->uv;
      vv   +=    m[j] * node->vv;
    }

    if( H <= 0.0 )  H = project->hmin;


    // -----------------------------------------------------------------------------------

    double U    = 0.0;
    double dUdt = 0.0;
    double dUdx = 0.0;
    double dUdy = 0.0;

    double V    = 0.0;
    double dVdt = 0.0;
    double dVdx = 0.0;
    double dVdy = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE*  node = elem->nd[j];

      double ndU  = node->v.U;
      double ndV  = node->v.V;

      U    +=    b[j] * ndU;
      dUdx += dbdx[j] * ndU;
      dUdy += dbdy[j] * ndU;
      dUdt +=    b[j] * node->v.dUdt;

      V    +=    b[j] * ndV;
      dVdx += dbdx[j] * ndV;
      dVdy += dbdy[j] * ndV;
      dVdt +=    b[j] * node->v.dVdt;
    }

    U    +=    b[6] * elem->U;
    dUdx += dbdx[6] * elem->U;
    dUdy += dbdy[6] * elem->U;
    //dUdt +=    b[6] * elem->dUdt;

    V    +=    b[6] * elem->V;
    dVdx += dbdx[6] * elem->V;
    dVdy += dbdy[6] * elem->V;
    //dVdt +=    b[6] * elem->dVdt;


    // -----------------------------------------------------------------------------------

    double vt = 0.0;

    if( isFS(project->actualTurb, BCONSET::kVtConstant) )
    {
      vt = type->vt;
    }
    else
    {
      for( int j=0; j<ncn; j++ )
      {
        vt += m[j] * elem->nd[j]->vt;
      }
    }

    if( isFS(project->actualTurb, BCONSET::kVtMin) )
    {
      if( vt < type->vt )  vt = type->vt;
    }


    // add eddy viscosity to kinematic viscosity -----------------------------------------

    vt += project->vk;


    // -----------------------------------------------------------------------------------
    // compute UVH-equation and coefficients of NEWTON-RAPHSON matrix

    double Ures = sqrt( U*U + V*V );           // absolute velocity


    if( force )
    {
      double  f, fx, fy;


      // ---------------------------------------------------------------------------------
      // compute x-momentum

      f   = H * dUdt;                          // time

      f  += H * (U * dUdx  +  V * dUdy);       // convection

      fx  = H * vt * (dUdx + dUdx);            // eddy viscosity
      fy  = H * vt * (dUdy + dVdx);

      fx -= H * uu;                            // turbulence
      fy -= H * uv;

      f  -= gravity * dHdx * P;                // gravity
      fx -= gravity * H * P;

      f  += cf * Ures * U;                     // bottom friction

      f  *= weight;
      fx *= weight;
      fy *= weight;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];
        force[uid[i]] -= b[k] * f  +  dbdx[k] * fx  +  dbdy[k] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute y-momentum

      f   = H * dVdt;                          // time

      f  += H * (U * dVdx  +  V * dVdy);       // convection

      fx  = H * vt * (dVdx + dUdy);            // eddy viscosity
      fy  = H * vt * (dVdy + dVdy);

      fx -= H * uv;                            // turbulence
      fy -= H * vv;

      f  -= gravity * dHdy * P;                // gravity
      fy -= gravity * H * P;

      f  += cf * Ures * V;                     // bottom friction

      f  *= weight;
      fx *= weight;
      fy *= weight;


      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];
        force[vid[i]] -= b[k] * f  +  dbdx[k] * fx  +  dbdy[k] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute continuity

      f  = dHdt; //dSdt;

#     ifndef kFiniteVolume

      f += H * (dUdx + dVdy)  +  U * dHdx  +  V * dHdy;

#     endif

      f *= weight;

      force[sid] -= f;
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

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];
        t[i]  = df__ * b[k]  +  df_x * dbdx[k]  +  df_y * dbdy[k];
        tx[i] = dfxx * dbdx[k];
        ty[i] = dfyy * dbdy[k];
      }

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        for( int j=0; j<nbn; j++ )
        {
          estifm[uid[i]][uid[j]] += b[k]*t[j] + dbdx[k]*tx[j] + dbdy[k]*ty[j];
        }
      }


      // --- V-derivative of x-momentum --------------------------------------------------

      df__  =  weight * H * dUdy;
      dfyx  =  weight * H * vt;
      df__ +=  weight * cf * iUres * U * V;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];
        t[i]  = df__ * b[k];
        ty[i] = dfyx * dbdx[k];
      }

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        for( int j=0; j<nbn; j++ )
        {
          estifm[uid[i]][vid[j]] += b[k]*t[j] + dbdy[k]*ty[j];
        }
      }


      // --- S-derivative of x-momentum --------------------------------------------------

      df__  = -weight * gravity * dHdx;
      dfx_  = -weight * gravity * H;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        estifm[uid[i]][sid] += b[k]*df__ + dbdx[k]*dfx_;
      }


      // --- U-derivative of y-momentum --------------------------------------------------

      df__  =  weight * H * dVdx;
      dfxy  =  weight * H * vt;
      df__ +=  weight * cf * iUres * U * V;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];
        t[i]  = df__ * b[k];
        tx[i] = dfxy * dbdy[k];
      }

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        for( int j=0; j<nbn; j++ )
        {
          estifm[vid[i]][uid[j]] += b[k]*t[j] + dbdx[k]*tx[j];
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

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];
        t[i]  = df__ * b[k]  +  df_x * dbdx[k]  +  df_y * dbdy[k];
        tx[i] = dfxx * dbdx[k];
        ty[i] = dfyy * dbdy[k];
      }

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        for( int j=0; j<nbn; j++ )
        {
          estifm[vid[i]][vid[j]] += b[k]*t[j] + dbdx[k]*tx[j] + dbdy[k]*ty[j];
        }
      }


      // --- S-derivative of y-momentum --------------------------------------------------

      df__  = -weight * gravity * dHdy;
      dfy_  = -weight * gravity * H;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        estifm[vid[i]][sid] += b[k]*df__ + dbdy[k]*dfy_;
      }


#     ifndef kFiniteVolume

      // --- U-derivative of continuity --------------------------------------------------

      df__ = weight * dHdx;
      df_x = weight * H;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        estifm[sid][uid[i]] += df__*b[k] + df_x*dbdx[k];
      }


      // --- V-derivative of continuity --------------------------------------------------

      df__ = weight * dHdy;
      df_y = weight * H;

      for( int i=0; i<nbn; i++ )
      {
        int k = bid[i];

        estifm[sid][vid[i]] += df__*b[k] + df_y*dbdy[k];
      }

#     endif


      // --- S-derivative of continuity --------------------------------------------------

      df__  = weight * relaxThdt_H;

      estifm[sid][sid] += df__;
    }
  }


# ifdef kFiniteVolume

  // -------------------------------------------------------------------------------------
  // compute flow through element boundary

  SHAPE* lineShape = SHAPE::get( kLine, 2 );

  for( int i=0; i<ncn; i++ )
  {
    int j = (i+1) % ncn;

    NODE* ndi = elem->nd[i];
    NODE* ndj = elem->nd[j];

    double Hi = ndi->v.S - ndi->z;
    double Ui = ndi->v.U;
    double Vi = ndi->v.V;

    double Hj = ndj->v.S - ndj->z;
    double Uj = ndj->v.U;
    double Vj = ndj->v.V;

    for( int g=0; g<lineShape->ngp; g++ )
    {
      double* m      = lineShape->f[g];
      double* dm     = lineShape->dfdx[g];
      double  weight = lineShape->weight[g];

      double nx =  dm[0]*y[i] + dm[1]*y[j];
      double ny = -dm[0]*x[i] - dm[1]*x[j];

      double Hg = m[0] * Hi  +  m[1] * Hj;
      double Ug = m[0] * Ui  +  m[1] * Uj;
      double Vg = m[0] * Vi  +  m[1] * Vj;

      if( force )
      {
        // -------------------------------------------------------------------------------
        // continuity

        force[sid] -= weight * Hg * (Ug*nx + Vg*ny);
      }


      if( estifm )
      {
        // -------------------------------------------------------------------------------
        // continuity

        // U-derivative of continuity
        estifm[sid][uid[i]] += weight * Hg * nx * m[0];
        estifm[sid][uid[j]] += weight * Hg * nx * m[1];

        // V-derivative of continuity
        estifm[sid][vid[i]] += weight * Hg * ny * m[0];
        estifm[sid][vid[j]] += weight * Hg * ny * m[1];
      }
    }
  }

# endif


  // -------------------------------------------------------------------------------------
  // add estifm and force for adjacent boundary elements
  // note: this is not done from method Coefs() like it is done in other equations
  //       (see e.g. EQS_UVS2D). this is necessary because the element equation S
  //       will be eliminated within this method (see below).

  for( int i=ncn; i<nnd; i++ )
  {
    NODE* node = elem->nd[i];

    if( isFS(node->flag, NODE::kBound) )
    {
      ELEM* bdel = project->M2D->Getbound(node->Getno());

      // row index for force vector and element stiffness matrix
      // and  column index for element stiffness matrix

      int r0U = uid[i-ncn];   int r1U = (i-ncn + 1)%ncn;
      int r0V = vid[i-ncn];   int r1V = vid[(i-ncn + 1)%ncn];

      int cH  = sid;

      Bound( bdel, project, estifm, force, r0U, r0V, r1U, r1V, cH );
    }
  }


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
        estifm[i][i]   = area * H;
        estifm[i][sid] = area * Un;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // eliminate surface equation S for element centers (partial Gauss elimination)

  double* SPtr = estifm[sid];

  for( int i=0; i<nbn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in x-momentum equation

    MPtr = estifm[uid[i]];

    factor = MPtr[sid] / SPtr[sid];

    force[uid[i]] -= factor * force[sid];

    for( int j=0; j<nbn; j++ )
    {
      MPtr[uid[j]] -= factor * SPtr[uid[j]];
      MPtr[vid[j]] -= factor * SPtr[vid[j]];
    }

    MPtr[sid] = 0.0;


    // eliminate in y-momentum equation

    MPtr = estifm[vid[i]];

    factor = MPtr[sid] / SPtr[sid];

    force[vid[i]] -= factor * force[sid];

    for( int j=0; j<nbn; j++ )
    {
      MPtr[uid[j]] -= factor * SPtr[uid[j]];
      MPtr[vid[j]] -= factor * SPtr[vid[j]];
    }

    MPtr[sid] = 0.0;
  }


  // save elimination equation -----------------------------------------------------------

  double* PElimEqPtr = PElimEq[elem->Getno()];

  for( int i=0; i<nbn; i++ )
  {
    PElimEqPtr[i]     = SPtr[uid[i]];
    PElimEqPtr[i+nbn] = SPtr[vid[i]];
  }

  PElimEqPtr[2*nbn]   = SPtr[sid];      // diagonal entry
  PElimEqPtr[2*nbn+1] = force[sid];     // RHS


  // -------------------------------------------------------------------------------------
  // eliminate momentum equation V for element centers (partial Gauss elimination)

  double* VPtr = estifm[vid[3]];

  for( int i=0; i<nbn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in x-momentum equation

    MPtr = estifm[uid[i]];

    factor = MPtr[vid[3]] / VPtr[vid[3]];

    force[uid[i]] -= factor * force[vid[3]];

    for( int j=0; j<nbn; j++ )
    {
      MPtr[uid[j]] -= factor * VPtr[uid[j]];
      MPtr[vid[j]] -= factor * VPtr[vid[j]];
    }

    MPtr[vid[3]] = 0.0;
  }

  for( int i=0; i<ncn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in y-momentum equation

    MPtr = estifm[vid[i]];

    factor = MPtr[vid[3]] / VPtr[vid[3]];

    force[vid[i]] -= factor * force[vid[3]];

    for( int j=0; j<nbn; j++ )
    {
      MPtr[uid[j]] -= factor * VPtr[uid[j]];
      MPtr[vid[j]] -= factor * VPtr[vid[j]];
    }

    MPtr[vid[3]] = 0.0;
  }


  // save elimination equation -----------------------------------------------------------

  double* VElimEqPtr = VElimEq[elem->Getno()];

  for( int i=0; i<nbn; i++ )
  {
    VElimEqPtr[i]     = VPtr[uid[i]];
    VElimEqPtr[i+nbn] = VPtr[vid[i]];
  }

  VElimEqPtr[2*nbn]   = VPtr[vid[3]];   // diagonal entry
  VElimEqPtr[2*nbn+1] = force[vid[3]];  // RHS


  // -------------------------------------------------------------------------------------
  // eliminate momentum equation U for element centers (partial Gauss elimination)

  double* UPtr = estifm[uid[3]];

  for( int i=0; i<ncn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in x-momentum equation

    MPtr = estifm[uid[i]];

    factor = MPtr[uid[3]] / UPtr[uid[3]];

    force[uid[i]] -= factor * force[uid[3]];

    for( int j=0; j<ncn; j++ )
    {
      MPtr[uid[j]] -= factor * UPtr[uid[j]];
      MPtr[vid[j]] -= factor * UPtr[vid[j]];
    }

    MPtr[uid[3]] = 0.0;


    // eliminate in y-momentum equation

    MPtr = estifm[vid[i]];

    factor = MPtr[uid[3]] / UPtr[uid[3]];

    force[vid[i]] -= factor * force[uid[3]];

    for( int j=0; j<ncn; j++ )
    {
      MPtr[uid[j]] -= factor * UPtr[uid[j]];
      MPtr[vid[j]] -= factor * UPtr[vid[j]];
    }

    MPtr[uid[3]] = 0.0;
  }


  // save elimination equation -----------------------------------------------------------

  double* UElimEqPtr = UElimEq[elem->Getno()];

  for( int i=0; i<nbn; i++ )
  {
    UElimEqPtr[i]     = UPtr[uid[i]];
    UElimEqPtr[i+nbn] = UPtr[vid[i]];
  }

  UElimEqPtr[2*nbn]   = UPtr[uid[3]];   // diagonal entry
  UElimEqPtr[2*nbn+1] = force[uid[3]];  // RHS
}


//////////////////////////////////////////////////////////////////////////////////////////

void EQS_UVS2D_LV::RegionAI( ELEM*    elem,
                             PROJECT* project,
                             double** estifm,
                             double*  force )
{
  // -------------------------------------------------------------------------------------
  // initializations

  SHAPE* lShape = elem->GetLShape();
  SHAPE* qShape = elem->GetQShape();

  TYPE* type = TYPE::Getid( elem->type );

  int ngp = qShape->ngp;         // number of GAUSS points
  int nnd = qShape->nnd;         // number of nodes in all
  int ncn = lShape->nnd;         // number of corner nodes

  double gravity = project->g;
  double S       = elem->P;

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
      for( int j=0; j<maxEleq; j++ )  estifm[i][j] = 0.0;
    }
  }


  // -------------------------------------------------------------------------------------
  // compute coordinates relative to first node

  double xmin, xmax, ymin, ymax;
  double x[kMaxNodes2D], y[kMaxNodes2D];

  xmin = xmax = x[0] = elem->nd[0]->x;
  ymin = ymax = y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    if( elem->nd[i]->x < xmin )  xmin = elem->nd[i]->x;
    if( elem->nd[i]->x > xmax )  xmax = elem->nd[i]->x;

    if( elem->nd[i]->y < ymin )  ymin = elem->nd[i]->y;
    if( elem->nd[i]->y > ymax )  ymax = elem->nd[i]->y;

    x[i] = elem->nd[i]->x - x[0];
    y[i] = elem->nd[i]->y - y[0];
  }

  x[0] = y[0] = 0.0;

  double dl = sqrt( (xmax-xmin)*(ymax-ymin) );
  double ls = type->lm * dl;


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

    for( int j=0; j<ncn; j++ )
    {
      dmdx[j] = trafo[0][0] * dfdxPtr[j] + trafo[0][1] * dfdyPtr[j];
      dmdy[j] = trafo[1][0] * dfdxPtr[j] + trafo[1][1] * dfdyPtr[j];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at the GP i

    double dHdt = 0.0;

    double H    = 0.0;
    double dHdx = 0.0;
    double dHdy = 0.0;

    double U    = 0.0;
    double dUdt = 0.0;
    double dUdx = 0.0;
    double dUdy = 0.0;

    double V    = 0.0;
    double dVdt = 0.0;
    double dVdx = 0.0;
    double dVdy = 0.0;

    double cf   = 0.0;

    double uu   = 0.0;
    double uv   = 0.0;
    double vv   = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE*  node = elem->nd[j];

      double ndZ  = node->z;
      double ndH  = node->v.S - ndZ;
      double ndU  = node->v.U;
      double ndV  = node->v.V;

      dHdt +=    m[j] * node->v.dSdt;

      H    +=    m[j] * ndH;
      dHdx += dmdx[j] * ndH;
      dHdy += dmdy[j] * ndH;

      U    +=    m[j] * ndU;
      dUdx += dmdx[j] * ndU;
      dUdy += dmdy[j] * ndU;
      dUdt +=    m[j] * node->v.dUdt;

      V    +=    m[j] * ndV;
      dVdx += dmdx[j] * ndV;
      dVdy += dmdy[j] * ndV;
      dVdt +=    m[j] * node->v.dVdt;

      cf   +=    m[j] * node->cf;

      uu   +=    m[j] * node->uu;
      uv   +=    m[j] * node->uv;
      vv   +=    m[j] * node->vv;
    }


    if( H <= 0.0 )  H = project->hmin;


    // compute eddy viscosity according to ELDER's assumption ------------------------------

    double vtxx = 0.0;
    double vtxy = 0.0;
    double vtyy = 0.0;

    for( int j=0; j<ncn; j++ )
    {
      NODE* node = elem->nd[j];

      vtxx += m[j] * node->exx * node->vt;
      vtxy += m[j] * node->exy * node->vt;
      vtyy += m[j] * node->eyy * node->vt;
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

    double Ures = sqrt( U*U + V*V );           // absolute velocity


    if( force )
    {
      double  f, fx, fy;
      double* forcePtr;


      // ---------------------------------------------------------------------------------
      // compute x-momentum

      f   = H * dUdt;                          // time

      f  += H * (U * dUdx  +  V * dUdy);       // convection

      fx  = H * (vtxx * dUdx + vtxy * dUdy);   // eddy viscosity (approximation)
      fy  = H * (vtxy * dUdx + vtyy * dUdy);

      fx -= H * uu;                            // turbulence
      fy -= H * uv;

      f  -= gravity * dHdx * S;                // gravity
      fx -= gravity * H * S;

      f  += cf * Ures * U;                     // bottom friction

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f  +  dmdx[j] * fx  +  dmdy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute y-momentum

      f   = H * dVdt;                          // time

      f  += H * (U * dVdx  +  V * dVdy);       // convection

      fx  = H * (vtxx * dVdx + vtxy * dVdy);   // eddy viscosity (approximation)
      fy  = H * (vtxy * dVdx + vtyy * dVdy);

      fx -= H * uv;                            // turbulence
      fy -= H * vv;

      f  -= gravity * dHdy * S;                // gravity
      fy -= gravity * H * S;

      f  += cf * Ures * V;                     // bottom friction

      f  *= weight;
      fx *= weight;
      fy *= weight;

      forcePtr = force + startV;

      for( int j=0; j<ncn; j++ )
      {
        forcePtr[j] -= m[j] * f  +  dmdx[j] * fx  +  dmdy[j] * fy;
      }


      // ---------------------------------------------------------------------------------
      // compute continuity

      f  = dHdt; //dSdt;

#     ifndef kFiniteVolume

      f += H * (dUdx + dVdy)  +  U * dHdx  +  V * dHdy;

#     endif

      f *= weight;

      force[startS] -= f;
    }


    if( estifm )
    {
      double* estifmPtr;
      double  t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];
      double  df__, df_x, df_y, dfx_, dfy_, dfxx, dfxy, dfyx, dfyy;


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

      dfxx  =  weight * H * vtxx;
      dfxy  =  weight * H * vtxy;

      dfyx  =  weight * H * vtxy;
      dfyy  =  weight * H * vtyy;

      df__ +=  weight * cf * (iUres * U*U  +  Ures);

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j]  +  dfxy * dmdy[j];
        ty[j] =                 dfyx * dmdx[j]  +  dfyy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j];

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // --- V-derivative of x-momentum --------------------------------------------------

      df__  =  weight * H * dUdy;
      df__ +=  weight * cf * iUres * U * V;

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j] + startV;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }


      // --- S-derivative of x-momentum --------------------------------------------------

      df__  = -weight * gravity * dHdx;
      dfx_  = -weight * gravity * H;

      for( int j=0; j<ncn; j++ )
      {
        estifm[j][startS] += m[j]*df__ + dmdx[j]*dfx_;
      }


      // --- U-derivative of y-momentum --------------------------------------------------

      df__  =  weight * H * dVdx;
      df__ +=  weight * cf * iUres * U * V;

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startV];

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k];
        }
      }


      // --- V-derivative of y-momentum --------------------------------------------------

      df__  =  weight * H * relaxThdt_UV;

      df__ +=  weight * H * dVdy;
      df_x  =  weight * H * U;
      df_y  =  weight * H * V;

      dfxx  =  weight * H * vtxx;
      dfxy  =  weight * H * vtxy;

      dfyx  =  weight * H * vtxy;
      dfyy  =  weight * H * vtyy;

      df__ +=  weight * cf * (iUres * V*V  +  Ures);

      for( int j=0; j<ncn; j++ )
      {
        t[j]  = df__ * m[j]  +  df_x * dmdx[j]  +  df_y * dmdy[j];
        tx[j] =                 dfxx * dmdx[j]  +  dfxy * dmdy[j];
        ty[j] =                 dfyx * dmdx[j]  +  dfyy * dmdy[j];
      }

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr = estifm[j + startV] + startV;

        for( int k=0; k<ncn; k++ )
        {
          estifmPtr[k] += m[j]*t[k] + dmdx[j]*tx[k] + dmdy[j]*ty[k];
        }
      }


      // --- S-derivative of y-momentum --------------------------------------------------

      df__  = -weight * gravity * dHdy;
      dfy_  = -weight * gravity * H;

      for( int j=0; j<ncn; j++ )
      {
        estifm[j + startV][startS] += m[j]*df__ + dmdy[j]*dfy_;
      }


#     ifndef kFiniteVolume

      // --- U-derivative of continuity --------------------------------------------------

      df__ = weight * dHdx;
      df_x = weight * H;

      estifmPtr = estifm[startS];

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr[j] += df__*m[j] + df_x*dmdx[j];
      }


      // --- V-derivative of continuity --------------------------------------------------

      df__ = weight * dHdy;
      df_y = weight * H;


      estifmPtr = estifm[startS] + startV;

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr[j] += df__*m[j] + df_y*dmdy[j];
      }

#     endif


      // --- S-derivative of continuity --------------------------------------------------

      df__  = weight * relaxThdt_H;

      estifm[startS][startS] += df__;
    }
  }


# ifdef kFiniteVolume

  // -------------------------------------------------------------------------------------
  // compute flow through element boundary

  SHAPE* lineShape = SHAPE::get( kLine, 2 );

  for( int i=0; i<ncn; i++ )
  {
    int j = (i+1) % ncn;

    NODE* ndi = elem->nd[i];
    NODE* ndj = elem->nd[j];

    double Hi = ndi->v.S - ndi->z;
    double Ui = ndi->v.U;
    double Vi = ndi->v.V;

    double Hj = ndj->v.S - ndj->z;
    double Uj = ndj->v.U;
    double Vj = ndj->v.V;

    for( int g=0; g<lineShape->ngp; g++ )
    {
      double* m      = lineShape->f[g];
      double* dm     = lineShape->dfdx[g];
      double  weight = lineShape->weight[g];

      double nx =  dm[0]*y[i] + dm[1]*y[j];
      double ny = -dm[0]*x[i] - dm[1]*x[j];

      double Hg = m[0] * Hi  +  m[1] * Hj;
      double Ug = m[0] * Ui  +  m[1] * Uj;
      double Vg = m[0] * Vi  +  m[1] * Vj;

      if( force )
      {
        // -------------------------------------------------------------------------------
        // continuity

        force[startS] -= weight * Hg * (Ug*nx + Vg*ny);
      }


      if( estifm )
      {
        // -------------------------------------------------------------------------------
        // continuity

        // U-derivative of continuity
        estifm[startS][i] += weight * Hg * nx * m[0];
        estifm[startS][j] += weight * Hg * nx * m[1];


        // V-derivative of continuity
        estifm[startS][i + startV] += weight * Hg * ny * m[0];
        estifm[startS][j + startV] += weight * Hg * ny * m[1];
      }
    }
  }

# endif


  // -------------------------------------------------------------------------------------
  // add estifm and force for adjacent boundary elements

  for( int i=ncn; i<nnd; i++ )
  {
    NODE* node = elem->nd[i];

    if( isFS(node->flag, NODE::kBound) )
    {
      ELEM* bdel = project->M2D->Getbound(node->Getno());

      // row index for force vector and element stiffness matrix
      // and  column index for element stiffness matrix

      int r0U =          i-ncn;   int r1U =          (i-ncn + 1)%ncn;
      int r0V = startV + i-ncn;   int r1V = startV + (i-ncn + 1)%ncn;

      int cH  = startS;

      Bound( bdel, project, estifm, force, r0U, r0V, r1U, r1V, cH );
    }
  }


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
        estifm[i][startS] = area * Un;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // eliminate equation S (partial Gauss elimination)

  double* SPtr = estifm[startS];

  for( int i=0; i<ncn; i++ )
  {
    double  factor;
    double* MPtr;

    // eliminate in x-momentum equation --------------------------------------------------

    MPtr = estifm[i];

    factor = MPtr[startS] / SPtr[startS];

    force[i] -= factor * force[startS];

    for( int j=0; j<ncn; j++ )
    {
      MPtr[j]          -= factor * SPtr[j];
      MPtr[j + startV] -= factor * SPtr[j + startV];
    }

    MPtr[startS] = 0.0;


    // eliminate in y-momentum equation --------------------------------------------------

    MPtr = estifm[i + startV];

    factor = MPtr[startS] / SPtr[startS];

    force[i+startV] -= factor * force[startS];

    for( int j=0; j<ncn; j++ )
    {
      MPtr[j]          -= factor * SPtr[j];
      MPtr[j + startV] -= factor * SPtr[j + startV];
    }

    MPtr[startS] = 0.0;
  }


  // -------------------------------------------------------------------------------------
  // save elimination equation

  double* PElimEqPtr = PElimEq[elem->Getno()];

  for( int i=0; i<ncn; i++ )
  {
    PElimEqPtr[i]     = SPtr[i];
    PElimEqPtr[i+ncn] = SPtr[i + startV];
  }

  PElimEqPtr[2*ncn]   = SPtr[startS];
  PElimEqPtr[2*ncn+1] = force[startS];
}
