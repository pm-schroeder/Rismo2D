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
// ---------------------------------------------------------------------------------------
// compute element stiffness matrix
//
// methods:   EQS_BL2D::Coefs()
//            EQS_BL2D::Bound()   one-dimensional boundary elements
//            EQS_BL2D::Region()  two-dimensional region elements
//
// Michael Schroeder in Mai 2005
// ---------------------------------------------------------------------------------------

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

#include "EqsBL2D.h"


int EQS_BL2D::Coefs( ELEM*    elem,
                     PROJECT* project,
                     double** estifm,
                     double*  force )
{
  if( isFS(elem->flag, ELEM::kDry) )   return false;


  // -------------------------------------------------------------------------------------

  if( isFS(coefs, kQuadratic ) )
  {
    if( isFS(elem->flag,ELEM::kBound) )
      Bound( elem, project, estifm, force, elem->GetQShape() );
    else
      Region( elem, project, estifm, force, elem->GetQShape() );
  }
  else
  {
    if( isFS(elem->flag,ELEM::kBound) )
      Bound( elem, project, estifm, force, elem->GetLShape() );
    else
      Region( elem, project, estifm, force, elem->GetLShape() );
  }

  // -------------------------------------------------------------------------------------

  return true;
}


void EQS_BL2D::Bound( ELEM*    elem,
                      PROJECT* project,
                      double** estifm,
                      double*  force,
                      SHAPE*   shape )
{
  SED* sed = &project->sed;

  int nnd = shape->nnd;

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

  if( !isFS(node[2]->bc.kind, BCON::kOutlet) )  return;


  for( int g=0; g<shape->ngp; g++ )     // loop on GAUSS points
  {
    double* n  = shape->f[g];           // quadratic shape
    double* dn = shape->dfdx[g];


    // weight of Gauss point g -----------------------------------------------------------

    double weight = shape->weight[g];


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


    // compute value of some variables at Gauss point ------------------------------------

    double U = 0.0;
    double V = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      U += n[i] * node[i]->v.U;
      V += n[i] * node[i]->v.V;
    }

    double Us = sqrt( U*U + V*V );

    double sx = 0.0;
    double sy = 0.0;

    if( Us > 1.0e-6 )
    {
      sx = U/Us;
      sy = V/Us;
    }


    // compute sediment transport at gauss point -----------------------------------------

    double Qb = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      int no = node[i]->Getno();
      Qb += n[i] * node[i]->v.Qb;
    }

    // force vector ----------------------------------------------------------------------

    double f = weight * Qb *( sx*nx + sy*ny );

    for( int i=0; i<nnd; i++ )  force[i] -= n[i] * f;


    // compute estifm --------------------------------------------------------------------

    if( estifm  &&  isFS(node[2]->bc.kind, BCON::kOutlet) )
    {
      double df = weight * ( sx*nx + sy*ny );

      for( int i=0; i<nnd; i++ )
      {
        for( int j=0; j<nnd; j++ )
        {
          estifm[i][j] += n[i] * df * n[j];
        }
      }
    }
  }

  for( int i=0; i<nnd; i++ )
  {
    NODE* nd = elem->Getnode(i);
    int   no = nd->Getno();


    // set specified flow conditions -----------------------------------------------------

    if( nd->fixqb  ||  isFS(nd->bc.kind, BCON::kRateC) )
    {
      force[i] = 0.0;

      if( estifm )
      {
        double* estifmPtr = estifm[i];
        for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
      }
    }
  }
}


void EQS_BL2D::Region( ELEM*    elem,
                       PROJECT* project,
                       double** estifm,
                       double*  force,
                       SHAPE*   shape )
{
  // -------------------------------------------------------------------------------------
  // initializations

  SED* sed = &project->sed;

  int ngp = shape->ngp;            // number of GAUSS points
  int nnd = shape->nnd;            // number of nodes

  if( force )  for( int i=0; i<maxEleq; i++ )  force[i] = 0.0;

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      double* estifmPtr = estifm[i];
      for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


  // compute coordinates relative to first node ------------------------------------------

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    x[i] = elem->nd[i]->x - *x;
    y[i] = elem->nd[i]->y - *y;
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // GAUSS point integration

  double area = 0.0;

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double  trafo[2][2];

    double* dfdxPtr = shape->dfdx[g];
    double* dfdyPtr = shape->dfdy[g];

    double  detj    = shape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );
    double  weight  = detj * shape->weight[g];

    area += weight;


    // -----------------------------------------------------------------------------------
    // compute values of specified shape function at GP g

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = shape->f[g];

    for( int i=0; i<nnd; i++ )
    {
      dndx[i] = trafo[0][0] * dfdxPtr[i] + trafo[0][1] * dfdyPtr[i];
      dndy[i] = trafo[1][0] * dfdxPtr[i] + trafo[1][1] * dfdyPtr[i];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP g

    // integrate with specified shape ----------------------------------------------------

    double sx  = 0.0;
    double sy  = 0.0;
    double Qb  = 0.0;
    double Qbc = 0.0;
    double PLs = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      NODE* node = elem->nd[i];
      int   no   = node->Getno();

      sx    += n[i] * sed->sx[no];
      sy    += n[i] * sed->sy[no];

      Qb    += n[i] * node->v.Qb;
      Qbc   += n[i] * sed->qbc[no];
      PLs   += n[i] * sed->PLs[no];
    }


    // -----------------------------------------------------------------------------------
    // compute right side
    // -----------------------------------------------------------------------------------

    double f, fx, fy;

    f  =  weight * (Qb - Qbc) * PLs;

    fx = -weight * Qb * sx;
    fy = -weight * Qb * sy;

    for( int i=0; i<nnd; i++ )
    {
      force[i] -= n[i] * f  +  dndx[i] * fx  +  dndy[i] * fy;
    }


    // -----------------------------------------------------------------------------------
    // compute components of NEWTON-RAPHSON Jacobi matrix
    // -----------------------------------------------------------------------------------

    if( estifm )
    {
      double df__, dfx_, dfy_;

      df__  =  weight * PLs;
      dfx_  = -weight * sx;
      dfy_  = -weight * sy;

      double t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];

      for( int i=0; i<nnd; i++ )
      {
        t[i]  = df__ * n[i];
        tx[i] = dfx_ * n[i];
        ty[i] = dfy_ * n[i];
      }

      for( int i=0; i<nnd; i++ )
      {
        double* estifmPtr = estifm[i];

        for( int j=0; j<nnd; j++ )
        {
          estifmPtr[j] += n[i]*t[j] + dndx[i]*tx[j] + dndy[i]*ty[j];
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////

  for( int i=0; i<nnd; i++ )
  {
    NODE* nd = elem->Getnode(i);
    int   no = nd->Getno();

    if( nd->fixqb )
    {
      force[i] = 0.0;

      if( estifm )
      {
        double* estifmPtr = estifm[i];
        for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;

        estifmPtr[i] = area;
      }
    }

    // set specified flow conditions -----------------------------------------------------

    else if( isFS(nd->bc.kind, BCON::kRateC) )
    {
      // compute transport on inlet boundary -----------------------------------------------

      double qbin = 0.0;

      if( nd->bc.val->Qb[0] < 0.0 )
      {
         qbin = fabs(nd->bc.val->Qb[0]) * sed->qbc[no];
      }
      else
      {
        if( nd->bc.val->Qb[0] <= sed->qbc[no] )  qbin = nd->bc.val->Qb[0];
        else                                     qbin = sed->qbc[no];
      }

      force[i] = area * (qbin - nd->v.Qb);

      if( estifm )
      {
        double* estifmPtr = estifm[i];
        for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;

        estifmPtr[i] = area;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
/*
void EQS_BL2D::BoundDiffusion( ELEM*    elem,
                               PROJECT* project,
                               double** estifm,
                               double*  force,
                               SHAPE*   shape )
{
  SED* sed = &project->sed;

  int nnd = shape->nnd;

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

  if( !isFS(node[2]->bc->kind, BCON::kOutlet) )  return;


  for( int g=0; g<shape->ngp; g++ )     // loop on GAUSS points
  {
    double* n  = shape->f[g];           // quadratic shape
    double* dn = shape->dfdx[g];


    // weight of Gauss point g -----------------------------------------------------------

    double weight = shape->weight[g];


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


    // compute value of some variables at Gauss point ------------------------------------

    double U = 0.0;
    double V = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      U += n[i] * node[i]->v.U;
      V += n[i] * node[i]->v.V;
    }

    double Us = sqrt( U*U + V*V );

    double sx = 0.0;
    double sy = 0.0;

    if( Us > 1.0e-6 )
    {
      sx = U/Us;
      sy = V/Us;
    }


    // compute sediment transport at gauss point -----------------------------------------

    double qb = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      int no = node[i]->Getno();
      qb += n[i] * node[i]->v.qb;
    }

    // force vector ----------------------------------------------------------------------

    double f = weight * qb *( sx*nx + sy*ny );

    for( int i=0; i<nnd; i++ )  force[i] -= n[i] * f;


    // compute estifm --------------------------------------------------------------------

    if( estifm  &&  isFS(node[2]->bc->kind, BCON::kOutlet) )
    {
      double df = weight * Us * ( sx*nx + sy*ny );

      for( int i=0; i<nnd; i++ )
      {
        for( int j=0; j<nnd; j++ )
        {
          estifm[i][j] += n[i] * df * n[j];
        }
      }
    }
  }

  for( int i=0; i<nnd; i++ )
  {
    NODE* nd = elem->Getnode(i);
    int   no = nd->Getno();


    // set specified flow conditions -----------------------------------------------------

    if( nd->fixqb  ||  (nd->bc && isFS(nd->bc->kind, BCON::kFlowC)) )
    {
      force[i] = 0.0;

      if( estifm )
      {
        double* estifmPtr = estifm[i];
        for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
      }
    }
  }
}


void EQS_BL2D::RegionDiffusion( ELEM*    elem,
                                PROJECT* project,
                                double** estifm,
                                double*  force,
                                SHAPE*   shape )
{
  // -------------------------------------------------------------------------------------
  // initializations

  SED* sed = &project->sed;

  int ngp = shape->ngp;            // number of GAUSS points
  int nnd = shape->nnd;            // number of nodes

  if( force )  for( int i=0; i<maxEleq; i++ )  force[i] = 0.0;

  if( estifm )
  {
    for( int i=0; i<maxEleq; i++ )
    {
      double* estifmPtr = estifm[i];
      for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;
    }
  }


  // compute coordinates relative to first node ------------------------------------------

  double x[kMaxNodes2D], y[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    x[i] = elem->nd[i]->x - *x;
    y[i] = elem->nd[i]->y - *y;
  }

  x[0] = y[0] = 0.0;


  // -------------------------------------------------------------------------------------
  // GAUSS point integration

  double area = 0.0;

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with linear shape functions

    double  trafo[2][2];

    double* dfdxPtr = shape->dfdx[g];
    double* dfdyPtr = shape->dfdy[g];

    double  detj    = shape->jacobi2D( nnd, dfdxPtr, dfdyPtr, x, y, trafo );
    double  weight  = detj * shape->weight[g];

    area += weight;


    // -----------------------------------------------------------------------------------
    // compute values of specified shape function at GP g

    double dndx[kMaxNodes2D], dndy[kMaxNodes2D];

    double* n = shape->f[g];

    for( int i=0; i<nnd; i++ )
    {
      dndx[i] = trafo[0][0] * dfdxPtr[i] + trafo[0][1] * dfdyPtr[i];
      dndy[i] = trafo[1][0] * dfdxPtr[i] + trafo[1][1] * dfdyPtr[i];
    }


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP g

    // integrate with specified shape ----------------------------------------------------

    double sx    = 0.0;
    double sy    = 0.0;
    double qb    = 0.0;
    double qbr   = 0.0;
    double PLL   = 0.0;

    double Exx   = 0.0;
    double Exy   = 0.0;
    double Eyy   = 0.0;

    double U     = 0.0;
    double V     = 0.0;
    double dCbdx = 0.0;
    double dCbdy = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      NODE* node = elem->nd[i];
      int   no   = node->Getno();

      sx  += n[i] * sed->sx[no];
      sy  += n[i] * sed->sy[no];

      qb  += n[i] * node->v.qb;
      qbr += n[i] * sed->qbr[no];
      PLL += n[i] * sed->iLs[no];

      Exx += n[i] * node->exx * node->vt;
      Exy += n[i] * node->exy * node->vt;
      Eyy += n[i] * node->eyy * node->vt;

      U   += n[i] * node->v.U;
      V   += n[i] * node->v.V;

      double Us = sqrt(node->v.U*node->v.U + node->v.V*node->v.V);

      if( Us > 1.0e-6 )
      {
        dCbdx += dndx[i] * (node->v.qb / Us);
        dCbdy += dndy[i] * (node->v.qb / Us);
      }
    }

    double Us = sqrt( U*U + V*V );

    // -----------------------------------------------------------------------------------
    // compute right side
    // -----------------------------------------------------------------------------------

    double f, fx, fy;

    f   =  weight * (qb - qbr) * PLL;

    fx  = -weight * qb * sx;                                // advection
    fy  = -weight * qb * sy;

    fx +=  weight * (Exx * dCbdx   +  Exy * dCbdy);         // diffusion
    fy +=  weight * (Exy * dCbdx   +  Eyy * dCbdy);

    for( int i=0; i<nnd; i++ )
    {
      force[i] -= n[i] * f  +  dndx[i] * fx  +  dndy[i] * fy;
    }


    // -----------------------------------------------------------------------------------
    // compute components of NEWTON-RAPHSON Jacobi matrix
    // -----------------------------------------------------------------------------------

    if( estifm )
    {
      double df__, dfx_, dfy_, dfxx, dfxy, dfyx, dfyy;

      df__  =  weight * Us * PLL;
      dfx_  = -weight * Us * sx;
      dfy_  = -weight * Us * sy;

      dfxx  =  weight * Exx;
      dfyy  =  weight * Eyy;
      dfxy  =
      dfyx  =  weight * Exy;

      double t[kMaxNodes2D], tx[kMaxNodes2D], ty[kMaxNodes2D];

      for( int i=0; i<nnd; i++ )
      {
        t[i]  = df__ * n[i];
        tx[i] = dfx_ * n[i]  +  dfxx * dndx[i]  +  dfxy * dndy[i];
        ty[i] = dfy_ * n[i]  +  dfyx * dndx[i]  +  dfyy * dndy[i];
      }

      for( int i=0; i<nnd; i++ )
      {
        double* estifmPtr = estifm[i];

        for( int j=0; j<nnd; j++ )
        {
          estifmPtr[j] += n[i]*t[j] + dndx[i]*tx[j] + dndy[i]*ty[j];
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////

  for( int i=0; i<nnd; i++ )
  {
    NODE* nd = elem->Getnode(i);
    int   no = nd->Getno();

    if( nd->fixqb )
    {
      force[i] = 0.0;

      if( estifm )
      {
        double* estifmPtr = estifm[i];
        for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;

        estifmPtr[i] = area;
      }
    }

    // set specified flow conditions -----------------------------------------------------

    else if( nd->bc  &&  isFS(nd->bc->kind, BCON::kFlowC) )
    {
      // compute transport on inlet boundary -----------------------------------------------

      double qbin = 0.0;

      if( nd->bc->spec.qb < 0.0 )
      {
         qbin = fabs(nd->bc->spec.qb) * sed->qbe[no];
      }
      else
      {
        if( nd->bc->spec.qb <= sed->qbe[no] )  qbin = nd->bc->spec.qb;
        else                                   qbin = sed->qbr[no];
      }

      force[i] = area * (qbin - nd->v.qb);

      if( estifm )
      {
        double* estifmPtr = estifm[i];
        for( int j=0; j<maxEleq; j++ )  estifmPtr[j] = 0.0;

        double Us = sqrt( nd->v.U*nd->v.U + nd->v.V*nd->v.V );
        estifmPtr[i] = area * Us;
      }
    }
  }
}
*/
