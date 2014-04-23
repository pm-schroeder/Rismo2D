// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_Dz
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

#include "EqsBL2D.h"


int EQS_DZ::Coefs( ELEM*    elem,
                   PROJECT* project,
                   double** estifm,
                   double*  force )
{
  if( isFS(elem->flag, ELEM::kDry) ) return false;

  // -------------------------------------------------------------------------------------

  if( coefs == kBottomEvol )
  {
    if( isFS(elem->flag,ELEM::kBound) )
      Bound( elem, project, estifm, force, elem->GetQShape() );
    else
      Region( elem, project, estifm, force, elem->GetQShape() );
  }
  else if( coefs == kBottomEvol_L )
  {
    if( isFS(elem->flag,ELEM::kBound) )
      Bound( elem, project, estifm, force, elem->GetLShape() );
    else
      Region( elem, project, estifm, force, elem->GetLShape() );
  }

  // -------------------------------------------------------------------------------------

  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////

void EQS_DZ::Bound( ELEM*    elem,
                    PROJECT* project,
                    double** estifm,
                    double*  force,
                    SHAPE*   qShape )
{
  SED*   sed    = &project->sed;
  SHAPE* lShape = elem->GetLShape();

  int ncn = lShape->nnd;
  int nnd = qShape->nnd;

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


  if(    !isFS(node[2]->bc.kind, BCON::kInlet)
      && !isFS(node[2]->bc.kind, BCON::kOutlet) )  return;


  for( int g=0; g<qShape->ngp; g++ )     // loop on GAUSS points
  {
    double* m  = lShape->f[g];           // linear shape
    double* n  = qShape->f[g];           // specified (quadratic) shape
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

    // compute sediment transport --------------------------------------------------------

    double Qb = 0.0;

    //if( isFS(node[2]->bc->kind, BCON::kFlowC) )
    //{
    //  // use experimental upstream sediment transport ... --------------------------------
    //  for( int i=0; i<nnd; i++ )  qb += n[i] * node[i]->bc->spec.qb;
    //}
    //else
    if( isFS(node[2]->bc.kind, BCON::kOutlet) )
    {
     // ... or compute sediment transport at gauss point ---------------------------------
      for( int i=0; i<nnd; i++ )  Qb += n[i] * node[i]->v.Qb;
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


    // force vector ----------------------------------------------------------------------

    if( force )
    {
      double f = weight * Qb *( sx*nx + sy*ny );

      for( int i=0; i<ncn; i++ )  force[i] -= m[i] * f;
    }
  }
}


void EQS_DZ::Region( ELEM*    elem,
                     PROJECT* project,
                     double** estifm,
                     double*  force,
                     SHAPE*   qShape )
{
  // -------------------------------------------------------------------------------------
  // initializations

  SED*   sed    = &project->sed;
  SHAPE* lShape = elem->GetLShape();

  double por    = 1.0 - project->sed.por;


  int ncn = lShape->nnd;
  int nnd = qShape->nnd;            // number of nodes
  int ngp = qShape->ngp;            // number of GAUSS points

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

  double x[kMaxNodes2D], y[kMaxNodes2D], adz[kMaxNodes2D];

  x[0] = elem->nd[0]->x;
  y[0] = elem->nd[0]->y;

  for( int i=1; i<nnd; i++ )
  {
    x[i]   = elem->nd[i]->x - *x;
    y[i]   = elem->nd[i]->y - *y;
  }

  x[0] = y[0] = 0.0;

  for( int i=0; i<nnd; i++ )
  {
    int no = elem->nd[i]->Getno();
    adz[i] = etaQb[no];
  }

  // -------------------------------------------------------------------------------------
  // GAUSS point integration

  double area = 0.0;

  for( int g=0; g<ngp; g++ )
  {
    // -----------------------------------------------------------------------------------
    // form JACOBIAN transformation matrix with shape functions

    double  trafo[2][2];

    double* dfdxPtr = lShape->dfdx[g];
    double* dfdyPtr = lShape->dfdy[g];

    double  detj    = lShape->jacobi2D( ncn, dfdxPtr, dfdyPtr, x, y, trafo );
    double  weight  = detj * lShape->weight[g];

    area += weight;


    // -----------------------------------------------------------------------------------
    // compute values of linear shape function at GP g

    double dmdx[kMaxNodes2D], dmdy[kMaxNodes2D];

    double* m = lShape->f[g];

    for( int i=0; i<ncn; i++ )
    {
      dmdx[i] = trafo[0][0] * dfdxPtr[i] + trafo[0][1] * dfdyPtr[i];
      dmdy[i] = trafo[1][0] * dfdxPtr[i] + trafo[1][1] * dfdyPtr[i];
    }


    // -----------------------------------------------------------------------------------
    // compute values of specified shape function at GP g

    double* n = qShape->f[g];


    // -----------------------------------------------------------------------------------
    // compute flow parameters and their derivatives at GP g

    // integrate with linear shape -------------------------------------------------------

    double dzdt = 0.0;

    for( int i=0; i<ncn; i++ )
    {
      NODE* node = elem->nd[i];
      int   no   = node->Getno();

      dzdt += m[i] * this->dzdt[no];
    }


    // integrate with specified shape ----------------------------------------------------

    double sx = 0.0;
    double sy = 0.0;
    double Qb = 0.0;

    for( int i=0; i<nnd; i++ )
    {
      NODE* node = elem->nd[i];
      int   no   = node->Getno();

      double ndsx = sed->sx[no];
      sx     +=    n[i] * ndsx;

      double ndsy = sed->sy[no];
      sy     +=    n[i] * ndsy;

      double ndqb = node->v.Qb;
      Qb     +=    n[i] * adz[i] * ndqb;
    }


    // erosion/deposition due to mud exchange with water body ----------------------------

    double EDs = 0.0;

    if( sed->M > 0.0 )
    {
      double H = 0.0;
      double M = 0.0;

      for( int i=0; i<ncn; i++ )
      {
        NODE* nd = elem->nd[i];

        H += m[i] * (nd->v.S - nd->z);

        if( nd->zor > nd->zero )  M += m[i] * sed->M;
      }

      if( H < 0.0 )  H = 0.0;

      double U = 0.0;
      double V = 0.0;
      double C = 0.0;

      for( int i=0; i<ncn; i++ )
      {
        NODE* nd = elem->nd[i];

        U += n[i] * nd->v.U;
        V += n[i] * nd->v.V;
        C += n[i] * nd->v.C;
      }

      double Us   = sqrt( U*U + V*V );

      double Utau = sed->GetUtau( Us, H, 0.0, project );
      double tau  = project->rho * Utau * Utau;

      double Es = 0.0;     // erosion
      double Ds = 0.0;     // deposition

      if( sed->tauc > 0.0  &&  tau > sed->tauc  &&  M > 0.0 )
        Es = M * (tau/sed->tauc - 1.0);

      if( sed->taus > 0.0  &&  tau < sed->taus  &&  C > 0.0 )
        Ds = sed->us * (1.0 - tau/sed->taus) * C;

      EDs += (Es - Ds) / sed->rhob;
    }


    // -----------------------------------------------------------------------------------
    // compute right side
    // -----------------------------------------------------------------------------------

    double f  =  weight * (por*dzdt - EDs);
    double fx = -weight * Qb * sx;
    double fy = -weight * Qb * sy;

    for( int i=0; i<ncn; i++ )
    {
      force[i] -= m[i] * f  +  dmdx[i] * fx  +  dmdy[i] * fy;
    }


    // -----------------------------------------------------------------------------------
    // compute components of NEWTON-RAPHSON Jacobi matrix
    // -----------------------------------------------------------------------------------

    double df__ =  weight * por;

    for( int i=0; i<ncn; i++ )
    {
      double* estifmPtr = estifm[i];

      for( int j=0; j<ncn; j++ )
      {
        estifmPtr[j] +=  m[i] * df__ * m[j];
      }
    }
  }
}
