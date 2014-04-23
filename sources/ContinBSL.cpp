// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class MODEL
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

// ---------------------------------------------------------------------------------------
// determine bed and suspended transport through control lines
// ---------------------------------------------------------------------------------------
//
#include "Defs.h"
#include "Report.h"
#include "Shape.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Subdom.h"


void MODEL::ContinuityBSL()
{
  static int     theFirstCall = true;
  static int     nofSect;
  static double* C;                          // concentration in control sections

  char   text[200];

# ifdef _MPI_
  static double* mpi_C;
#endif


  // determine number of continuity sections on the first call ---------------------------

  if( theFirstCall )
  {
    theFirstCall = false;
    nofSect = -1;

    for( int e=0; e<control->Getne(); e++ )
    {
      ELEM* el = control->Getelem(e);

      int n = el->type;

      if ( n <= 0 )
        REPORT::rpt.Error( "number of continuity section < 1 (MODEL::ContinuityBSL - 1)" );

      if( n > nofSect )  nofSect = n;
    }


    // MPI: maximum of nofSect -----------------------------------------------------------

#   ifdef _MPI_
    nofSect = subdom->Mpi_max( nofSect );
#   endif


    // allocate dynamic memory and initialization ----------------------------------------

    if( nofSect > 0 )
    {
      C = new double[nofSect];

      if( !C )
        REPORT::rpt.Error( kMemoryFault, "can not allocate memory (MODEL::ContinuityBSL - 2)" );
    }
  }

  for( int i=0; i<nofSect; i++ )  C[i] = 0.0;


  // determine concentration at continuity sections --------------------------------------

  for( int e=0; e<control->Getne(); e++ )             // loop on control elements
  {
    ELEM* el = control->Getelem(e);

    int n = el->type;

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
      double* dN = qShape->dfdx[g];      // derivatives of shape fkt


      // compute normal vector at Gauss point g ------------------------------------------
      // since the normal is not reduced to unit length it
      // implies the transformation of the integrand

      double nx =  dN[0]*node[0]->y + dN[1]*node[1]->y + dN[2]*node[2]->y;
      double ny = -dN[0]*node[0]->x - dN[1]*node[1]->x - dN[2]*node[2]->x;


      // compute concentration Cgp and depth Hgh at Gauss point g ------------------------
      // using shape functions

      double Cgp = N[0]*node[0]->v.C + N[1]*node[1]->v.C + N[2]*node[2]->v.C;
      double Ugp = N[0]*node[0]->v.U + N[1]*node[1]->v.U + N[2]*node[2]->v.U;
      double Vgp = N[0]*node[0]->v.V + N[1]*node[1]->v.V + N[2]*node[2]->v.V;
      double Hgp = M[0] * (node[0]->v.S  -  node[0]->z)
                 + M[1] * (node[1]->v.S  -  node[1]->z);

      if( Hgp < 0.0 ) Hgp = 0.0;


      // Gauss point integration ---------------------------------------------------------

      double weight = qShape->weight[g];
      C[n-1] += weight * ( Ugp*nx +Vgp*ny ) * Cgp * Hgp;
    }
  }

  // MPI: broadcasting of C --------------------------------------------------------------

# ifdef _MPI_

  if( subdom->npr > 1 )
  {
    if( subdom->pid == 0 )
    {
      MPI_Status status;

      for( int s=1; s<subdom->npr; s++ )
      {
        MPI_Recv( mpi_C, nofSect, MPI_DOUBLE, s, 1, MPI_COMM_WORLD, &status );

        for( int i=0; i<nofSect; i++ )  C[i] += mpi_C[i];
      }
    }
    else
    {
      MPI_Send( C, nofSect, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
    }
  }

  MPI_Barrier( MPI_COMM_WORLD );

# endif


  sprintf( text, "\n\n (ContinuityBSL)   %s at %2d sections\n\n  %s\n\n",
                 "check of sediment continuity", nofSect,
                 "No      total load    percent     " );
  REPORT::rpt.Output( text, 1 );


  double reference = C[0];
  if( fabs(reference) < 1.0e-10 )  reference = 1.0;

  for( int i=0; i<nofSect; i++ )
  {
    double percentage = 100.0 * (C[i] / reference);

    sprintf( text, "  %2d   %12.3le   %7.3lf\n", i+1, C[i], percentage );
    REPORT::rpt.Output( text, 1 );
  }
}
