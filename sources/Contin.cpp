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

#include "Defs.h"
#include "Report.h"
#include "Shape.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Subdom.h"


void MODEL::Continuity()
{
  static int     theFirstCall = true;
  static int     nofSect = 0;
  static double* Qx;
  static double* Qy;

  char   text[200];

# ifdef _MPI_
  static double* mpi_Qx;
  static double* mpi_Qy;
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

      if( n > nofSect )  nofSect = n;
    }


    // MPI: maximum of nofSect -----------------------------------------------------------

#   ifdef _MPI_
    nofSect = subdom->Mpi_max( nofSect );
#   endif


    // allocate dynamic memory and initialization ----------------------------------------

    if( nofSect > 0 )
    {
      Qx = new double[nofSect];
      Qy = new double[nofSect];

      if( !Qx  ||  !Qy )
        REPORT::rpt.Error( kMemoryFault, "can not allocate memory (MODEL::continuity - 2)" );

#     ifdef _MPI_
      if( subdom->npr > 1  &&  subdom->pid == 0 )
      {
        mpi_Qx = new double[nofSect];
        mpi_Qy = new double[nofSect];

        if( !mpi_Qx  ||  !mpi_Qy )
          REPORT::rpt.Error( kMemoryFault, "can not allocate memory (MODEL::continuity - 2a)" );
      }
#     endif
    }
  }

  sprintf( text, "\n\n (MODEL::Continuity)     %s at %2d sections\n",
                 "check of continuity", nofSect );
  REPORT::rpt.Output( text, 1 );


  for( int i=0; i<nofSect; i++ )
  {
    Qx[i] = Qy[i] = 0.0;
  }


  // determine discharge through continuity sections ---------------------------------------

  for( int e=0; e<control->Getne(); e++ )
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
      double* dN = qShape->dfdx[g];


      // compute normal vector at Gauss point g --------------------------------------------
      // since the normal is not reduced to unit length it
      // implies the transformation of the integrand

      double nx =  dN[0]*node[0]->y + dN[1]*node[1]->y + dN[2]*node[2]->y;
      double ny = -dN[0]*node[0]->x - dN[1]*node[1]->x - dN[2]*node[2]->x;


      double H =    M[0] * (node[0]->v.S  -  node[0]->z)
          +  M[1] * (node[1]->v.S  -  node[1]->z);

      if ( H <= 0.0 ) H = 0.0;

      double U = N[0]*node[0]->v.U + N[1]*node[1]->v.U + N[2]*node[2]->v.U;
      double V = N[0]*node[0]->v.V + N[1]*node[1]->v.V + N[2]*node[2]->v.V;

      double weight = qShape->weight[g];

      // ---------------------------------------------------------------------------------
      // sum discharge through continuity section [n-1]
      // MPI: If continuity line is on MPI-interface reduce discharge to one half

      if( isFS(node[2]->flag, NODE::kInface) )
      {
        Qx[n-1] += weight * H * U * nx / 2.0;
        Qy[n-1] += weight * H * V * ny / 2.0;
      }
      else
      {
        Qx[n-1] += weight * H * U * nx;
        Qy[n-1] += weight * H * V * ny;
      }
    }
  }


  // MPI: receive Qx and Qy on master process (pid==0) ---------------------------------------------

# ifdef _MPI_

  if( subdom->npr > 1 )
  {
    if( subdom->pid == 0 )
    {
      MPI_Status status;

      for( int s=1; s<subdom->npr; s++ )
      {
        MPI_Recv( mpi_Qx, nofSect, MPI_DOUBLE, s, 1, MPI_COMM_WORLD, &status );
        MPI_Recv( mpi_Qy, nofSect, MPI_DOUBLE, s, 2, MPI_COMM_WORLD, &status );

        for( int i=0; i<nofSect; i++ )
        {
          Qx[i] += mpi_Qx[i];
          Qy[i] += mpi_Qy[i];
        }
      }
    }
    else
    {
      MPI_Send( Qx, nofSect, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
      MPI_Send( Qy, nofSect, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD );
    }
  }

  MPI_Barrier( MPI_COMM_WORLD );

# endif


  if( subdom->pid == 0 )
  {
    sprintf( text, "\n  %s\n\n",
                   "No       total Q         Q in X          Q in Y        percent" );
    REPORT::rpt.Output( text, 1 );


    double reference = 1.0;

    for( int i=0; i<nofSect; i++ )
    {
      double Q = Qx[i] + Qy[i];

      if( i == 0 )
      {
        reference = Q;

        if( fabs (reference) < 1.0e-10)  reference = 1.0;
      }

      double percentage = 100.0 * (Q / reference);

      sprintf( text, "  %2d   %14.5le  %14.5le  %14.5le   %8.4lf\n",
                     i+1, Q, Qx[i], Qy[i], percentage );
      REPORT::rpt.Output( text, 1 );
    }
  }
}
