// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version)
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================
//
// ---------------------------------------------------------------------------------------
// functions:
//     void   Statistic()     compute statistics of correction vector
//     double GetRelax()      determine relaxation factor
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Node.h"
#include "Model.h"
#include "Subdom.h"

#include "Eqs.h"

void EQS::Update( MODEL*  model,
                  SUBDOM* subdom,
                  double* X,
                  int     ind,
                  int     varInd,
                  double* maxAbs,
                  double* maxPer,
                  double* avAbs,
                  double* avPer,
                  int*    noAbs,
                  int*    noPer )
{
  *maxAbs = 0.0;
  *maxPer = 0.0;
  *avAbs  = 0.0;
  *avPer  = 0.0;

  *noAbs = -1;
  *noPer = -1;

  int countAbs = 0;
  int countPer = 0;

  for( int n=0; n<model->np; n++ )
  {
    NODE* nd = model->node[n];

    int eqno = GetEqno( nd, ind );

    if( eqno >= 0 )
    {
      double chAbs = fabs( X[eqno] );


      // absolute changes, maximum and on average

      *avAbs += chAbs;
      countAbs++;

      if( chAbs > *maxAbs )
      {
        *maxAbs = chAbs;
        *noAbs  = nd->Getname();
      }


      // percentage changes, maximum and on average

      double H = nd->v.S - nd->z;

      double val = 0.0;

      switch( varInd )
      {
        case kVarU:   val = nd->v.U;       break;
        case kVarV:   val = nd->v.V;       break;
        case kVarH:
        case kVarS:   val = H;             break;
        case kVarK:   val = nd->v.K;       break;
        case kVarD:   val = nd->v.D;       break;
        case kVarC:   val = nd->v.C;       break;
        case kVarQb:  val = nd->v.Qb;      break;
        case kVarDz:  val = nd->dz;        break;

        case kVarUH:  val = H * nd->v.U;   break;
        case kVarVH:  val = H * nd->v.V;   break;
      }

      if( fabs(val) > 1.0e-9 )
      {
        double chPer = chAbs / fabs(val);

        *avPer += fabs(chPer);
        countPer++;

        if( fabs(chPer) > fabs(*maxPer) )
        {
          *maxPer = chPer;
          *noPer  = nd->Getname();
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////
# ifdef _MPI_

  int    mpi_noAbs;
  int    mpi_noPer;
  double mpi_maxAbs;
  double mpi_maxPer;

  if( subdom->pid == 0 )
  {
    MPI_Status status;

    for( int s=1; s<subdom->npr; s++ )
    {
      MPI_Recv( &mpi_noAbs,  1, MPI_INT,    s, 1, MPI_COMM_WORLD, &status );
      MPI_Recv( &mpi_maxAbs, 1, MPI_DOUBLE, s, 2, MPI_COMM_WORLD, &status );

      if( mpi_maxAbs > *maxAbs )
      {
        *noAbs  = mpi_noAbs;
        *maxAbs = mpi_maxAbs;
      }

      MPI_Recv( &mpi_noPer,  1, MPI_INT,    s, 1, MPI_COMM_WORLD, &status );
      MPI_Recv( &mpi_maxPer, 1, MPI_DOUBLE, s, 2, MPI_COMM_WORLD, &status );

      if( fabs(mpi_maxPer) > fabs(*maxPer) )
      {
        *noPer  = mpi_noPer;
        *maxPer = mpi_maxPer;
      }
    }
  }
  else
  {
    MPI_Send( noAbs,  1, MPI_INT,    0, 1, MPI_COMM_WORLD );
    MPI_Send( maxAbs, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD );

    MPI_Send( noPer,  1, MPI_INT,    0, 1, MPI_COMM_WORLD );
    MPI_Send( maxPer, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD );
  }

  MPI_Bcast( noAbs,  1, MPI_INT,    0, MPI_COMM_WORLD );
  MPI_Bcast( maxAbs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  MPI_Bcast( noPer,  1, MPI_INT,    0, MPI_COMM_WORLD );
  MPI_Bcast( maxPer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );


  //*maxAbs  = subdom->Mpi_max( *maxAbs );
  countAbs = subdom->Mpi_sum( countAbs );
  *avAbs   = subdom->Mpi_sum( *avAbs );
# endif
  ////////////////////////////////////////////////////////////////////////////////////////

  *avAbs /= countAbs;

  if( countPer )  *avPer /= countPer;

  *avPer  *= 100.0;
  *maxPer *= 100.0;
}
