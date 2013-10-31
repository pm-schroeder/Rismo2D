// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
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
//  distributed version of flexible GMRES.
//  Zhongze Li, Dec. 2000
// ---------------------------------------------------------------------------------------
//
//  ON ENTRY
// ==========
//     crsm  =  matrix in compact row storage
//   precon  =  a pointer to the preconditioner
//        b  =  right hand side vector
//
//  ON RETURN
// ===========
//        x  =  solution
//
//  Zhongze Li, Dec. 2000
// ---------------------------------------------------------------------------------------
//   code changes
//
//   date        author
//   ----------  --------------
//   2004-08-14  P.M. Schroeder  adaption for use in Finite ELement Program Rismo2D
//
//                               changed variable names ----------------------------------
//                               n    -> neq                    | number of equations
//                               sol  -> x                      | solution
//                               rhs  -> b                      | right hand side
// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Project.h"
#include "Memory.h"
#include "CRSMat.h"
#include "Eqs.h"
#include "Precon.h"

#include "P_fgmresd.h"


#define ZERO      0.0
#define EPSILON   1.0e-20
#define EPSMAC    1.0e-16


P_FGMRESD::P_FGMRESD()
{
  solverType = kParmsFgmresd;
}


P_FGMRESD::~P_FGMRESD()
{
}


int P_FGMRESD::Iterate( PROJECT* project, CRSMAT* crsm, double* b, double* x,
                        PRECON* precon )
{
  // -------------------------------------------------------------------------------------
  // PART 1:  Initialization.
  // -------------------------------------------------------------------------------------

  // retrieve information from input data structure
  int im     = this->mkyrl;                // dimension of Krylov subspace

  int neq    = crsm->m_neq;                // dimension of subdomain
  int neq_dn = crsm->m_neq_dn;

  double eps = this->maxDiff;              // precision of computing

  // allocate memory for working local arrays
  double** vv = new double*[im+1];
  for( int i=0; i<im+1; i++ )
  {
    vv[i] = new double[neq];
  }

  double** hh = new double*[im];
  for( int i=0; i<im; i++ )
  {
    hh[i] = new double[im+1];
  }

  double** w = new double*[im];
  for( int i=0; i<im; i++ )
  {
    w[i] = new double[neq];
  }

  double* c     = (double*) MEMORY::memo.Array_eq( neq );
  double* s     = (double*) MEMORY::memo.Array_eq( neq );
  double* rs    = (double*) MEMORY::memo.Array_eq( neq );
  double* hhloc = (double*) MEMORY::memo.Array_eq( neq );

  double* xs    = (double*) MEMORY::memo.Array_eq( neq );   // best solution up to iters


  // -------------------------------------------------------------------------------------
  // PART 2:  Iteration.
  // -------------------------------------------------------------------------------------

  this->accuracy = 1.0;

  double alpha      = -1.0;
  int    out_flag   = true;
  int    convergent = true;
  int    iters      = 0;      // counter for iterations
  int    itac       = 0;      // iters where the best solution occured

  double eps1   = 0.0;
  double resid  = 0.0;
  double relacc = 0.0;

  // outer loop starts here
  while( out_flag )
  {
    // vv[0] = A*x
    crsm->MulVec( x, vv[0], project, eqs );

    for( int j=0; j<neq; j++ )
    {
      vv[0][j] = b[j] - vv[0][j];
    }

    double ro = ddot( neq_dn, vv[0], vv[0] );
#   ifdef _MPI_
    ro = project->subdom.Mpi_sum( ro );
#   endif

    ro = sqrt( ro );

    if( fabs(ro-ZERO) <= EPSILON )
    {
      convergent = true;
      out_flag   = false;
      break;
    }

    double t = 1.0 / ro;
    dscal( neq, t, vv[0] );

    if( iters == 0 )
    {
      resid  = ro;
      relacc = 1.0;
      eps1   = eps * ro;

      // output the original residual norm -----------------------------------------------
      REPORT::rpt.Message( 2, "\n (P_FGMRESD::Iterate)    original L2-Norm ||r|| = %10.4le\n", ro );
    }

    // initialize first term of rhs of hessenberg system
    rs[0] = ro;

    int i = -1;

    int in_flag = true;

    while( in_flag )
    {
      i++;
      iters++;
      int i1 = i + 1;

      // matrix * vector product / preconditioning operation, w[i] = M^{-1} * vv[i]
      if( precon )
      {
        precon->Solve( project, eqs, vv[i], w[i] );
        crsm->MulVec( w[i], vv[i1], project, eqs );
      }
      else
      {
        crsm->MulVec( vv[i], vv[i1], project, eqs );
      }

      // classical gram - schmidt
      for( int j=0; j<=i; j++ )
      {
        hhloc[j] = ddot( neq_dn, vv[j], vv[i1] );

#       ifdef _MPI_
        hh[i][j] = project->subdom.Mpi_sum( hhloc[j] );
#       endif
      }

      for( int j=0; j<=i; j++ )
      {
        alpha = -hh[i][j];
        daxpy( neq, alpha, vv[j], vv[i1] );
      }

      t = ddot( neq_dn, vv[i1], vv[i1] );
#     ifdef _MPI_
      t = project->subdom.Mpi_sum( t );
#     endif

      t = sqrt( t );

      hh[i][i1] = t;

      if( fabs(t-ZERO) > EPSILON )
      {
        t = 1.0 / t;
        dscal( neq, t, vv[i1] );
      }

      // done with classical gram schimd and arnoldi step.
      // now update factorization of hh
      if( i != 0 )
      {
        for( int k=1; k<=i; k++ )
        {
          int k1 = k - 1;

          t = hh[i][k1];

          hh[i][k1] =  c[k1] * t  +  s[k1] * hh[i][k];
          hh[i][k]  = -s[k1] * t  +  c[k1] * hh[i][k];
        }
      }

      double gam = sqrt( hh[i][i]*hh[i][i] + hh[i][i1]*hh[i][i1] );
      if( fabs(gam-ZERO) <= EPSILON )  gam = EPSMAC;

      // determine-next-plance-rotation
      c[i]   = hh[i][i]/gam;
      s[i]   = hh[i][i1]/gam;
      rs[i1] = -s[i]*rs[i];
      rs[i]  =  c[i]*rs[i];

      // determine res. norm and test for convergence
      hh[i][i] = c[i]*hh[i][i] + s[i]*hh[i][i1];

      ro = fabs( rs[i1] );
#     ifdef _MPI_
      ro = project->subdom.Mpi_max( ro );
#     endif

      if( (i+1 >= im) || (ro <= eps1) || iters >= maxIter )
      {
        in_flag = false;
      }

      relacc = ro / resid;

#     ifdef kIteratCount
      REPORT::rpt.Screen( "\r (P_FGMRESD::Iterate)    %d. %s = %10.4le",
              iters, "Iteration | accuracy", relacc );
      REPORT::rpt.PrintTime( -1 );
#     endif

      // ---------------------------------------------------------------------------------
      // check for the best solution up to now

      //if( relacc < accuracy )
      //{
        itac = iters;
        accuracy = relacc;
      //}
    }

    // now compute solution, first solve upper triangular system
    rs[i] = rs[i]/hh[i][i];

    for( int ii=2; ii<=i+1; ii++ )
    {
      int k  = i - ii + 1;
      int k1 = k + 1;

      t = rs[k];
      for( int j=k1; j<=i; j++ )
      {
        t -= hh[j][k] * rs[j];
      }
      rs[k] = t/hh[k][k];
    }

    // done with back substitution. now form linear combination to get solution
    for( int j=0; j<=i; j++ )
    {
      t = rs[j];
      daxpy( neq, t, w[j], x );
    }

    // -----------------------------------------------------------------------------------
    // save the best solution up to now

    //if( iters == itac )
      dcopy( neq, x, xs );


    // -----------------------------------------------------------------------------------
    // test for return

    if( ro <= eps1 )
    {
      char text[300];
      sprintf( text, "\n (P_FGMRESD::Iterate)    %d. %s %10.4le\n",
                     iters, "Iteration | accuracy =", relacc );
      REPORT::rpt.Message( 2, text );

      convergent = true;
      out_flag   = false;
    }
    else if( iters >= maxIter )
    {
      dcopy( neq, xs, x );

      char text[300];
      sprintf( text, "\n (P_FGMRESD::Iterate)    %d. %s %10.4le\n",
                     itac, "Iteration | accuracy =", accuracy );
      REPORT::rpt.Message( 2, text );

      convergent = false;
      out_flag   = false;
    }
  }

  // -------------------------------------------------------------------------------------
  // PART 3:  Release the working buffer
  // -------------------------------------------------------------------------------------

  for( int i=0; i<im+1; i++ )
  {
    delete[] vv[i];
  }
  delete[] vv;

  for( int i=0; i<im; i++ )
  {
    delete[] w[i];
  }
  delete w;

  for( int i=0; i<im; i++ )
  {
    delete[] hh[i];
  }
  delete[] hh;

  MEMORY::memo.Detach( c );
  MEMORY::memo.Detach( s );
  MEMORY::memo.Detach( rs );
  MEMORY::memo.Detach( hhloc );
  MEMORY::memo.Detach( xs );

  iterCountCG = iters;

  // exchange integer value convergent (secure)
# ifdef _MPI_
  convergent = project->subdom.Mpi_max( convergent );
# endif

  return convergent;
}
