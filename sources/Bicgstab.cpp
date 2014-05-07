// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software.
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================

//////////////////////////////////////////////////////////////////////////////////////////
// conjugate gradient solver Bi-CGSTAB (van der VORST, 1992)
//////////////////////////////////////////////////////////////////////////////////////////

#include "Defs.h"
#include "Report.h"
#include "Project.h"
#include "CRSMat.h"
#include "Eqs.h"
#include "Memory.h"
#include "Precon.h"
#include "Bicgstab.h"

//#define kCGWarning
//#define kDebug

#define kBiCGStabRestart  1.0e-30


BICGSTAB::BICGSTAB()
{
  solverType = kBicgstab;
}


BICGSTAB::~BICGSTAB()
{
}


//////////////////////////////////////////////////////////////////////////////////////////

int BICGSTAB::Iterate( PROJECT* project,
                       CRSMAT*  crsmat,
                       double*  b,
                       double*  x,
                       PRECON*  precon )
{
  // allocate memory for vectors ---------------------------------------------------------

  int neq    = crsmat->m_neq;
  int neq_dn = crsmat->m_neq_dn;

  double* r0 = (double*) MEMORY::memo.Array_eq( neq );
  double* ri = (double*) MEMORY::memo.Array_eq( neq );
  double* vi = (double*) MEMORY::memo.Array_eq( neq );
  double* pi = (double*) MEMORY::memo.Array_eq( neq );
  double* yi = (double*) MEMORY::memo.Array_eq( neq );
  double* zi = (double*) MEMORY::memo.Array_eq( neq );
  double* ti = (double*) MEMORY::memo.Array_eq( neq );
  double* t_ = (double*) MEMORY::memo.Array_eq( neq );
  double* xm = (double*) MEMORY::memo.Array_eq( neq );      // best solution for X

  memcpy( xm, x, neq*sizeof(double) );


  // initialisation ----------------------------------------------------------------------

  double alfa  = 1.0;
  double rho   = 1.0;
  double omega = 1.0;
  double beta;

  for( int j=0; j<neq; j++ )
  {
    vi[j] = 0.0;
    pi[j] = 0.0;
  }


  // L2-Norm of right side (to check for convergence) ------------------------------------

//  double l2r0 = 0.0;
//
//  for( int j=0; j<neq_dn; j++ )  l2r0 += b[j] * b[j];
//
//# ifdef _MPI_
//  l2r0 = project->subdom.Mpi_sum( l2r0 );
//# endif
//
//  l2r0 = sqrt( l2r0 );
//
//  REPORT::rpt.Message( "\n (BICGSTAB::Iterate)     %s = %10.4le\n",
//                       "L2-norm of RHS ||B||", l2r0 );
//
//  accuracy = 1.0;


  // -------------------------------------------------------------------------------------

  crsmat->MulVec( x, t_, project, eqs );

  for( int j=0; j<neq; j++ )  r0[j] = ri[j] = b[j] - t_[j];


  // L2-Norm -----------------------------------------------------------------------------

  double l2r0 = 0.0;

  for( int j=0; j<neq_dn; j++ ) l2r0 += r0[j] * r0[j];

# ifdef _MPI_
  l2r0 = project->subdom.Mpi_sum( l2r0 );
# endif

  l2r0 = sqrt( l2r0 );

  REPORT::rpt.Message( 3, "\n (BICGSTAB::Iterate)     %s = %10.4le\n",
                          "original L2-Norm ||r||", l2r0 );

  accuracy = 1.0;


  // -------------------------------------------------------------------------------------
  // start of iterations

  double diff;
  int    itac = 0;
  int    rest = 0;

  int it = 0;

  for( ;; )
  {
    it++;

    double temp1, temp2;

    temp1 = 0.0;
    for( int j=0; j<neq_dn; j++ )  temp1 += r0[j] * ri[j];

#   ifdef _MPI_
    temp1 = project->subdom.Mpi_sum( temp1 );
#   endif


//**REPORT::rpt.Message( " *** temp1 = %12.4le\n", temp1 );

    if( fabs(temp1) < kBiCGStabRestart )
    {
      rest++;

      crsmat->MulVec( x, ri, project, eqs );

      for( int j=0; j<neq; j++ )
      {
        t_[j] = b[j] - ri[j];
        ri[j] = t_[j];
        r0[j] = t_[j];
        vi[j] = 0.0;
        pi[j] = 0.0;
      }

      alfa  = 1.0;
      rho   = 1.0;
      omega = 1.0;

      temp1 = 0.0;
      for( int j=0; j<neq_dn; j++ )  temp1 += r0[j] * ri[j];

#     ifdef _MPI_
      temp1 = project->subdom.Mpi_sum( temp1 );
#     endif
    }


    beta = temp1 / rho  *  alfa / omega;
    rho  = temp1;


    // -----------------------------------------------------------------------------------

    for( int j=0; j<neq; j++ )  pi[j] = ri[j]  +  beta * (pi[j] - omega * vi[j]);


    // -----------------------------------------------------------------------------------

    if( precon )
    {
      precon->Solve( project, eqs, pi, yi );
      crsmat->MulVec( yi, vi, project, eqs );
    }
    else
    {
      crsmat->MulVec( pi, vi, project, eqs );
    }

//**REPORT::rpt.Message( 5, " *** ILU_solver finished\n" );

    // -----------------------------------------------------------------------------------

    temp1 = 0.0;
    for( int j=0; j<neq_dn; j++ )  temp1 += r0[j] * vi[j];

#   ifdef _MPI_
    temp1 = project->subdom.Mpi_sum( temp1 );
#   endif

//**REPORT::rpt.Message( 5, " *** temp1 = %12.4le\n", temp1 );

#   ifdef kCGWarning
    if( fabs(temp1) < kZero )
    {
      REPORT::rpt.Screen( "\n WARNING: (r0,vi) is very small (%le)\n", temp1 );
    }
#   endif

    alfa = rho / temp1;


    // -----------------------------------------------------------------------------------

    for( int j=0; j<neq; j++ )  ri[j] -= alfa * vi[j];


    // check for convergence -------------------------------------------------------------

    diff = 0.0;

    for( int j=0; j<neq_dn; j++ )  diff += ri[j] * ri[j];

#   ifdef _MPI_
    diff = project->subdom.Mpi_sum( diff );
#   endif

    diff = sqrt( diff ) / l2r0;

//**REPORT::rpt.Message( 5, " *** diff = %12.4le\n", diff );

#   ifdef kIteratCount
    REPORT::rpt.Screen( "\r (BICGSTAB::Iterate)     %d. %s %s %10.4le",
            it+1, "Iteration ",
            "| accuracy =", diff );
#   endif

    if( diff < accuracy )
    {
      itac = it;
      accuracy = diff;

      if( precon )  for( int j=0; j<neq; j++ )  xm[j] = x[j] + alfa * yi[j];
      else          for( int j=0; j<neq; j++ )  xm[j] = x[j] + alfa * pi[j];
    }

    if( diff < maxDiff )
    {
      REPORT::rpt.Message( 3, "\n (BICGSTAB::Iterate)     %d. %s %s %10.4le %s %d\n",
                              it+1, "Iteration ",
                              "| accuracy =", diff,
                              "| restarts =", rest );

      if( precon )  for( int j=0; j<neq; j++ )  x[j] += alfa * yi[j];
      else          for( int j=0; j<neq; j++ )  x[j] += alfa * pi[j];
      break;
    }

    else if( it >= maxIter )
    {
      memcpy( x, xm, neq*sizeof(double) );

      REPORT::rpt.Message( 3, "\n (BICGSTAB::Iterate)     %d. %s %s %10.4le %s %d\n",
                              itac+1, "Iteration ",
                              "| accuracy =", accuracy,
                              "| restarts =", rest );
      break;
    }


    // -----------------------------------------------------------------------------------

    if( precon )
    {
      precon->Solve( project, eqs, ri, zi );
      crsmat->MulVec( zi, ti, project, eqs );
    }
    else
    {
      crsmat->MulVec( ri, ti, project, eqs );
    }

    // -----------------------------------------------------------------------------------

    temp1 = 0.0;
    temp2 = 0.0;

    for( int j=0; j<neq_dn; j++ )  temp1 += ti[j] * ri[j];
    for( int j=0; j<neq_dn; j++ )  temp2 += ti[j] * ti[j];

#   ifdef _MPI_
    temp1 = project->subdom.Mpi_sum( temp1 );
    temp2 = project->subdom.Mpi_sum( temp2 );
#   endif

#   ifdef kCGWarning
    if( fabs(temp2) < kZero )
    {
      REPORT::rpt.Screen( "\n WARNING: (ti,ti) is very small (%le)\n", temp2 );
    }
#   endif

    omega = temp1 / temp2;


    // -----------------------------------------------------------------------------------

    if( precon )  for( int j=0; j<neq; j++ )  x[j] += alfa * yi[j]  +  omega * zi[j];
    else          for( int j=0; j<neq; j++ )  x[j] += alfa * pi[j]  +  omega * ri[j];

    // -----------------------------------------------------------------------------------

    for( int j=0; j<neq; j++ )  ri[j] -= omega * ti[j];
  }


  // finish: detach temporary used memory ------------------------------------------------

  MEMORY::memo.Detach( r0 );
  MEMORY::memo.Detach( ri );
  MEMORY::memo.Detach( vi );
  MEMORY::memo.Detach( pi );
  MEMORY::memo.Detach( yi );
  MEMORY::memo.Detach( zi );
  MEMORY::memo.Detach( ti );
  MEMORY::memo.Detach( t_ );
  MEMORY::memo.Detach( xm );

  iterCountCG = it;

  int convergent;

  if( it < maxIter  &&  diff < maxDiff )  convergent = true;
  else                                    convergent = false;

  // exchange integer value convergent (secure)
# ifdef _MPI_
  convergent = project->subdom.Mpi_max( convergent );
# endif

  return convergent;
}
