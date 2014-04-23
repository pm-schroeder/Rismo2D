// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class P_BCGSTABD
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

//////////////////////////////////////////////////////////////////////////////////////////
// P_BCGSTABD: Bi Conjugate Gradient stabilized (BCGSTAB)
//             with right-hand preconditioner M^{-1}.
//             This is an improved BCG routine:
//             (1) no matrix transpose is involved;
//             (2) the convergence is smoother.
// ---------------------------------------------------------------------------------------
//   AUTHOR: Yu Liang, Dept. of CS, University of Minnesota
//   Date:  9/5/2001.
//
// ---------------------------------------------------------------------------------------
//   code changes
//
//   date        author
//   ----------  --------------
//   2004-08-14  P.M. Schroeder  adaption for use in Finite ELement Program Rismo2D
//
//                               changed variable names ----------------------------------
//                               nloc -> neq                    | number of equations
//                               sol  -> x                      | solution
//                               rhs  -> b                      | right hand side
//
//                               removed variables ---------------------------------------
//                               comm, wk1, wk2
//
//                               changed functions ---------------------------------------
//                               amxdis                 -> crsm->MulVec
//                               precon->prec_operation -> precon->Solve
//
//                               BLAS-functions ------------------------------------------
//                               DSCAL(n,a,x,1)         -> dscal(n,a,x)
//                                                         x[1,n} *= a
//                               DAXPY(n,a,x,1,y,1)     -> daxpy(n,a,x,y)
//                                                         y[1,n] += a * x[1,n]
//
//                               removed functions ---------------------------------------
//                               PARMS_malloc, CreateVec, VecAssign, CopyComm
//
// ///////////////////////////////////////////////////////////////////////////////////////

#include "Defs.h"
#include "Report.h"
#include "Project.h"
#include "Memory.h"
#include "CRSMat.h"
#include "Eqs.h"
#include "Precon.h"

#include "P_bcgstabd.h"

//#define kCGWarning


P_BCGSTABD::P_BCGSTABD()
{
  solverType = kParmsBcgstabd;
};


P_BCGSTABD::~P_BCGSTABD()
{
};


int P_BCGSTABD::Iterate( PROJECT* project, CRSMAT* crsm, double* b, double* x,
                         PRECON* precon )
{
  // -------------------------------------------------------------------------------------
  //  preconditioned BCGSTAB  -- Distributed version
  // -------------------------------------------------------------------------------------
  //
  //  original definition:
  //  void bcgstabd( DistMatrix dm, PreCon precon, IterPar ipar, Vec rhs, Vec x )
  //
  //  ON ENTRY
  // ==========
  //     crsm  =  (distributed) matrix in compact row storage format
  //   precon  =  a pointer to precon structure
  //      b    =  right hand side vector
  //
  //  ON RETURN
  // ===========
  //        x  =  local solution vector handler
  // -------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------
  //  Algorithm: ("Iterative Methods for Sparse Linear Systems",Saad96)
  //      ------------------------------------------------------------
  //      Initialization: r = b - AM^{-1}x, r0=r, p=r, rho=(r0, r),
  //      Iterate:
  //      (1) v = AM^{-1} p
  //      (2) alpha = rho / (r0, v)
  //      (3) s = r - alpha v
  //      (4) t = AM^{-1} s
  //      (5) omega = (t, s) / (t, t)
  //      (6) x = x + alpha * p + omega * s
  //      (7) r = s - omega * t
  //      convergence test goes here
  //      (8) beta = rho, rho = (r0, r), beta = rho*alpha/(beta*omega)
  //          p = r + beta * (p - omega * v)
  // -------------------------------------------------------------------------------------

  double alpha, beta, tscal;
  double rho, rho_r0v;
  double omega, omega_ts, omega_tt;
  double relacc;

  // -------------------------------------------------------------------------------------
  //                    PART 1: Initialization
  // -------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------
  //  section 1.1: retrieve information from input data structure
  // -------------------------------------------------------------------------------------

  int neq    = crsm->m_neq;             // dimension of subdomain
  int neq_dn = crsm->m_neq_dn;

  double eps = this->maxDiff;           // precision of computing

  // -------------------------------------------------------------------------------------
  //  section 1.3: allocate working buffer
  // -------------------------------------------------------------------------------------

  double* r0  = (double*) MEMORY::memo.Array_eq( neq );     // initial residual
  double* r   = (double*) MEMORY::memo.Array_eq( neq );     // current residual
  double* s   = (double*) MEMORY::memo.Array_eq( neq );
  double* t   = (double*) MEMORY::memo.Array_eq( neq );
  double* v   = (double*) MEMORY::memo.Array_eq( neq );
  double* p   = (double*) MEMORY::memo.Array_eq( neq );
  double* tmp = (double*) MEMORY::memo.Array_eq( neq );
  double* xs  = (double*) MEMORY::memo.Array_eq( neq );     // best solution up to iters


  // -------------------------------------------------------------------------------------
  //  section 1.4: compute the initial residual: r_0=b-AM^{-1}x_0
  // -------------------------------------------------------------------------------------

  // preconditioning operation  ----------------------------------------------------------

  if( precon )
  {
    precon->Solve( project, eqs, x, tmp );                  // tmp = M^{-1}*x
    crsm->MulVec( tmp, r, project, eqs );                   // r   = A*tmp
  }
  else
  {
    crsm->MulVec( x, r, project, eqs );                     // r = A*x
  }

  tscal = -1.0;
  dscal( neq, tscal, r );

  tscal = -tscal;
  daxpy( neq, tscal, b, r );                                // r = b - AM^{-1}*x

  dcopy( neq, r, r0 );                                      // r0 = r
  dcopy( neq, r, p );                                       // p  = r

  rho = ddot( neq_dn, r0, r );                              // rho = (r0,r)
# ifdef _MPI_
  rho = project->subdom.Mpi_sum( rho );
# endif

  double resid = sqrt( rho );

  if( resid > 1.001 )
  {
    for( int i=0; i<neq; i++ )  x[i] = 0.0;

    dcopy( neq, b, r );                                     // r = b

    dcopy( neq, r, r0 );                                    // r0 = r
    dcopy( neq, r, p );                                     // p  = r

    rho = ddot( neq_dn, r0, r );                            // rho = (r0,r)
#   ifdef _MPI_
    rho = project->subdom.Mpi_sum( rho );
#   endif

    resid = sqrt( rho );
  }

  // output the original residual norm ---------------------------------------------------

  REPORT::rpt.Message( 2, "\n (P_BCGSTABD::Iterate)   %s = %10.4le\n",
                          "original L2-Norm ||r||", resid );


  // -------------------------------------------------------------------------------------
  //                    PART 2: Iteration
  // -------------------------------------------------------------------------------------

  this->accuracy = 1.0;

  int convergent = false;
  int iters      = 0;         // counter for iterations
  int itac       = 0;         // iters where the best solution occured

  while( !convergent )
  {
    iters++;

    // -----------------------------------------------------------------------------------
    //  section 2.1:    v = A*M^{-1}*p
    // -----------------------------------------------------------------------------------

    if( precon )
    {
      precon->Solve( project, eqs, p, tmp );                // tmp = M^{-1}*p
      crsm->MulVec( tmp, v, project, eqs );                 // v   = A*tmp
    }
    else
    {
      crsm->MulVec( p, v, project, eqs );                   // v = A*p
    }

    // -----------------------------------------------------------------------------------
    //  section 2.2: alpha=rho/(r0, v)
    // -----------------------------------------------------------------------------------

    rho_r0v = ddot( neq_dn, r0, v );
#   ifdef _MPI_
    rho_r0v = project->subdom.Mpi_sum( rho_r0v );
#   endif

#   ifdef kCGWarning
    if( fabs(rho_r0v) < kZero )
    {
      REPORT::rpt.Screen( "\n WARNING: (r0, v) is very small (%le)\n", rho_r0v );
    }
#   endif

    alpha = rho / rho_r0v;

    // -----------------------------------------------------------------------------------
    // section 2.3: s = r - alpha*v
    // -----------------------------------------------------------------------------------

    dcopy( neq, r, s );
    tscal = -alpha;
    daxpy( neq, tscal, v, s );

    // -----------------------------------------------------------------------------------
    // section 2.4: t = A*M^{-1} * s
    // -----------------------------------------------------------------------------------

    if( precon )
    {
      precon->Solve( project, eqs, s, tmp );                // tmp = M^{-1}*s
      crsm->MulVec( tmp, t, project, eqs );                 // t = A*tmp
    }
    else
    {
      crsm->MulVec( s, t, project, eqs );                   // t = A*s
    }

    // -----------------------------------------------------------------------------------
    // section 2.5: omega = (t,s)/(t,t)
    //    test breakdown according to omega_j
    // -----------------------------------------------------------------------------------

    omega_ts = ddot( neq_dn, t, s );
#   ifdef _MPI_
    omega_ts = project->subdom.Mpi_sum( omega_ts );
#   endif

    omega_tt = ddot( neq_dn, t, t );
#   ifdef _MPI_
    omega_tt = project->subdom.Mpi_sum( omega_tt );
#   endif

#   ifdef kCGWarning
    if( fabs(omega_tt) < kZero )
    {
      REPORT::rpt.Screen( "\n WARNING: (t,t) is very small (%le)\n", omega_tt );
    }
#   endif

    omega = omega_ts / omega_tt;

    // -----------------------------------------------------------------------------------
    // section 2.6: x = x  +  alpha * p  +  omega * s
    // -----------------------------------------------------------------------------------

    daxpy( neq, alpha, p, x );
    daxpy( neq, omega, s, x );

    // -----------------------------------------------------------------------------------
    //   section 2.7:  r = s - omega * t
    // -----------------------------------------------------------------------------------

    dcopy( neq, s, r );
    tscal = -omega;
    daxpy( neq, tscal, t, r );

    // -----------------------------------------------------------------------------------
    //  section 2.8: p = beta * p  +  r  -  beta * omega * v
    // -----------------------------------------------------------------------------------

    beta = rho;

    rho = ddot( neq_dn, r0, r );
#   ifdef _MPI_
    rho = project->subdom.Mpi_sum( rho );
#   endif

    beta = rho * alpha / (beta * omega);

    dscal( neq, beta, p );

    tscal = 1.0;
    daxpy( neq, tscal, r, p);

    tscal = -beta * omega;
    daxpy( neq, tscal, v, p );

    // -----------------------------------------------------------------------------------
    //  section 2.9:  Convergence Test
    // -----------------------------------------------------------------------------------

    relacc = ddot( neq_dn, r, r );
#   ifdef _MPI_
    relacc = project->subdom.Mpi_sum( relacc );
#   endif

    relacc = sqrt( relacc ) / resid;

    // output residual norm --------------------------------------------------------------

#   ifdef kIteratCount
    REPORT::rpt.Screen( "\r (P_BCGSTABD::Iterate)   %d. %s = %10.4le",
            iters, "Iteration | accuracy", relacc );
    REPORT::rpt.PrintTime( -1 );
#   endif

//  if( relacc < accuracy )
    {
      itac = iters;
      accuracy = relacc;
      dcopy( neq, x, xs );
    }

    if( relacc < eps )
    {
      REPORT::rpt.Message( 2, "\n (P_BCGSTABD::Iterate)   %d. %s %10.4le\n",
                              iters, "Iteration | accuracy =", relacc );
      convergent = true;
    }

    else if( iters >= maxIter )
    {
      dcopy( neq, xs, x );

      REPORT::rpt.Message( 2, "\n (P_BCGSTABD::Iterate)   %d. %s %10.4le\n",
                              itac, "Iteration | accuracy =", accuracy );
      break;
    }
  } // END OF while


  // -------------------------------------------------------------------------------------
  //    PART 3: Obtain the very solution of original system
  // -------------------------------------------------------------------------------------

  if( precon )
  {
    dcopy( neq, x, tmp );
    precon->Solve( project, eqs, tmp, x );
  }

  // -------------------------------------------------------------------------------------
  //    PART 4: finish: detach temporary used memory
  // -------------------------------------------------------------------------------------

  MEMORY::memo.Detach( r0 );
  MEMORY::memo.Detach( r );
  MEMORY::memo.Detach( s );
  MEMORY::memo.Detach( t );
  MEMORY::memo.Detach( v );
  MEMORY::memo.Detach( p );
  MEMORY::memo.Detach( tmp );
  MEMORY::memo.Detach( xs );

  iterCountCG = iters;

  // exchange integer value convergent (secure)
# ifdef _MPI_
  convergent = project->subdom.Mpi_max( convergent );
# endif

  return convergent;
}
