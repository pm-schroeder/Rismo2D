// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class FRONTM
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
#include "Memory.h"
#include "Elem.h"
#include "Node.h"
#include "Project.h"
#include "CRSMat.h"
#include "Fromat.h"
#include "Solver.h"
#include "P_bcgstabd.h"

#include "Frontm.h"

//#define kEqCounter

#define kMinPivot  1.0e-20


FRONTM::FRONTM()
{
  id = NULL;
  solverType = kFrontm;
}


FRONTM::~FRONTM()
{
}


//////////////////////////////////////////////////////////////////////////////////////////
// MPI:
// LHS matrix "crsm" and  RHS vector "b" are expected to be in local storage
// global assembling of "vec" is performed after elimination of internal nodes

void FRONTM::Direct( PROJECT* project, CRSMAT* crsm, double* b, double* vec )
{
  REPORT::rpt.PrintTime( 1 );
  REPORT::rpt.Message( 2, "\n" );


  ////////////////////////////////////////////////////////////////////////////////////////
  // 1.  Initializations

  //REPORT::rpt.Message( 5, " (FRONTM::Direct)        1.  Initializations\n" );

  char     text[250];

  int      neq    = crsm->m_neq;                  // total number of equations
  int      neq_up = crsm->m_neq_up;               // start number of upstream equations
  int      neq_dn = crsm->m_neq_dn;               // start number of downstream equations

  int*     width  = crsm->m_width;
  int**    index  = crsm->m_index;
  REALPR** A      = crsm->m_A;

  int  maximumFW  = 0;

  memcpy( vec, b, neq*sizeof(double) );           // copy vector "b" to "vec"


  // -------------------------------------------------------------------------------------
  // allocate dynamic front memory

  // matrix to hold lower and upper part of the factorized equation system
  CRSMAT* solve = new CRSMAT( neq );


  // -------------------------------------------------------------------------------------

  FROMAT* fromat = new FROMAT( mfw, neq );

  if( !fromat )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory - FRONTM::Direct(3)" );

  // -------------------------------------------------------------------------------------

  EQL* eql = NULL;
  if( neq > 0 )  eql = &fromat->m_eql[0];


  ////////////////////////////////////////////////////////////////////////////////////////
  // 2.  forward factorization: build upper diagonal of equation matrix
  //
  // MPI: restricted to interior subdomain nodes: e < neq_up

  //REPORT::rpt.Message( 5, " (FRONTM::Direct)        2.  forward factorization\n" );

  int eqCounter = 0;

  while( eql )
  {
    int e = eql->no;

    if( e >= neq_up )  break;

    // -----------------------------------------------------------------------------------
    // insert equation "e"
    fromat->Insert( e, width, index, A );

    if( fromat->m_actFW > maximumFW )  maximumFW = fromat->m_actFW;

    eqCounter++;

#   ifdef kEqCounter
    if( eqCounter%100 == 0 )
    {
      REPORT::rpt.Screen( "\r (FRONTM::Direct)        %d (maximum FW = %d)",
                          eqCounter, maximumFW );
    }
#   endif

    // -----------------------------------------------------------------------------------
    // ... and eliminate (only equations with eqno < neq_up)

    eql = Eliminate( fromat, eql, neq_up, solve, width, index, A, vec );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 3.  insert equations from up- and downstream interface without factorization

  //REPORT::rpt.Message( 5, " (FRONTM::Direct)        3.  interface equations\n" );

  for( int e=neq_up; e<neq; e++ )
  {
    // insert equation "e" ---------------------------------------------------------------
    fromat->Insert( e, width, index, A );

    if( fromat->m_actFW > maximumFW )  maximumFW = fromat->m_actFW;
  }

  // check for singularity, other errors and give size message ---------------------------
  REPORT::rpt.Screen( 2, "\r (FRONTM::Direct)        %5d (maximum FW = %5d)",
                         neq_up+1, maximumFW );

  sprintf( text, " (FRONTM::Direct)        %d was maximum front width\n", maximumFW );
  REPORT::rpt.Output( text, 2 );

  for( int i=0; i<fromat->m_actFW; i++ )
  {
    if( fromat->m_frow[i].no < (unsigned int) neq_up )
      REPORT::rpt.Error( "singular matrix - FRONTM::Direct(6)" );
  }

  if( eqCounter != neq_up )
  {
    REPORT::rpt.Error( "bad number of equations - FRONTM::Direct(7)" );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // 4.  (MPI) solve iteratively for interface nodes
  //
# ifdef _MPI_

  //REPORT::rpt.Message( " (FRONTM::Direct)        4.  MPI\n" );

  // assemble local vector "vec" from all adjacent subdomains ----------------------------

  eqs->Mpi_assemble( vec, project );


  // solve the remaining matrix (interface nodes) with iterative solver BICGSTAB ---------

  if( project->subdom.npr > 1 )
  {
    // check number of equations
    int n = fromat->m_actFW;

    if( n != neq - neq_up )
      REPORT::rpt.Error( "singular matrix - FRONTM::Direct(4)" );

    double* frB = (double*) MEMORY::memo.Array_eq( n );
    double* frX = (double*) MEMORY::memo.Array_eq( n );

    if( !frB || !frX )
      REPORT::rpt.Error( kMemoryFault, "can not allocate memory - FRONTM::Direct(5)" );

    for( int i=0; i<n; i++ )
    {
      int eqno = fromat->m_frow[i].no;

      frB[i] = vec[eqno];
      frX[i] = 0.0;
    }


    // determine the L2-Norm of vec = scaling factor -------------------------------------

    double scale = 0.0;
    for( int i=0; i<n; i++ )  scale += frB[i] * frB[i];

    scale = project->subdom.Mpi_sum( scale );

    REPORT::rpt.Message( 2, "\n (FRONTM::Direct)        %s = %10.4le\n",
                        "L2-Norm of RHS ||b||", sqrt(scale) );


    // diagonal scaling ------------------------------------------------------------------
    if( fabs(scale) > 0.0 )
    {
      for( int i=0; i<n; i++ )
      {
        frB[i] /= scale;
        for( int j=0; j<n; j++ )  fromat->m_frow[i].eq[j] /= scale;
      }
    }

    // solve with conjugate gradient solver
    BiCGStab( project, fromat, frB, frX );

    // copy solution to vector
    for( int i=0; i<n; i++ )
    {
      int eqno = fromat->m_frow[i].no;
      vec[eqno] = frX[i];
    }

    MEMORY::memo.Detach( frB );
    MEMORY::memo.Detach( frX );
  }

# endif
  ////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////
  // 5.  backward solve for interior nodes

  //REPORT::rpt.Message( 5, " (FRONTM::Direct)        5.  backward solve\n" );

  eql = NULL;

  if( neq > 0 )  eql = fromat->m_eql[0].Last();

  while( eql )
  {
    int e = eql->no;

    if( e < neq_up )
    {
      int     w = solve->m_width[e];
      int*    I = solve->m_index[e];
      REALPR* A = solve->m_A[e];

      for( int j=1; j<w; j++ )
      {
        vec[e] -= A[j] * vec[ I[j] ];
      }

      vec[e] /= A[0];
    }

    eql = eql->prev;
  }


  ////////////////////////////////////////////////////////////////////////////////////////

  REPORT::rpt.Screen( 2, "\n" );


  // free temporary allocated memory -----------------------------------------------------

  delete solve;
  delete fromat;
}


//////////////////////////////////////////////////////////////////////////////////////////
// Eliminate equation eql or next equation from format and
// append elimination to matrix solve.
//////////////////////////////////////////////////////////////////////////////////////////

EQL* FRONTM::Eliminate( FROMAT*  fromat,
                        EQL*     eql,
                        int      neq,
                        CRSMAT*  solve,
                        int*     width,
                        int**    index,
                        REALPR** A,
                        double*  vec )
{
  int eqno  = eql->no;
  int actFW = fromat->m_actFW;

  // -------------------------------------------------------------------------------------
  // search for pivot
  int    ixe   = fromat->m_eql[eqno].ind;
  int    ixpiv = ixe;
  double pivot = fromat->m_frow[ixpiv].eq[ixe];

  for( int r=0; r<actFW; r++ )
  {
    double diag = fromat->m_frow[r].eq[r];
    double ecol = fromat->m_frow[r].eq[ixe];

    if( fromat->m_frow[r].no < (unsigned int) neq )
    {
      if( fabs(ecol) > 10.0 * fabs(pivot)  &&  fabs(diag) > 10.0 * fabs(pivot) )
      {
        pivot = diag;
        ixpiv = r;
      }
    }
  }

  int eleq = fromat->m_frow[ixpiv].no;

  // -------------------------------------------------------------------------------------
  // proceed to next equation from equation list eql
  if( eleq == eqno )
  {
    eql = eql->next;
  }
  else
  {
    // insert missing columns of elimination equation "eleq" (if necessary) --------------
    fromat->Insert( eleq, width, index, A );

    // change equation order in EQ(uation)L(ist)
    fromat->m_eql[eleq].Remove();
    fromat->m_eql[eleq].Insert( eql );
  }

  // eliminate elimination equation ------------------------------------------------------
  fromat->Eliminate( eleq, vec );

  // -------------------------------------------------------------------------------------
  // save elimination equation "e" to CRS-matrix
  solve->Append( eleq, pivot, fromat );

  // -------------------------------------------------------------------------------------

  return eql;
}


#define kBiCGStabRestart  1.0e-30
//#define kIteratCount

int FRONTM::BiCGStab( PROJECT* project, FROMAT* fromat, double* b, double* x )
{
  maxIter = 5000;
  maxDiff = 1.0e-4;

  // allocate memory for vectors ---------------------------------------------------------

  int neq    = fromat->m_actFW;

  double* r0 = (double*) MEMORY::memo.Array_eq( neq );
  double* ri = (double*) MEMORY::memo.Array_eq( neq );
  double* vi = (double*) MEMORY::memo.Array_eq( neq );
  double* pi = (double*) MEMORY::memo.Array_eq( neq );
  double* ti = (double*) MEMORY::memo.Array_eq( neq );
  double* t_ = (double*) MEMORY::memo.Array_eq( neq );
  double* xm = (double*) MEMORY::memo.Array_eq( neq );      // best solution for x

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
//  for( int j=0; j<neq; j++ )
  //{
  //  if( fromat->m_frow[j].no < eqs->neq_dn )  l2r0 += b[j] * b[j];
  //}
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

  fromat->MulVec( x, t_, eqs, &project->subdom );

  for( int j=0; j<neq; j++ )  r0[j] = ri[j] = b[j] - t_[j];


  // L2-Norm -----------------------------------------------------------------------------

  double l2r0 = 0.0;

  for( int j=0; j<neq; j++ )
  {
    if( fromat->m_frow[j].no < (unsigned)eqs->neq_dn )  l2r0 += r0[j] * r0[j];
  }
# ifdef _MPI_
  l2r0 = project->subdom.Mpi_sum( l2r0 );
# endif

  l2r0 = sqrt( l2r0 );

  REPORT::rpt.Message( 2, "\n (BICGSTAB::Iterate)     %s = %10.4le\n",
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
    for( int j=0; j<neq; j++ )
    {
      if( fromat->m_frow[j].no < (unsigned)eqs->neq_dn )  temp1 += r0[j] * ri[j];
    }
#   ifdef _MPI_
    temp1 = project->subdom.Mpi_sum( temp1 );
#   endif


    if( fabs(temp1) < kBiCGStabRestart )
    {
      rest++;

      fromat->MulVec( x, ri, eqs, &project->subdom );

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
      for( int j=0; j<neq; j++ )
      {
        if( fromat->m_frow[j].no < (unsigned)eqs->neq_dn )  temp1 += r0[j] * ri[j];
      }
#     ifdef _MPI_
      temp1 = project->subdom.Mpi_sum( temp1 );
#     endif
    }


    beta = temp1 / rho  *  alfa / omega;
    rho  = temp1;


    // -----------------------------------------------------------------------------------

    for( int j=0; j<neq; j++ )  pi[j] = ri[j]  +  beta * (pi[j] - omega * vi[j]);


    // -----------------------------------------------------------------------------------

    fromat->MulVec( pi, vi, eqs, &project->subdom );


    // -----------------------------------------------------------------------------------

    temp1 = 0.0;
    for( int j=0; j<neq; j++ )
    {
      if( fromat->m_frow[j].no < (unsigned)eqs->neq_dn )  temp1 += r0[j] * vi[j];
    }
#   ifdef _MPI_
    temp1 = project->subdom.Mpi_sum( temp1 );
#   endif

//**REPORT::rpt.Message( " *** temp1 = %12.4le\n", temp1 );

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

    for( int j=0; j<neq; j++ )
    {
      if( fromat->m_frow[j].no < (unsigned)eqs->neq_dn )  diff += ri[j] * ri[j];
    }
#   ifdef _MPI_
    diff = project->subdom.Mpi_sum( diff );
#   endif

    diff = sqrt( diff ) / l2r0;

//**REPORT::rpt.Message( " *** diff = %12.4le\n", diff );

#   ifdef kIteratCount
    REPORT::rpt.Screen( "\r (BICGSTAB::Iterate)     %d. %s %s %10.4le",
            it+1, "Iteration ",
            "| accuracy =", diff );
#   endif

    if( diff < accuracy )
    {
      itac = it;
      accuracy = diff;
      for( int j=0; j<neq; j++ )  xm[j] = x[j] + alfa * pi[j];
    }

    //for( int j=0; j<neq; j++ )
    //{
    //  char text[200];
    //  sprintf( text, "      xm[%04d] = %.4le\n", j, xm[j] );
    //  if( j % 50 == 0 )  REPORT::rpt.Output( text );
    //}

    if( diff < maxDiff )
    {
      REPORT::rpt.Message( 2, "\n (BICGSTAB::Iterate)     %d. %s %s %10.4le %s %d\n",
                              it+1, "Iteration ",
                              "| accuracy =", diff,
                              "| restarts =", rest );

      for( int j=0; j<neq; j++ )  x[j] += alfa * pi[j];
      break;
    }

    else if( it >= maxIter )
    {
      memcpy( x, xm, neq*sizeof(double) );

      REPORT::rpt.Message( 2, "\n (BICGSTAB::Iterate)     %d. %s %s %10.4le %s %d\n",
                              itac+1, "Iteration ",
                              "| accuracy =", accuracy,
                              "| restarts =", rest );
      break;
    }


    // -----------------------------------------------------------------------------------

    fromat->MulVec( ri, ti, eqs, &project->subdom );

    // -----------------------------------------------------------------------------------

    temp1 = 0.0;
    temp2 = 0.0;

    for( int j=0; j<neq; j++ )
    {
      if( fromat->m_frow[j].no < (unsigned)eqs->neq_dn )
      {
        temp1 += ti[j] * ri[j];
        temp2 += ti[j] * ti[j];
      }
    }
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

    for( int j=0; j<neq; j++ )  x[j] += alfa * pi[j]  +  omega * ri[j];


    // -----------------------------------------------------------------------------------

    for( int j=0; j<neq; j++ )  ri[j] -= omega * ti[j];
  }


  // finish: detach temporary used memory ------------------------------------------------

  MEMORY::memo.Detach( r0 );
  MEMORY::memo.Detach( ri );
  MEMORY::memo.Detach( vi );
  MEMORY::memo.Detach( pi );
  MEMORY::memo.Detach( ti );
  MEMORY::memo.Detach( t_ );
  MEMORY::memo.Detach( xm );

  iterCountCG = it;

  if( it < maxIter  &&  diff < maxDiff )  return true;
  else                                    return false;
}
