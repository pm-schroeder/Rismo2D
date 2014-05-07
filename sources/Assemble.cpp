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

#include "Defs.h"
#include "Report.h"
#include "Elem.h"
#include "Node.h"
#include "Model.h"
#include "Project.h"
#include "Subdom.h"

#include "CRSMat.h"


void CRSMAT::AssembleEstifm_im( EQS*     eqs,
                                MODEL*   model,
                                PROJECT* project )
{
  int i, j, k, l, m;

  REALPR*  APtr;
  double*  estifmPtr;

  int dfcn = eqs->dfcn;
  int dfel = eqs->dfel;

  double** estifm = eqs->estifm;

  REPORT::rpt.Message( 3, "\n (CRSMAT::Assemble...)   %s (%d elements)\n",
                          "assembling estifm", model->ne );


  // -------------------------------------------------------------------------------------
  // assemble elements

  for( int e=0; e<model->ne; e++ )                     // loop on elements
  {
    ELEM* el = model->elem[e];

    int row,  col;
    int rind, cind;


    // -----------------------------------------------------------------------------------
    // compute element coefficients

    if( !eqs->Coefs(el, project, estifm, (double *) 0) )  continue;

    int nnd = el->Getnnd();


    // -----------------------------------------------------------------------------------
    // insert element stiffness matrix (estifm)

    for( i=0; i<nnd; i++ )                  // loop on nodes
    {
      for( j=0; j<dfcn; j++ )               // loop on node-equations
      {
        rind = i + j*nnd;
        row  = eqs->GetEqno( el->nd[i], j );

        if( row >= 0 )
        {
          APtr      = m_A[row];
          estifmPtr = estifm[rind];

          for( k=0; k<nnd; k++ )            // loop on nodes
          {
            for( l=0; l<dfcn; l++ )         // loop on node-equations
            {
              cind = k + l*nnd;
              col  = eqs->GetEqno( el->nd[k], l );

              if( col >= 0 )
              {
                for( m=0; m<m_width[row]; m++ )
                {
                  if( col == m_index[row][m] )
                  {
                    APtr[m] += (REALPR) estifmPtr[cind];
                    break;
                  }
                }
              }
            }
          }

          for( l=0; l<dfel; l++ )           // loop on element-equations
          {
            cind = dfcn*nnd + l;
            col  = eqs->GetEqno( el, l );

            if( col >= 0 )
            {
              for( m=0; m<m_width[row]; m++ )
              {
                if( col == m_index[row][m] )
                {
                  APtr[m] += (REALPR) estifmPtr[cind];
                  break;
                }
              }
            }
          }
        }
      }
    }


    for( j=0; j<dfel; j++ )                 // loop on element-equations
    {
      rind = dfcn*nnd + j;
      row  = eqs->GetEqno( el, j );

      if( row >= 0 )
      {
        APtr      = m_A[row];
        estifmPtr = estifm[rind];

        for( k=0; k<nnd; k++ )              // loop on nodes
        {
          for( l=0; l<dfcn; l++ )           // loop on node-equations
          {
            cind = k + l*nnd;
            col  = eqs->GetEqno( el->nd[k], l );

            if( col >= 0 )
            {
              for( m=0; m<m_width[row]; m++ )
              {
                if( col == m_index[row][m] )
                {
                  APtr[m] += (REALPR) estifmPtr[cind];
                  break;
                }
              }
            }
          }
        }


        for( l=0; l<dfel; l++ )             // loop on element-equations
        {
           cind = dfcn*nnd + l;
           col  = eqs->GetEqno( el, l );

           if( col >= 0 )
           {
            for( m=0; m<m_width[row]; m++ )
            {
              if( col == m_index[row][m] )
              {
                APtr[m] += (REALPR) estifmPtr[cind];
                break;
              }
            }
          }
        }
      }
    }
  }
}


void CRSMAT::AssembleEqs_im( EQS*     eqs,
                             double*  vector,
                             MODEL*   model,
                             PROJECT* project )
{
  int neq  = eqs->neq;
  int dfcn = eqs->dfcn;
  int dfel = eqs->dfel;

  double*  force  = eqs->force;
  double** estifm = eqs->estifm;

  REPORT::rpt.Message( 3, "\n (CRSMAT::Assemble...)   %s (%d elements)\n",
                          "assembling eqs", model->ne );


  // initializations ---------------------------------------------------------------------

  for( int i=0; i<neq; i++ )  vector[i] = 0.0;


  // -------------------------------------------------------------------------------------
  // assemble elements

  for( int e=0; e<model->ne; e++ )
  {
    ELEM* el = model->elem[e];


    // compute element coefficients ------------------------------------------------------

    if( !eqs->Coefs(el, project, estifm, force) )  continue;

    int nnd = el->Getnnd();


    // insert element stiffness matrix (estifm) and (force) ------------------------------

    for( int i=0; i<nnd; i++ )               // loop on nodes
    {
      for( int j=0; j<dfcn; j++ )            // loop on node-equations
      {
        int rind = i + j*nnd;
        int row  = eqs->GetEqno( el->nd[i], j );

        if( row >= 0 )
        {
          REALPR* APtr      = m_A[row];
          double* estifmPtr = estifm[rind];

          vector[row] += force[rind];

          for( int k=0; k<nnd; k++ )         // loop on nodes
          {
            for( int l=0; l<dfcn; l++ )      // loop on node-equations
            {
              int cind = k + l*nnd;
              int col  = eqs->GetEqno( el->nd[k], l );

              if( col >= 0 )
              {
                for( int m=0; m<m_width[row]; m++ )
                {
                  if( col == m_index[row][m] )
                  {
                    APtr[m] += (REALPR) estifmPtr[cind];
                    break;
                  }
                }
              }
            }
          }

          for( int l=0; l<dfel; l++ )        // loop on element-equations
          {
            int cind = dfcn*nnd + l;
            int col  = eqs->GetEqno( el, l );

            if( col >= 0 )
            {
              for( int m=0; m<m_width[row]; m++ )
              {
                if( col == m_index[row][m] )
                {
                  APtr[m] += (REALPR) estifmPtr[cind];
                  break;
                }
              }
            }
          }
        }
      }
    }


    for( int j=0; j<dfel; j++ )              // loop on element-equations
    {
      int rind = dfcn*nnd + j;
      int row  = eqs->GetEqno( el, j );

      if( row >= 0 )
      {
        REALPR* APtr      = m_A[row];
        double* estifmPtr = estifm[rind];

        vector[row] += force[rind];

        for( int k=0; k<nnd; k++ )           // loop on nodes
        {
          for( int l=0; l<dfcn; l++ )        // loop on node-equations
          {
            int cind = k + l*nnd;
            int col  = eqs->GetEqno( el->nd[k], l );

            if( col >= 0 )
            {
              for( int m=0; m<m_width[row]; m++ )
              {
                if( col == m_index[row][m] )
                {
                  APtr[m] += (REALPR) estifmPtr[cind];
                  break;
                }
              }
            }
          }
        }


        for( int l=0; l<dfel; l++ )          // loop on element-equations
        {
          int cind = dfcn*nnd + l;
          int col  = eqs->GetEqno( el, l );

          if( col >= 0 )
          {
            for( int m=0; m<m_width[row]; m++ )
            {
              if( col == m_index[row][m] )
              {
                APtr[m] += (REALPR) estifmPtr[cind];
                break;
              }
            }
          }
        }
      }
    }
  }


  REPORT::rpt.Message( 3, "\n\n%-25s%s\n\n%15s %1s  %8s  %14s  %14s\n\n",
                          " (CRSMAT::Assemble...)", "Newton-Raphson-residuum / force vector ...",
                          " ", " ", "  node", "   average", "   maximum" );

  for( int e=0; e<dfcn; e++ )
  {
    int    no  = 0;
    double ave = 0.0;
    double max = 0.0;

    for( int n=0; n<model->np; n++ )
    {
      NODE* nd = model->node[n];
      int eqno = eqs->GetEqno( nd, e );

      if( eqno >= 0 )
      {
        double vec = vector[eqno];

        ave += vec;
        if( fabs(vec) > fabs(max) )
        {
          max = vec;
          no  = nd->Getname();
        }
      }
    }

    int tot = neq;

    //////////////////////////////////////////////////////////////////////////////////////
#   ifdef _MPI_
    max = project->subdom.Mpi_max( max );
    tot = project->subdom.Mpi_sum( eqs->neq_up );
    ave = project->subdom.Mpi_sum( ave );
#   endif
    //////////////////////////////////////////////////////////////////////////////////////

    if( tot )  ave /= tot;

    REPORT::rpt.Message( 3, " %15s %1d  %8d  %14.5le  %14.5le\n", " ", e+1, no, ave, max );
  }

  REPORT::rpt.Message( 3, "\n" );
}


void CRSMAT::AssembleForce( EQS*     eqs,
                            double*  vector,
                            MODEL*   model,
                            PROJECT* project )
{
  int i, j;

  int neq  = eqs->neq;
  int dfcn = eqs->dfcn;
  int dfel = eqs->dfel;

  double* force = eqs->force;

  REPORT::rpt.Message( 3, "\n (CRSMAT::Assemble...)   %s (%d elements)\n",
                          "assembling force", model->ne );


  // initializations ---------------------------------------------------------------------

  for( i=0; i<neq; i++ ) vector[i] = 0.0;


  // -------------------------------------------------------------------------------------
  // assemble elements

  for( int e=0; e<model->ne; e++ )
  {
    ELEM* el = model->elem[e];

    int row;
    int rind;

    int nnd = el->Getnnd();


    // compute element coefficients ------------------------------------------------------

    if( !eqs->Coefs(el, project, (double **) 0, force) )  continue;


    // insert force vector ---------------------------------------------------------------

    for( i=0; i<nnd; i++ )
    {
      for( j=0; j<dfcn; j++ )
      {
        rind = i + j*nnd;
        row  = eqs->GetEqno( el->nd[i], j );

        if( row >= 0 )
        {
          vector[row] += force[rind];
        }
      }
    }


    for( j=0; j<dfel; j++ )
    {
      rind = dfcn*nnd + j;
      row  = eqs->GetEqno( el, j );

      if( row >= 0 )
      {
        vector[row] += force[rind];
      }
    }
  }
}

