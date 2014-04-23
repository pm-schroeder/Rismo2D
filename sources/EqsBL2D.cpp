// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_BL2D
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
#include "Vars.h"
#include "Shape.h"
#include "Type.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Model.h"
#include "Project.h"

#include "EqsBL2D.h"


EQS_BL2D::EQS_BL2D() : EQS( 1, 1, 0 )
{
  neq  = 0;
}


EQS_BL2D::~EQS_BL2D()
{
}

// ---------------------------------------------------------------------------------------

void EQS_BL2D::Execute( PROJECT* project, int shape )
{
  MODEL*   model     = project->M2D;

  GRID*    rg        = model->region;
  int      rgnp      = rg->Getnp();
  int      rgne      = rg->Getne();

  SED*     sed       = &project->sed;
  TIMEINT* tmint     = &project->timeint;

  int      diverged_cg = 0;
  char     text[500];

  int      noAbs, noPer;
  double   maxAbs, maxPer, avAbs, avPer;

  model->Incinit();

  // print information on actual iteration -----------------------------------------------
  project->PrintTheCycle( 1 );
  REPORT::rpt.PrintTime( 1 );


  // use bedload transport capacity under equilibrium conditions -------------------------
  if( project->sed.lsType == 1 )
  {
    // set boundary conditions -----------------------------------------------------------
    sed->Bcon( project, model );

    // copy computed equilibrium transport capacity to NODE::v.Qb ------------------------
    // The (equilibrium) bed load transport capacity SED::qbc is computed during a call
    //  to SED::Initialize() in the cycle-switch in PROJECT::Compute().
    for( int n=0; n<rgnp; n++ )
    {
      NODE* nd = rg->Getnode(n);
      int   no = nd->Getno();

      nd->v.Qb = sed->qbc[no];
    }

    sprintf( text, "\n%-25s%s\n ",
                    " (EQS_BL2D::Execute)",
                    "equilibrium bedload transport rates determined" );
    REPORT::rpt.Output( text, 1 );
  }


  // ... use loading law to compute bedload transport capacity ---------------------------
  else
  {
    // compute friction coefficients -----------------------------------------------------
    model->DoFriction( project );

    // initialize Reynolds stresses and eddy viscosity -----------------------------------
    rg->Turbulence( project );

    // initialize turbulent diffusivity/dispersion ---------------------------------------
    rg->EddyDisp();

    // set up equation numbers -----------------------------------------------------------
    project->fix[0]   = 0;
    project->elemKind = ELEM::kRegion | ELEM::kBound;

    if( shape == kLinear )
    {
      SetEqno( model, 1, 0, 0, project->fix, project->elemKind );
      coefs = kLinear;
    }
    else
    {
      SetEqno( model, 1, 1, 0, project->fix, project->elemKind );
      coefs = kQuadratic;
    }

    initStructure = true;

    // allocate some temporary used arrays -----------------------------------------------
    //      B = right hand side of equation system
    //      X = solution of equation system

    double* B = (double*) MEMORY::memo.Array_nd( rgnp );
    double* X = (double*) MEMORY::memo.Array_nd( rgnp );

    // set boundary conditions -----------------------------------------------------------
    sed->Bcon( project, model );

    // initialize qb ---------------------------------------------------------------------
    for( int n=0; n<rgnp; n++ )
    {
      NODE* nd = rg->Getnode(n);
      int   no = nd->Getno();

      nd->v.Qb  = sed->qbc[no];
      nd->fixqb = false;
    }

    // -----------------------------------------------------------------------------------
    // solve loading law equation
    for( int n=0; n<rgnp; n++ )  X[n] = 0.0;

    diverged_cg = Solve( model, neq, B, X, project );

    // statistics of correction vector ---------------------------------------------------
    sprintf( text, "\n\n%-25s%s\n\n %s\n\n",
                  " (EQS_BL2D::Execute)", "convergence parameters ...",
                  " variable   node          average        maximum" );
    REPORT::rpt.Message( 1, text );

    Update( model, &project->subdom,
            X, 0, kVarQb, &maxAbs, &maxPer, &avAbs, &avPer, &noAbs, &noPer );

    sprintf( text, "      %2s    %7d     %12.5le   %12.5le %s\n",
                  "qb", noAbs, avAbs, maxAbs, "     (abs)" );
    REPORT::rpt.Message( 1, text );

    sprintf( text, "      %2s    %7d     %12.5lf   %12.5lf %s\n\n",
                  "  ", noPer, avPer, maxPer, "     ( %% )" );
    REPORT::rpt.Message( 1, text );

    // update node variables -------------------------------------------------------------
    for( int n=0; n<rgnp; n++ )
    {
      NODE* nd = rg->Getnode(n);
      int   no = nd->Getno();
      int eqno = GetEqno( nd, 0 );

      nd->v.Qb += X[eqno];
    }

    // -----------------------------------------------------------------------------------
    for( int n=0; n<rgnp; n++ )
    {
      NODE* nd = rg->Getnode(n);
      int   no = nd->Getno();

      // check for minimum allowed transport capacity ... --------------------------------
      if( nd->v.Qb < sed->minQb )
      {
        nd->v.Qb = 0.0;
      }
    }

    // interpolate transport rate qb at midside nodes ------------------------------------
    if( shape == kLinear )
    {
      for( int e=0; e<rg->Getne(); e++ )
      {
        ELEM* el = rg->Getelem(e);

        int ncn = el->Getncn();
        int nnd = el->Getnnd();

        for( int i=ncn; i<nnd; i++ )
        {
          int il, ir;
          el->GetQShape()->getCornerNodes( i, &il, &ir );

          NODE* ndm = el->Getnode(i);
          NODE* ndl = el->Getnode(il);
          NODE* ndr = el->Getnode(ir);

          ndm->v.Qb = 0.5 * (ndl->v.Qb + ndr->v.Qb);
        }
      }
    }

    // finalizing ------------------------------------------------------------------------
    MEMORY::memo.Detach( X );
    MEMORY::memo.Detach( B );

    sprintf( text, "\n%-25s%s\n ",
                    " (EQS_BL2D::Execute)",
                    "non-equilibrium bedload transport rates determined" );
    REPORT::rpt.Output( text, 1 );

    // copy computed non-equilibrium transport capacity to sed->qbc[] --------------------
    for( int n=0; n<rgnp; n++ )
    {
      NODE* nd = rg->Getnode(n);
      int   no = nd->Getno();

      sed->qbc[no] = nd->v.Qb;
    }
  }
}
