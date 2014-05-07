// ======================================================================================
//                                      G R I D
// ======================================================================================
// The method GRID::Dispersion() determines the dispersion coefficients
// due to streamline curvature.
// The method is replaced by class EQS_DISP and not longer used.
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
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 03.04.2006     sc     first implementation / first concept
//
// ======================================================================================

// ---------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Node.h"
#include "Elem.h"
#include "Memory.h"
#include "Model.h"
#include "Project.h"

#include "Grid.h"


void GRID::Dispersion( PROJECT* project,
                       double*  Duu,
                       double*  Dvv,
                       double*  Duv,
                       double*  Vsec )
{
  // determine curvature of stream lines -------------------------------------------------

  double* curv = project->M2D->Curv2D();


  // compute dispersion coefficients -----------------------------------------------------

  for( int n=0; n<np; n++ )
  {
    double U  = node[n].v.U;
    double V  = node[n].v.V;
    double H  = node[n].v.S - node[n].z;
    double Us = sqrt( U*U + V*V );

    double cf = node[n].cf;

    if( Us > 1.0e-6 )
    {
      double m   = project->kappa / sqrt( cf );
      double alf = (m + 0.5) / project->kappa / project->kappa / m;
      double acu = alf * curv[n];
      double HUs = H * Us;

      double Dxx = Us * HUs / m / m;
      double Dyy = acu * acu * HUs * HUs * H / 3.0;
      //double Dxy = acu * HUs * HUs / m / 2.0;
      double Dxy = acu * HUs * HUs / (2.0*m + 1.0);

      // rotation to mean flow direction
      double sina =  V / Us;
      double cosa =  U / Us;

      Duu[n] = cosa*cosa*Dxx - 2.0*sina*cosa*Dxy + sina*sina*Dyy;
      Dvv[n] = sina*sina*Dxx + 2.0*sina*cosa*Dxy + cosa*cosa*Dyy;
      Duv[n] = sina*cosa*(Dxx-Dyy) + (cosa*cosa - sina*sina)*Dxy;

      if( Vsec )  Vsec[n] = acu * HUs;
    }
  }

  MEMORY::memo.Detach( curv );


  char text[80];

  sprintf( text, "\n (GRID::Dispersion)      %s\n",
                 "Dispersion coefficients determined" );
  REPORT::rpt.Output( text, 4 );
}
