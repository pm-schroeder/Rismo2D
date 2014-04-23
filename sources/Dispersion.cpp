// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class GRID
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
