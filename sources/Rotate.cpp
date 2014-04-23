// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS
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
#include "Bcon.h"
#include "Node.h"

#include "Eqs.h"


void EQS::Rotate2D( int nnd, NODE** nd, int df, double** estifm, double* force )
{
  for( int i=0; i<nnd; i++ )
  {
    BCON* bcon = &nd[i]->bc;

    if( isFS(nd[i]->flag, NODE::kRotat) )
    {
      double rot[2][2];

      rot[0][0] = bcon->Getrot(0,0);
      rot[0][1] = bcon->Getrot(0,1);
      rot[1][0] = bcon->Getrot(1,0);
      rot[1][1] = bcon->Getrot(1,1);

      double det = rot[0][0]*rot[1][1] - rot[0][1]*rot[1][0];

      if( fabs(det) < 1.0e-12 )
      {
        if( nnd > 3 )
        {
          // print this warning only for region elements, to reduce the amount of output------------
          REPORT::rpt.Warning( kGridFault, "determinant of rotation matrix very small (%le)", det );
          REPORT::rpt.Warning( 0, "no rotation of equation system for node %d", nd[i]->Getname() );
        }

        // do not rotate estifm and force, jump to the next node of elem ---------------------------
        continue;
      }


      int indU = i;
      int indV = nnd  +  i;


      if( estifm )
      {
        // multiply estifm[][] * rot[][] -----------------------------------------------------------

        for( int j=0; j<maxEleq; j++ )
        {
          double tmp      = estifm[j][indU] * rot[0][0]  +
                            estifm[j][indV] * rot[1][0];
          estifm[j][indV] = estifm[j][indU] * rot[0][1]  +
                            estifm[j][indV] * rot[1][1];
          estifm[j][indU] = tmp;
        }


        // multiply transpose(rot[][]) * estifm[][] ------------------------------------------------

        for( int j=0; j<maxEleq; j++ )
        {
          double tmp      = rot[0][0] * estifm[indU][j]  +
                            rot[1][0] * estifm[indV][j];
          estifm[indV][j] = rot[0][1] * estifm[indU][j]  +
                            rot[1][1] * estifm[indV][j];
          estifm[indU][j] = tmp;
        }
      }


      if( force )
      {
        // force vector: multiply transpose(rot[][]) * force[] -------------------------------------

        double tmp  = rot[0][0] * force[indU]  +  rot[1][0] * force[indV];
        force[indV] = rot[0][1] * force[indU]  +  rot[1][1] * force[indV];
        force[indU] = tmp;
      }
    }
  }
}
