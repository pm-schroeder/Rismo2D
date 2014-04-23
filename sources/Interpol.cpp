// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class DRYREW
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
#include "Node.h"
#include "Elem.h"

#include "Model.h"


double DRYREW::interpolate( NODE *node, int var, int depth )
{
  double val = 0.0;
  double dis = 0.0;

  int k    = 0;
  int noel = node->noel;

  for( int i=0; i<noel; i++ )
  {
    ELEM* el = (ELEM *) node->el[i];

    int ncn = el->Getncn();

    if( isFS(el->flag, ELEM::kRegion) )
    {
      for( int j=0; j<ncn; j++ )
      {
        if( depth > 0 )
        {
          if( node != el->nd[j] )
          {
            double dx = node->x - el->nd[j]->x;
            double dy = node->y - el->nd[j]->y;

            double d = sqrt( dx*dx + dy*dy );

            double v = interpolate( el->nd[j], var, depth-1 );

            if( fabs(v) > 1.0e-30 )
            {
              val +=   v / d;
              dis += 1.0 / d;

              k++;
            }
          }
        }

        else
        {
          double U = el->nd[j]->v.U;
          double V = el->nd[j]->v.V;

          double S = el->nd[j]->v.S;
          double A = el->nd[j]->z;
          double H = S - A;

          double v = 0.0;

          switch( var )
          {
            case kVarU:   v = U;                    break;
            case kVarV:   v = V;                    break;
            case kVarH:   v = H;                    break;
            case kVarS:   v = S;                    break;
            case kVarK:   v = el->nd[j]->v.K;       break;
            case kVarD:   v = el->nd[j]->v.D;       break;

            case kVarUH:  v = U * H;                break;
            case kVarVH:  v = V * H;                break;
          }

          if( node != el->nd[j]  &&  !isFS(el->nd[j]->flag, NODE::kDry) )
          {
            k++;

            double dx = node->x - el->nd[j]->x;
            double dy = node->y - el->nd[j]->y;

            double d = sqrt( dx*dx + dy*dy );

            val +=   v / d;
            dis += 1.0 / d;
          }
        }
      }
    }
  }

  if( k )  return val / dis;
  else     return 0.0;
}
