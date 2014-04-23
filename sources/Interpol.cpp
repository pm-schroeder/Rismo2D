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
