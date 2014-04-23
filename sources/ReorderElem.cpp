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
#include "Report.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Section.h"
#include "Reorder.h"

#include "Model.h"


ELEM* MODEL::ReorderElem( int ns, SECTION* section )
{
  // compute the center of all elements

  double* xcenter = (double*) MEMORY::memo.Array_el( region->Getne() );
  double* ycenter = (double*) MEMORY::memo.Array_el( region->Getne() );

  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* rElem = region->Getelem(e);

    rElem->center( &xcenter[rElem->Getno()], &ycenter[rElem->Getno()] );
    rElem->mark = false;
  }


  // reorder elements

  REPORT::rpt.Screen( 5, "\n\n" );

  int count = 1;

  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* rElem = region->Getelem(e);

    if( (count%50) == 0  ||  e == (region->Getne()-1) )
      REPORT::rpt.Screen( 5, "  working on element %7d\r", count );

    count++;

    double extent = 1.5 * rElem->area() / rElem->perimeter();

    if( isFS(rElem->flag, ELEM::kRegion) )
    {
      // check, if element is upstream to the first section

      double x = xcenter[rElem->Getno()];
      double y = ycenter[rElem->Getno()];

      double xl = section[0].tangentialDistance(x, y);
      double yl = section[0].normalDistance(x, y);

//    if( yl < 0.0 )
      if( xl>=0.0 && xl<=section[0].length() && yl<0.0  )
      {
        ELEM* listE = section[0].elem;
        ELEM* prevE = (ELEM*) 0;

        while( listE )
        {
          double xe, ye;

          xe = xcenter[listE->Getno()];
          ye = ycenter[listE->Getno()];

          double xel = section[0].tangentialDistance(xe, ye);
          double yel = section[0].normalDistance(xe, ye);

          if(     (yl < yel-extent)
              ||  (fabs(yl-yel) < extent  &&  xl < xel) )  break;

          prevE = listE;
          listE = listE->link;
        }

        if( prevE )
        {
          rElem->link = listE;
          prevE->link = rElem;
        }

        else
        {
          section[0].elem = rElem;
          rElem->link     = listE;
        }

        goto endOfLoop;
      }


      // search in following sections

      for( int i=1; i<ns; i++ )
      {
        // check, if element is upstream to section i

        double d = section[i].normalDistance(x, y);
        double L = section[i].tangentialDistance(x, y);


//      if( d < 0.0 )
        if( L>=0.0 && L<=section[i].length() && d<0.0 )
        {
          xl = section[i-1].tangentialDistance(x, y);
          yl = section[i-1].normalDistance(x, y);

          ELEM* listE = section[i-1].elem;
          ELEM* prevE = (ELEM*) 0;

          while( listE )
          {
            double xe = xcenter[listE->Getno()];
            double ye = ycenter[listE->Getno()];

            double xel = section[i-1].tangentialDistance(xe, ye);
            double yel = section[i-1].normalDistance(xe, ye);

            if(     (yl < yel-extent)
                ||  (fabs(yl-yel) < extent  &&  xl < xel) )  break;

            prevE = listE;
            listE = listE->link;
          }

          if( prevE )
          {
            rElem->link = listE;
            prevE->link = rElem;
          }

          else
          {
            section[i-1].elem = rElem;
            rElem->link       = listE;
          }

          goto endOfLoop;
        }
      }


      // element is downstream to the last section

      xl = section[ns-1].tangentialDistance(x, y);
      yl = section[ns-1].normalDistance(x, y);

      ELEM* listE = section[ns-1].elem;
      ELEM* prevE = (ELEM*) 0;

      while( listE )
      {
        double xe, ye;

        xe = xcenter[listE->Getno()];
        ye = ycenter[listE->Getno()];

        double xel = section[ns-1].tangentialDistance(xe, ye);
        double yel = section[ns-1].normalDistance(xe, ye);

        if(     (yl < yel-extent)
            ||  (fabs(yl-yel) < extent  &&  xl < xel) )  break;

        prevE = listE;
        listE = listE->link;
      }

      if( prevE )
      {
        rElem->link = listE;
        prevE->link = rElem;
      }

      else
      {
        section[ns-1].elem = rElem;
        rElem->link        = listE;
      }
    }

    endOfLoop:  continue;
  }


  // create new list of elements

  ELEM* first = (ELEM*) 0;
  ELEM* listE = (ELEM*) 0;

  for( int i=0; i<ns; i++ )
  {
    if( listE )  listE->link   = section[i].elem;
    else         first = listE = section[i].elem;

    if( listE )  while( listE->link )  listE = listE->link;
  }


  // renumber elements

  int eno = 0;
  listE = first;

  while( listE )
  {
    listE->Setno( eno );
    listE->Setname( eno + 1 );
    eno++;

    listE = listE->link;
  }


  REPORT::rpt.Screen( 5, "\n\n" );


  // reorder elements according to KING's procedure

  this->list = first;

  REORDER reorder( region );

  reorder.initFront( first );
  reorder.start();


  // renumber elements

  eno   = 0;
  listE = first;

  while( listE )
  {
    listE->Setno( eno );
    listE->Setname( eno + 1 );
    eno++;

    listE = listE->link;
  }


  // create a new element list

  int new_ne = region->Getne();

  ELEM* newElem = new ELEM [new_ne];

  if( !newElem )
    REPORT::rpt.Error( "can not allocate memory - reorderElem - 2" );

  listE = first;

  for( int i=0; i<new_ne; i++ )
  {
    newElem[i] = *listE;
    listE = listE->link;
  }

  region->KillElem();
  region->Setelem( new_ne, newElem );

  Initialize();

  region->Connection( 0l );


  // determine averaged and maximum front width of corner nodes

  REPORT::rpt.Screen( 5, "\n\n" );

  LastNode();

  for( int i=0; i<np; i++ )  node[i]->mark = false;

  int    fw    = 0;
  int    maxFW = 0;
  double sumFW = 0.0;

  for( int i=0; i<ne; i++ )
  {
    ELEM* E = elem[i];

    int ncn = E->Getncn();

    for( int j=0; j<ncn; j++ )
    {
      if( !E->nd[j]->mark )
      {
        E->nd[j]->mark = true;
        fw++;
      }
    }

    if( fw > maxFW )  maxFW = fw;

    sumFW += fw;

    REPORT::rpt.Screen( 1, "  front width is       %5d;  maximum: %5d\r", fw, maxFW );

    for( int j=0; j<ncn; j++ )
    {
      if( E->isLast[0][j] )
      {
        E->nd[j]->mark = false;
        fw--;
      }
    }
  }

  sumFW /= region->Getne();

  REPORT::rpt.Screen( 1, "\n\n" );

  REPORT::rpt.Message( 1, "\n" );
  REPORT::rpt.Message( 1, "%-25s%s\n", " (MODEL::ReorderElem)",
                                    "front width of corner nodes" );
  REPORT::rpt.Message( 1, "%25s%s%d\n",      " ", "maximum: ", maxFW );
  REPORT::rpt.Message( 1, "%25s%s%-8.2lf\n", " ", "average: ", sumFW );

  REPORT::rpt.Screen( 1, "\n\n" );


  MEMORY::memo.Detach( xcenter );
  MEMORY::memo.Detach( ycenter );

  return first;
}
