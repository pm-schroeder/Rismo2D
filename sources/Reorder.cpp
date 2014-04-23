// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class REORDER
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

// -------------------------------------------------------------------------------------------------
// Reordering algorithm for elements by KING (1970)

#include "Defs.h"
#include "Report.h"
#include "Node.h"
#include "Elem.h"
#include "Grid.h"

#include "Reorder.h"


REORDER::REORDER()
{
}


REORDER::REORDER( GRID* region )
{
  // initialization of public variables
  max      = 0;
  square   = 0.0;

  minDepth = 0;      // minDepth and maxDepth may be
  maxDepth = 0;      // changed after construction

  // initialization of private variables
  frnt     = NULL;
  num      = 0;
  hist_cnt = 1;
  prev_el  = NULL;

  // allocate memory for FRONT structure
  frontPtr = new REOFRONT [region->Getnp()];

  if ( !frontPtr )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (REORDER::REORDER - 1)" );

  // set up node pointers in FRONT structure
  for( int n=0; n<region->Getnp(); n++ )
  {
    frontPtr[n].node = region->Getnode(n);
  }

  // initialize elements
  nel = 0;

  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* el = region->Getelem(e);

    nel++;
    el->mark = false;
  }
}


REORDER::~REORDER()
{
  delete[] frontPtr;
}


void REORDER::initFront( ELEM* elem )
{
  num = updateFront( elem, &frnt, num, 1 );
  if( num > max )  max = num;
  square += num*num;
  nel--;
  elem->mark = true;
  if( prev_el ) prev_el = elem;
  prev_el = elem;
  elem->link = (ELEM*) 0;
}


void REORDER::start()
{
  int       min;             // minimale Anzahl aktiver Knoten
  int       size;            // Anzahl der Elemente an der Front
  REOFRONT* fr;              // Listenzeiger fuer aktiven Knoten
  ELEM**    list;            // Liste der Elemente an der Front
  ELEM*     elem;            // ausgewaehltes Element aus *list


  if( !prev_el )
  {
    REPORT::rpt.Error( "no starting element for reordering (REORDER::start - 1)" );
  }


  for( ;; )
  {
     if( num > max )  max = num;

     if( nel == 0 )
       REPORT::rpt.Error( "unexpected error (nel == 0) (REORDER::start - 2)" );

     if( num == 0 )
     {
       REPORT::rpt.Error( "\nunexpected error (num == 0) (REORDER::start - 3)\n" );
     }


     // find number of elements at actual front

     size = 0;
     fr   = frnt;

     while ( fr )
     {
        size += fr->node->noel;
        fr = fr->next;
     }

     list = new ELEM* [size];

     if ( !list )
       REPORT::rpt.Error( "can not allocate memory (REORDER::start - 4)" );


     // create list of front elements

     size = frontListOfElems (frnt, list);

     if ( size == 0 )
     {
       delete[] list;
       REPORT::rpt.Error ("Unexpected error - reorder (3)");
     }

     elem = list[0];

     if ( nel == 1 )
     {
       prev_el->link = elem;
       elem->link = (ELEM *)0;
       delete[] list;

       REPORT::rpt.Screen( 5, "  element count down %7d;  front width: 0\r", nel );
       break;
     }

     if ( size > 1 )
     {
       register int i;
       int depth, m, found;
       double sqr, min_sqr;     // square sum of front width
       ELEM *el;

       depth = minDepth;

       do
       {
         found   = -1;
         m       = 0;
         min_sqr = 0.0;

         for (i=0; i<size; i++)
         {
           int n;

           if (!list[i])  continue;

           found++;

           // Front an die Hinzunahme des Elementes list[i] anpassen
           min = updateFront (list[i], &frnt, num, 0);

           // maximale Frontweite bis zur Rekursionstiefe depth
           sqr = square + min*min;
           n = frontWidth (&frnt, min, depth, &sqr);

           removeFront (list[i], &frnt, min);

           if (min > n)     n = min;
           if (!m || n < m) m = n;

           if (min_sqr>0.1 && sqr>min_sqr)
           {
             list[i] = (ELEM *)0;
           }

           else if (min_sqr<0.1 || sqr<min_sqr)
           {
             register int j;

             // wenn list[i] zu einer guenstigeren Frontweite fuehrt ...
             min_sqr = sqr;
             elem = list[i];
             for (j=0; j<i; j++)
               list[j] = (ELEM *)0;
           }
         }
         depth++;

       } while ( found > 0  &&  depth <= maxDepth);

       // aus allen guenstigen Frontelementen dasjenige waehlen, das
       // am laengsten in der Liste steht

       if( found > 1 )  el = chooseElem( list, size );

       if( found > 1  &&  el )  elem = el;
     }

     delete[] list;

     prev_el->link = elem;
     elem->link = (ELEM *)0;

     // Front an die Hinzunahme des ausgewaehlten Elementes elem anpassen
     min = updateFront (elem, &frnt, num, 1);

     REPORT::rpt.Screen( 5, "  element count down %7d;  front width: %d \r", nel, min );

     nel--;
     prev_el = elem;
     num     = min;
     square += num*num;


     // und dann weiter im 'reordering'
  }
}


int REORDER::frontListOfElems( REOFRONT *fr, ELEM  **list )
{
   int size = 0;

   while( fr )
   {
     // Schleife ueber alle Elemente die an dem Fronknoten fr liegen

     int noel = fr->node->noel;

     for( int i=0; i<noel; i++ )
     {
      ELEM *el = (ELEM *) fr->node->el[i];

      // Element bereits erfasst ?
      if( el->mark )  continue;

      // Element steht bereits in der Liste ?
      for( int j=0; j<size; j++ )
      {
        if( el == list[j] )
        {
          el = (ELEM *) 0;
          break;
        }
      }

      if( el )
      {
        list[size] = el;
        size++;
      }
    }

    fr = fr->next;
  }

  return size;
}


int REORDER::updateFront( ELEM *el, REOFRONT **frnt, int num, int histflg )
{
   int i, ncn;

   ncn = el->Getncn();

   el->mark = true;

   for( i=0; i<ncn; i++ )
   {
      NODE *nd;
      REOFRONT *fr, *fr_prev;

      nd = el->nd[i];

      fr      = *frnt;
      fr_prev = NULL;

      while( fr )
      {
         if( fr->node == nd )  break;

         fr_prev = fr;
         fr      = fr->next;
      }


      // wenn der Knoten bereits in der Liste steht ...

      if ( fr )
      {
         fr->miss--;

         if ( !fr->miss )
         {
            // Knoten ist mit dem Element 'el' komplett und
            // wird aus der Liste aktiver Knoten entfernt

            num--;
            if (histflg)
               fr->hist = 0;

            if (fr_prev)
               fr_prev->next = fr->next;
            else
               *frnt = fr->next;
         }
      }


      // ... andernfalls wird er an die Liste angehaengt

      else
      {
         if( nd->noel > 1)
         {
            fr = &frontPtr[nd->Getno()];

            if( fr_prev )
               fr_prev->next = fr;
            else
               *frnt = fr;

            fr->next = NULL;
            fr->miss = nd->noel - 1;

            if( histflg )
               fr->hist = hist_cnt++;

            num++;
         }
      }
   }

   return (num);
}


int REORDER::removeFront( ELEM *el, REOFRONT **frnt, int num )
{
   int i, ncn;

   ncn = el->Getncn();

   el->mark = false;

   for( i=0; i<ncn; i++ )
   {
      NODE *nd;
      REOFRONT *fr, *fr_prev;

      nd = el->nd[i];

      fr      = *frnt;
      fr_prev = NULL;

      while( fr )
      {
         if( fr->node == nd )  break;

         fr_prev = fr;
         fr      = fr->next;
      }

      // wenn der Knoten in der Liste steht ...

      if( fr )
      {
         fr->miss++;

         if( fr->miss == nd->noel )
         {
            // Knoten war mit dem Element 'el' neu hinzugekommen
            // und kann nun wieder entfernt werden
            num--;

            if( fr_prev )
               fr_prev->next = fr->next;
            else
               *frnt = fr->next;
         }
      }

      // ... andernfalls wird er nun wieder aktiv

      else
      {
         if( nd->noel > 1 )
         {
            fr = &frontPtr[nd->Getno()];

            if( fr_prev )
               fr_prev->next = fr;
            else
               *frnt = fr;

            fr->next = NULL;
            fr->miss = 1;
            num++;
         }
      }
   }

   return num;
}


int REORDER::frontWidth( REOFRONT **frnt, int num, int depth, double *square )
{
   int       i;
   int       dp, m, size, min;
   double    sqr, min_sqr;
   REOFRONT* fr;
   ELEM**    list;

   if( num == 0  ||  depth == 0 )
      return num;


   // Anzahl der Elemente an der Front bestimmen

   size = 0;
   fr   = *frnt;

   while( fr )
   {
      size += fr->node->noel;
      fr = fr->next;
   }

   list = new ELEM* [size];

   if( !list )
     REPORT::rpt.Error ("Could not allocate memory - frontWidth (1)");

   // Liste der Frontelemente erstellen, die noch nicht erfasst wurden
   size = frontListOfElems( *frnt, list );

   if( size==0 )
   {
      delete[] list;
      REPORT::rpt.Error ("Unexpected error -  frontWidth (2)");
   }

   dp = (minDepth<depth)? (minDepth):(depth-1);

   do
   {
      m       = 0;
      min_sqr = 0.0;

      for( i=0; i<size; i++ )
      {
         int n;

         if( !list[i] )  continue;

         // Front an die Hinzunahme des Elementes list[0] anpassen
         min = updateFront( list[i], frnt, num, 0 );

         // und dann rekursiv weiter
         sqr = *square + min*min;
         n = frontWidth( frnt, min, dp, &sqr );

         removeFront( list[i], frnt, min );

         if( min > n )     n = min;
         if( !m || n < m ) m = n;

         if( min_sqr>0.1 && sqr>min_sqr )
         {
            list[i] = (ELEM *) 0;
         }

         else if( min_sqr<0.1 || sqr<min_sqr )
         {
            register int j;

            // wenn list[i] zu einer guenstigeren Frontweite fuehrt ...
            min_sqr = sqr;

            for( j=0; j<i; j++ )
               list[j] = (ELEM *)0;
         }
      }
      dp++;

   } while( dp<depth );

   delete[] list;

   *square = min_sqr;

   return m;
}


// aus der Liste der naechstbesten Elemente wird nun dasjenige heraus-
// gesucht, dessen Knoten am laengsten in der Frontliste stehen
//
// Aenderung: ..., dessen Knoten die kleinste Nummer hat

ELEM* REORDER::chooseElem( ELEM **list, int size )
{
   ELEM* elem = (ELEM*) 0;
   int   less = 0;


   // Schleife ueber die moeglichen Elemente

   for( int i=0; i<size; i++ )
   {
      if( !list[i] ) continue;


      // Element mit der kleinsten Nummer suchen

      if( !less  ||  list[i]->Getname() < less )
      {
        elem = list[i];
        less = elem->Getname();
      }

/*
      int ncn = list[i]->shsp->sh[0]->nnd;

      for( int j=0; j<ncn; j++ )
      {
         // wenn's den Knoten [j] gibt und er in der Frontliste steht ...
         // (Knoten in der Frontliste haben ein fr->hist!=0 !!)

         NODE* nd = list[i]->nd[j];

         REOFRONT* fr = &frontPtr[nd->no];

         if( fr->hist )
         {
            if( !less  ||  fr->hist < less )
            {
               less = fr->hist;
               elem = list[i];
            }
         }
      }
*/
   }

   return elem;
}
