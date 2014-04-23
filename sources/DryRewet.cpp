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

#include "Node.h"
#include "Shape.h"
#include "Type.h"
#include "Elem.h"
#include "Project.h"

#include "Grid.h"

//#define kDebug


// ---------------------------------------------------------------------------------------
// mark dry nodes and elements
// ---------------------------------------------------------------------------------------

int GRID::Dry( double  dryLimit,
               int     countDown )
{
  // initialization: mark all dry nodes and elements -------------------------------------

  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    if( isFS(nd->flag, NODE::kDry) )  nd->mark = true;
    else                              nd->mark = false;

    SF( nd->flag, NODE::kDry );
  }


  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if ( isFS(el->flag, ELEM::kDry) )  el->mark = true;
    else                               el->mark = false;

    CF( el->flag, ELEM::kDry );
  }


  // loop on all elements ----------------------------------------------------------------

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    int ncn = el->Getncn();

    for( int i=0; i<ncn; i++ )
    {
      // set element dry, if flow depth is less than dry limit ---------------------------
      // or the node was already dry

      double H = el->nd[i]->v.S - el->nd[i]->z;

      if( H < dryLimit  ||  el->nd[i]->mark )
      {
        SF( el->flag, ELEM::kDry );
      }
    }
  }


  // loop on all elements: all nodes at wet elements are wet -----------------------------

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( !isFS(el->flag, ELEM::kDry) )
    {
      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )  CF( el->nd[i]->flag, NODE::kDry );
    }
  }


  // report all nodes that have got newly dry --------------------------------------------

  REPORT::rpt.Output( "\n (dry)           the following nodes have got dry\n", 5 );

  int jdry = 0;

  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    if( isFS(nd->flag, NODE::kDry) )
    {
      nd->v.U = 0.0;
      nd->v.V = 0.0;
      nd->v.S = nd->z;

      nd->v.K = 0.0;
      nd->v.D = 0.0;

      nd->v.dUdt = 0.0;
      nd->v.dVdt = 0.0;

      nd->v.dSdt = 0.0;

      if( !nd->mark )
      {
        nd->countDown = countDown;

        char text[20];
        sprintf( text, "  %5d", nd->Getname() );
        REPORT::rpt.Output( text, 2 );

        jdry++;
        if( !( jdry % 10) )  REPORT::rpt.Output( "\n", 5 );
      }
    }

    nd->mark = false;
  }

  if( jdry % 10 )  REPORT::rpt.Output( "\n", 5 );


  // report all elements that have got newly dry -------------------------------------------

  REPORT::rpt.Output( "\n (dry)           the following elements have got dry\n", 5 );

  int idry = 0;

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem( e );

    if( isFS(el->flag, ELEM::kDry) && !el->mark )
    {
      char text[20];
      sprintf( text, "  %5d", el->Getname() );
      REPORT::rpt.Output( text, 5 );

      idry++;
      if( !( idry % 10) )  REPORT::rpt.Output( "\n", 5 );
    }

    el->mark = false;
  }

  if( idry % 10 )  REPORT::rpt.Output( "\n", 5 );


  return idry + jdry;
}


// ---------------------------------------------------------------------------------------
// mark nodes and elements to be rewetted
// ---------------------------------------------------------------------------------------

int GRID::Rewet( double   rewetLimit,
                 int      rewetPasses,
                 PROJECT* project )
{
  int i, j, k;
  int l, m, pass, ncn, nnd, rewetFlag;
  double wElev, K, D, cm, cd, vt;
  char text[80];

  cm = project->KD.cm;
  cd = project->KD.cd;


  // initialization: mark all dry elements -----------------------------------------------

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    el->mark = false;

    if( isFS(el->flag, ELEM::kDry) )
    {
      el->mark = true;
    }
  }


  REPORT::rpt.Output( "\n (rewet)         the following nodes have been rewetted\n", 5 );

  m = 0;

  for( pass=0; pass<rewetPasses; pass++ )
  {
    for( int n=0; n<Getnp(); n++ )
    {
      NODE* nd = Getnode(n);

      nd->mark = false;
    }

    sprintf( text, " pass %d ...\n", pass+1 );
    REPORT::rpt.Output( text, 5 );


    for( int e=0; e<Getne(); e++ )
    {
      ELEM * el = Getelem(e);

      if( isFS(el->flag, ELEM::kBound) )  continue;


      TYPE* type = TYPE::Getid( el->type );


      // look only for dry elements ------------------------------------------------------

      if( isFS(el->flag, ELEM::kDry) )
      {
        // loop on corner nodes: minimal water elevation ---------------------------------

        ncn = el->Getncn();        // number of corner nodes
        nnd = el->Getnnd();        // total number of nodes

        rewetFlag = false;

        wElev = 0.0;

        for( j=0, i=0; j<ncn; j++ )
        {
          if( !isFS(el->nd[j]->flag, NODE::kDry) )
          {
            // determine averaged water elevation at wet points --------------------------

            rewetFlag = true;

            wElev += el->nd[j]->v.S;
            i++;
          }
        }


        // loop on all nodes: check for rewetting ----------------------------------------

        if( i )
        {
          wElev /= i;

          for( j=0; j<nnd; j++ )
          {
            if( (wElev - el->nd[j]->z) < rewetLimit )  rewetFlag = false;
          }
        }


        // if all nodes have been rewetted, rewet whole element --------------------------

        if( rewetFlag )
        {
          for( j=0; j<ncn; j++ )
          {
            // re-initialize flow parameters ---------------------------------------------

            if( isFS(el->nd[j]->flag, NODE::kDry) )
            {
              el->nd[j]->mark = true;
              el->nd[j]->v.S = wElev;
            }
          }

          for( j=ncn; j<nnd; j++ )
          {
            if( isFS(el->nd[j]->flag, NODE::kDry) )
            {
              el->nd[j]->mark = true;

              // get left and right corner node to midside node j ------------------------
              el->GetLShape()->getCornerNodes( j, &i, &k );

              el->nd[j]->v.S = (el->nd[i]->v.S + el->nd[k]->v.S) / 2.0;
            }
          }

          for( j=0; j<nnd; j++ )
          {
            if( isFS(el->nd[j]->flag, NODE::kDry) )
            {
              el->nd[j]->v.U = 0.0;
              el->nd[j]->v.V = 0.0;

              el->nd[j]->v.dUdt = 0.0;
              el->nd[j]->v.dVdt = 0.0;

              el->nd[j]->v.dSdt = 0.0;


              K = el->nd[j]->v.K = dryRew.interpolate(el->nd[j], kVarK, 1);
              D = el->nd[j]->v.D = dryRew.interpolate(el->nd[j], kVarD, 1);

              if( isFS(project->actualTurb, BCONSET::kVtConstant) )
              {
                el->nd[j]->vt = type->vt;
              }
              else if( isFS(project->actualTurb, BCONSET::kVtPrandtlKol) )
              {
                if( fabs(D) > 0.0 )  vt = cm * cd * K * K / D;
                else                 vt = 0.0;

                el->nd[j]->vt = vt;
              }
            }
          }
        }
      }
    }


    // check count down for nodes to be rewetted and report ------------------------------

    l = 0;
    for( int n=0; n<Getnp(); n++ )
    {
      NODE* nd = Getnode(n);

      if( nd->mark )
      {
        nd->countDown--;

        if( nd->countDown <= 0 )
        {
          CF( nd->flag, NODE::kDry );


          // look for boundary conditions ------------------------------------------------

          BCON *bc = &nd->bc;


          // ... outflow boundary or fixed flow depth h --------------------------------

          if( isFS(bc->kind, BCON::kOutlet)  ||  isFS(bc->kind, BCON::kSetS) )
          {
            nd->v.S = bc->val->S;
          }


          // ... inflow boundary or fixed velocities U,V -------------------------------

          if( isFS(bc->kind, BCON::kInlet) )
          {
            nd->v.U = bc->val->U * bc->niox / (nd->v.S - nd->z);
            nd->v.V = bc->val->U * bc->nioy / (nd->v.S - nd->z);
          }

          else if( !isFS(bc->kind, BCON::kAutoSlip) )
          {
            if ( isFS(bc->kind, BCON::kFixU) )  nd->v.U = bc->val->U;
            if ( isFS(bc->kind, BCON::kFixV) )  nd->v.V = bc->val->V;
          }

#         ifdef kDebug
          sprintf( text,
                   "  %5d (UVhKDvt): %9.2le %9.2le %9.2le %9.2le %9.2le %9.2le\n",
                                  nd->Getname(),
                                  nd->v.U,
                                  nd->v.V,
                                  nd->v.S - nd->z,
                                  nd->v.K,
                                  nd->v.D,
                                  nd->vt );
          REPORT::rpt.Output( text, 5 );

#         else

          sprintf( text, "  %5d", nd->Getname() );
          REPORT::rpt.Output( text, 5 );

          l++;
          if( !( l % 10) )  REPORT::rpt.Output( "\n", 5 );
#         endif

          m++;
        }

        else
        {
          nd->v.U = 0.0;
          nd->v.V = 0.0;
          nd->v.K = 0.0;
          nd->v.D = 0.0;
          nd->v.S = nd->z;
        }
      }

      nd->mark = false;
    }

#   ifndef kDebug
    if( l % 10 )
      REPORT::rpt.Output( "\n", 5 );
#   endif


    // rewet all dry elements without dry nodes ------------------------------------------

    for( int e=0; e<Getne(); e++ )
    {
      ELEM* el = Getelem(e);

      if( isFS(el->flag, ELEM::kDry) )
      {
        nnd = el->Getnnd();

        rewetFlag = true;

        for( j=0; j<nnd; j++ )
        {
          if( isFS(el->nd[j]->flag, NODE::kDry) )
          {
            rewetFlag = false;
            break;
          }
        }

        if( rewetFlag )
        {
          CF( el->flag, ELEM::kDry );
        }
      }
    }
  }


  // loop on all elements: all nodes at wet elements are wet -----------------------------

  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    SF( nd->flag, NODE::kDry );

    nd->mark = false;
  }


  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( !isFS(el->flag, ELEM::kDry) )
    {
      nnd = el->Getnnd();

      for( j=0; j<nnd; j++ )  CF( el->nd[j]->flag, NODE::kDry );
    }
  }


  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    if( isFS(nd->flag, NODE::kDry) )
    {
      nd->v.U = 0.0;
      nd->v.V = 0.0;
      nd->v.K = 0.0;
      nd->v.D = 0.0;
      nd->v.S = nd->z;
    }
  }



  // report elements which have been rewetted --------------------------------------------

  REPORT::rpt.Output( "\n (rewet)         the following elements have been rewetted\n", 5 );

  l  = 0;
  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    nnd = el->Getnnd();

    if( el->mark  &&  !isFS(el->flag, ELEM::kDry) )
    {
      for( j=0; j<nnd; j++ )
      {
        if( isFS(el->nd[j]->flag, NODE::kDry) )
          REPORT::rpt.Error ("unexpected error - rewet (1)");
      }

      sprintf( text, "  %5d", el->Getname() );
      REPORT::rpt.Output( text, 5 );

      l++;
      if( !( l % 10) )  REPORT::rpt.Output( "\n", 5 );
    }

    el->mark = false;
  }

  if( l % 10 )  REPORT::rpt.Output( "\n", 5 );


# ifndef kDebug
  if( l % 10 )  REPORT::rpt.Output( "\n", 5 );
# endif


  return l;
}


void GRID::ReportDry( PROJECT* project,
                      double   dryLimit,
                      int      countDown )
{
  int dryNodes = false;


  // check for dry nodes -------------------------------------------------------------------

  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    nd->countDown = 0;

    CF( nd->flag, NODE::kDry );

    if( (nd->v.S - nd->z) < dryLimit )
    {
      SF( nd->flag, NODE::kDry );

      dryNodes = true;
    }
  }


  // report all dry nodes ------------------------------------------------------------------

  if( dryNodes )
  {
    REPORT::rpt.Output( "\n\n----------------------------------------", 5 );
    REPORT::rpt.Output( "------------------------------\n", 5 );
    REPORT::rpt.Output( "\n ... the following nodes are dry\n", 5 );

    int jdry = 0;

    for( int n=0; n<Getnp(); n++ )
    {
      NODE* nd = Getnode(n);

      if( isFS(nd->flag, NODE::kDry) )
      {
        char text[20];
        sprintf( text, "  %5d", nd->Getname() );
        REPORT::rpt.Output( text, 5 );

        jdry++;
        if( !( jdry % 10) ) REPORT::rpt.Output( "\n", 5 );
      }
    }

    if( jdry % 10 ) REPORT::rpt.Output( "\n", 5 );


    // loop on all elements ----------------------------------------------------------------

    for( int e=0; e<Getne(); e++ )
    {
      ELEM* el = Getelem(e);

      CF( el->flag, ELEM::kDry );

      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        if( isFS(el->nd[i]->flag, NODE::kDry) )
        {
          SF( el->flag, ELEM::kDry );
          break;
        }
      }
    }


    for( int n=0; n<Getnp(); n++ )
    {
      NODE* nd = Getnode(n);

      SF( nd->flag, NODE::kDry );
    }


    // loop on all elements: all nodes at wet elements are wet -----------------------------

    for( int e=0; e<Getne(); e++ )
    {
      ELEM* el = Getelem(e);

      if( !isFS(el->flag, ELEM::kDry) )
      {
        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ ) CF( el->nd[i]->flag, NODE::kDry );
      }
    }


    int idry = 0;

    REPORT::rpt.Output( "\n ... the following elements are dry\n", 5 );

    for( int e=0; e<Getne(); e++ )
    {
      ELEM* el = Getelem(e);

      if( isFS(el->flag, ELEM::kDry) )
      {
        char text[20];
        sprintf( text, "  %5d", el->Getname() );
        REPORT::rpt.Output( text, 5 );

        idry++;
        if( !( idry % 10) )  REPORT::rpt.Output( "\n", 5 );
      }
    }

    if( idry % 10 )  REPORT::rpt.Output( "\n", 5 );

    REPORT::rpt.Output( "\n", 5 );


    for( int n=0; n<Getnp(); n++ )
    {
      NODE* nd = Getnode(n);

      if( isFS(nd->flag, NODE::kDry) )
      {
        nd->v.U = 0.0;
        nd->v.V = 0.0;
        nd->v.S = nd->z;

        nd->v.K = 0.0;
        nd->v.D = 0.0;

        nd->v.dUdt = 0.0;
        nd->v.dVdt = 0.0;

        nd->v.dSdt = 0.0;
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////
// mark dry nodes and elements (method 2)
// this method will keep dried nodes dry until countDown counts to zero
//////////////////////////////////////////////////////////////////////////////////////////

void GRID::DryRewet( double  dryLimit,
                     double  rewetLimit,
                     int     countDown,
                     int*    dried,
                     int*    wetted )
{
  // -------------------------------------------------------------------------------------
  // check for dry nodes

  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    // mark recently dry nodes
    if( isFS(nd->flag, NODE::kDry) )  nd->mark = true;
    else                              nd->mark = false;

    CF( nd->flag, NODE::kMarsh );


    double H = nd->v.S - nd->zor;

    if( H < dryLimit )
    {
      SF( nd->flag, NODE::kDry );
      nd->countDown = countDown;
    }
    else
    {
      if( nd->countDown > 0 )  nd->countDown--;
      else                     CF( nd->flag, NODE::kDry );
    }

    nd->z = nd->zor;
  }


  // -------------------------------------------------------------------------------------
  // loop on all elements: check for dry elements

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    // mark elements that have been dry --------------------------------------------------
    if( isFS(el->flag, ELEM::kDry) )  el->mark = true;
    else                              el->mark = false;

    CF( el->flag, ELEM::kDry | ELEM::kMarsh );

    int ncn = el->Getncn();

    int ndry = 0;

    for( int i=0; i<ncn; i++ )
    {
      if( isFS(el->nd[i]->flag, NODE::kDry) )  ndry++;
    }


    // set dry flag if all nodes are dry
    // set marsh flag if at least one node is dry

    if( ndry == ncn )
    {
      SF( el->flag, ELEM::kDry );
    }

    else if( ndry )
    {
      SF( el->flag, ELEM::kMarsh );

      for( int i=0; i<ncn; i++ )
      {
        if( isFS(el->nd[i]->flag, NODE::kDry) )
        {
          SF( el->nd[i]->flag, NODE::kMarsh );
        }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  // loop on all elements: all nodes at wet elements are wet

  for( int n=0; n<Getnp(); n++ )
  {
    SF( Getnode(n)->flag, NODE::kDry );
  }


  // initialize water surface in marsh elements

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( !isFS(el->flag, ELEM::kDry) )
    {
      int ncn = el->Getncn();
      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )  CF( el->nd[i]->flag, NODE::kDry );

      for( int i=0; i<ncn; i++ )
      {
        double Si, Sj, Sk;

        if( isFS(el->nd[i]->flag, NODE::kMarsh) )
        {
          //----------------------------------------------------------------------------------------
          // 01.10.2004, sc
          // The following block of code will initialize the water surface in marsh elements.
          // This will be done only for elements that have been dry in the previous time step.
          //
          // improvement on 16.02.2013, sc
          // 1. A potential water elevation Si at the dry node is computed from adjacent nodes.
          // 2. The water elevation will be initialized according to the variable rewetLimit:
          //    a) Si > Z + rewetLimit                : S = Z + rewetLimit
          //    b) Si < Z + dryLimit                  : Z = S - dryLimit
          //    c) Z + dryLimit < Si < Z + rewetLimit : S = Si

          if( el->nd[i]->mark || firstDryRew )
          {
            int j = (i + 1)       % ncn;
            int k = (i + ncn - 1) % ncn;

            if(     !isFS(el->nd[j]->flag, NODE::kMarsh)
                &&  !isFS(el->nd[k]->flag, NODE::kMarsh) )
            {
              Sj = el->nd[j]->v.S;
              Sk = el->nd[k]->v.S;

              Si = (Sj + Sk) / 2.0;
            }
            else if( !isFS(el->nd[j]->flag, NODE::kMarsh) )
            {
              Si = el->nd[j]->v.S;
            }
            else if( !isFS(el->nd[k]->flag, NODE::kMarsh) )
            {
              Si = el->nd[k]->v.S;
            }
            else
            {
              j = (i + 2) % ncn;

              if( isFS(el->nd[j]->flag, NODE::kMarsh) )
                REPORT::rpt.Error( "unexpected marsh element - dryRewet (1)" );

              Si = el->nd[j]->v.S;
            }

            el->nd[i]->v.S = Si;

            // Adapt the water elevation S or bottom elevation z between dryLimit and rewetLimit.
            if( el->nd[i]->v.S > el->nd[i]->zor + rewetLimit )
            {
              el->nd[i]->v.S = el->nd[i]->zor + rewetLimit;
            }

            if( el->nd[i]->v.S < el->nd[i]->zor + dryLimit )
            {
              el->nd[i]->z = el->nd[i]->v.S - dryLimit;
            }
          }

          else
          {
            // Adapt bottom elevation of dry nodes.
            // changed on 18.04.2006, sc
            if( el->nd[i]->v.S < el->nd[i]->zor + dryLimit )
            {
              el->nd[i]->z = el->nd[i]->v.S - dryLimit;
            }
          }
        }
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // interpolate midside nodes

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    for( int i=ncn; i<nnd; i++ )
    {
      int    left, right;
      double leftS, rightS;
      double leftZ, rightZ;


      // get left and right corner node to midside node i

      el->GetLShape()->getCornerNodes( i, &left, &right );

      leftS  = el->nd[left]->v.S;
      rightS = el->nd[right]->v.S;

      leftZ  = el->nd[left]->z;
      rightZ = el->nd[right]->z;

      el->nd[i]->v.S =  0.5 * (leftS + rightS);
      el->nd[i]->z   =  0.5 * (leftZ + rightZ);

      if(     isFS(el->flag, ELEM::kMarsh)
          &&  el->nd[i]->z + dryLimit/10.0  <  el->nd[i]->zor )
      {
        SF( el->nd[i]->flag, NODE::kMarsh );
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // initialize dry nodes

  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    if( isFS(nd->flag, NODE::kDry) )
    {
      nd->v.U = 0.0;
      nd->v.V = 0.0;
      nd->v.S = nd->z = nd->zor;

      nd->v.K = 0.0;
      nd->v.D = 0.0;

      nd->v.dUdt = 0.0;
      nd->v.dVdt = 0.0;
      nd->v.dSdt = 0.0;
    }
  }


  // -------------------------------------------------------------------------------------
  // report elements

  char text[100];

  REPORT::rpt.Output( "\n (GRID::DryRewet)        the following elements have got dry\n", 5 );

  int cnt = 0;

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( isFS(el->flag, ELEM::kDry)  &&  !el->mark )
    {
      sprintf( text, "  %6d", el->Getname() );
      REPORT::rpt.Output( text, 5 );

      cnt++;
      if( !(cnt % 10) )  REPORT::rpt.Output( "\n", 5 );
    }
  }

  if( cnt % 10 )  REPORT::rpt.Output( "\n", 5 );

  *dried += cnt;


  REPORT::rpt.Output( "\n (GRID::DryRewet)        the following elements have got wet\n", 5 );

  cnt = 0;

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( !isFS(el->flag, ELEM::kDry)  &&  el->mark )
    {
      sprintf( text, "  %6d", el->Getname() );
      REPORT::rpt.Output( text, 5 );

      cnt++;
      if( !(cnt % 10) )  REPORT::rpt.Output( "\n", 5 );
    }
  }

  if( cnt % 10 )  REPORT::rpt.Output( "\n", 5 );

  *wetted += cnt;


  // -------------------------------------------------------------------------------------
  // following the complete list of dry / marsh elements

# ifdef kDebug

  REPORT::rpt.Output( "\n (GRID::DryRewet)        list of dry elements\n", 5 );

  int pos =  0;
  int E1  = -1;
  int E2  = -1;
  int E3  = -1;

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( isFS(el->flag, ELEM::kDry) )
    {
      if( E1 < 0 )
      {
        E1 = el->Getname();
      }

      else if( E2 < 0 )
      {
        E2 = el->Getname();

        if( E2 > E1 + 1 )
        {
          switch( pos )
          {
            case  5:  REPORT::rpt.Output( "\n", 5 );
                      pos = 0;

            default:  sprintf( text, "           %6d", E1 );
                      REPORT::rpt.Output( text, 5 );
                      break;
          }

          pos++;

          E1 = E2;
          E2 = -1;
        }
      }


      else
      {
        E3 = el->Getname();

        if( E3 > E2 + 1 )
        {
          switch( pos )
          {
            case  5:  REPORT::rpt.Output( "\n", 5 );
                      pos = 0;

            default:  sprintf( text, "  %6d - %6d", E1, E2 );
                      REPORT::rpt.Output( text, 5 );
                      break;
          }

          pos++;

          E1 = E3;
          E2 = -1;
        }

        else
        {
          E2 = E3;
        }
      }
    }
  }

  if( E1 > 0  &&  E2 > 0 )
  {
    switch( pos )
    {
      case  5:  REPORT::rpt.Output( "\n", 5 );

      default:  sprintf( text, "  %6d - %6d", E1, E2 );
                REPORT::rpt.Output( text, 5 );
                break;
    }
  }

  else if( E1 > 0 )
  {
    switch( pos )
    {
      case  5:  REPORT::rpt.Output( "\n", 5 );

      default:  sprintf( text, "           %6d", E1 );
                REPORT::rpt.Output( text, 5 );
                break;
    }
  }

  REPORT::rpt.Output( "\n", 5 );

# endif
}


//////////////////////////////////////////////////////////////////////////////////////////
// mark dry nodes and elements (method 3)
// this method will keep wetted nodes wet until countDown counts to zero
//////////////////////////////////////////////////////////////////////////////////////////

void GRID::RewetDry( double  dryLimit,
                     double  rewetLimit,
                     int     countDown,
                     int*    dried,
                     int*    wetted )
{
  // -------------------------------------------------------------------------------------
  // check for dry nodes

  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    // mark recently dry nodes
    if( isFS(nd->flag, NODE::kDry) )  nd->mark = true;
    else                              nd->mark = false;

    CF( nd->flag, NODE::kMarsh );

    double H = nd->v.S - nd->zor;

    if( H < dryLimit )
    {
      SF( nd->flag, NODE::kDry );
      if( nd->countDown > 0 )  nd->countDown--;
    }
    else
    {
      CF( nd->flag, NODE::kDry );
      if( !isFS(nd->flag, NODE::kInlet) )  nd->countDown = countDown;
      else                                  nd->countDown = 0;
    }

    nd->z = nd->zor;
  }


  // -------------------------------------------------------------------------------------
  // loop on all elements: check for dry elements

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    // mark elements that have been dry
    if( isFS(el->flag, ELEM::kDry) )  el->mark = true;
    else                              el->mark = false;

    CF( el->flag, ELEM::kDry | ELEM::kMarsh );

    int ncn = el->Getncn();

    int ndry   = 0;
    int marsh  = false;

    for( int i=0; i<ncn; i++ )
    {
      if( isFS(el->nd[i]->flag, NODE::kDry) )
      {
        marsh = true;
        if( el->nd[i]->countDown == 0 )  ndry++;
      }
    }

    // set dry flag if all nodes are dry
    // set marsh flag if at least one node is dry
    if( ndry == ncn )
    {
      SF( el->flag, ELEM::kDry );
    }

    else if( marsh )
    {
      SF( el->flag, ELEM::kMarsh );

      for( int i=0; i<ncn; i++ )
      {
        if( isFS(el->nd[i]->flag, NODE::kDry) )
        {
          SF( el->nd[i]->flag, NODE::kMarsh );
        }
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // loop on all elements: initialize water surface elevation

  for( int n=0; n<Getnp(); n++ )
  {
    SF( Getnode(n)->flag, NODE::kDry );
  }

  int again = true;

  while( again )
  {
    again = false;

    for( int e=0; e<Getne(); e++ )
    {
      ELEM* el = Getelem(e);

      if( !isFS(el->flag, ELEM::kDry) )
      {
        int ncn = el->Getncn();

        for( int i=0; i<ncn; i++ )
        {
          double Si, Sj, Sk;

          if( isFS(el->nd[i]->flag, NODE::kMarsh) && isFS(el->nd[i]->flag, NODE::kDry) )
          {
            again = true;

            CF( el->nd[i]->flag, NODE::kDry );

            //----------------------------------------------------------------------------------------
            // 01.10.2004, sc
            // The following block of code will initialize the water surface in marsh elements.
            // This will be done only for elements that have been dry in the previous time step.
            //
            // improvement on 16.02.2013, sc
            // 1. A potential water elevation Si at the dry node is computed from adjacent nodes.
            // 2. The water elevation will be initialized according to the variable rewetLimit:
            //    a) Si > Z + rewetLimit                : S = Z + rewetLimit
            //    b) Si < Z + dryLimit                  : Z = S - dryLimit
            //    c) Z + dryLimit < Si < Z + rewetLimit : S = Si

            if( el->nd[i]->mark  ||  firstDryRew )
            {
              int j = (i + 1)       % ncn;
              int k = (i + ncn - 1) % ncn;

              if(     !isFS(el->nd[j]->flag, NODE::kMarsh)
                  &&  !isFS(el->nd[k]->flag, NODE::kMarsh) )
              {
                Sj = el->nd[j]->v.S;
                Sk = el->nd[k]->v.S;

                Si = (Sj + Sk) / 2.0;
              }
              else if( !isFS(el->nd[j]->flag, NODE::kMarsh) )
              {
                Si = el->nd[j]->v.S;
              }
              else if( !isFS(el->nd[k]->flag, NODE::kMarsh) )
              {
                Si = el->nd[k]->v.S;
              }
              else
              {
                continue;
              }

              el->nd[i]->v.S = Si;

              // Adapt the water elevation S or bottom elevation z between dryLimit and rewetLimit.
              if( el->nd[i]->v.S > el->nd[i]->zor + rewetLimit )
              {
                el->nd[i]->v.S = el->nd[i]->zor + rewetLimit;
              }

              if( el->nd[i]->v.S < el->nd[i]->zor + dryLimit )
              {
                el->nd[i]->z = el->nd[i]->v.S - dryLimit;
              }
            }

            else
            {
              // Adapt bottom elevation of dry nodes.
              // changed on 18.04.2006, sc
              if( el->nd[i]->v.S < el->nd[i]->zor + dryLimit )
              {
                el->nd[i]->z = el->nd[i]->v.S - dryLimit;
              }
            }
          }
        }
      }
    }
  }


  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( !isFS(el->flag, ELEM::kDry) )
    {
      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )  CF( el->nd[i]->flag, NODE::kDry );
    }
  }


  // -------------------------------------------------------------------------------------
  // interpolate midside nodes

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    for( int i=ncn; i<nnd; i++ )
    {
      int    left, right;
      double leftS, rightS;
      double leftZ, rightZ;


      // get left and right corner node to midside node i

      el->GetLShape()->getCornerNodes( i, &left, &right );

      leftS  = el->nd[left]->v.S;
      rightS = el->nd[right]->v.S;

      leftZ  = el->nd[left]->z;
      rightZ = el->nd[right]->z;

      el->nd[i]->v.S =  0.5 * (leftS + rightS);
      el->nd[i]->z   =  0.5 * (leftZ + rightZ);

      if(     isFS(el->flag, ELEM::kMarsh)
          &&  el->nd[i]->z + dryLimit/10.0  <  el->nd[i]->zor )
      {
        SF( el->nd[i]->flag, NODE::kMarsh );
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // initialize dry nodes

  for( int n=0; n<Getnp(); n++ )
  {
    NODE* nd = Getnode(n);

    if( isFS(nd->flag, NODE::kDry) )
    {
      nd->v.U = 0.0;
      nd->v.V = 0.0;
      nd->v.S = nd->z = nd->zor;

      nd->v.K = 0.0;
      nd->v.D = 0.0;

      nd->v.dUdt = 0.0;
      nd->v.dVdt = 0.0;
      nd->v.dSdt = 0.0;
    }
  }


  // -------------------------------------------------------------------------------------
  // report elements

  char text[100];

  REPORT::rpt.Output( "\n (GRID::RewetDry)        the following elements have got dry\n", 5 );

  int cnt = 0;

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( isFS(el->flag, ELEM::kDry)  &&  !el->mark )
    {
      sprintf( text, "  %6d", el->Getname() );
      REPORT::rpt.Output( text, 5 );

      cnt++;
      if( !(cnt % 10) )  REPORT::rpt.Output( "\n", 5 );
    }
  }

  if( cnt % 10 )  REPORT::rpt.Output( "\n", 5 );

  *dried += cnt;

  REPORT::rpt.Output( "\n (GRID::RewetDry)        the following elements have got wet\n", 5 );

  cnt = 0;

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( !isFS(el->flag, ELEM::kDry)  &&  el->mark )
    {
      sprintf( text, "  %6d", el->Getname() );
      REPORT::rpt.Output( text, 5 );

      cnt++;
      if( !(cnt % 10) )  REPORT::rpt.Output( "\n", 5 );
    }
  }

  if( cnt % 10 )  REPORT::rpt.Output( "\n", 5 );

  *wetted += cnt;


  // -------------------------------------------------------------------------------------
  // following the complete list of dry / marsh elements

# ifdef kDebug

  REPORT::rpt.Output( "\n (GRID::RewetDry)        list of dry elements\n", 5 );

  int pos =  0;
  int E1  = -1;
  int E2  = -1;
  int E3  = -1;

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    if( isFS(el->flag, ELEM::kDry) )
    {
      if( E1 < 0 )
      {
        E1 = el->Getname();
      }

      else if( E2 < 0 )
      {
        E2 = el->Getname();

        if( E2 > E1 + 1 )
        {
          switch( pos )
          {
            case  5:  REPORT::rpt.Output( "\n", 5 );
                      pos = 0;

            default:  sprintf( text, "           %6d", E1 );
                      REPORT::rpt.Output( text, 5 );
                      break;
          }

          pos++;

          E1 = E2;
          E2 = -1;
        }
      }


      else
      {
        E3 = el->Getname();

        if( E3 > E2 + 1 )
        {
          switch( pos )
          {
            case  5:  REPORT::rpt.Output( "\n", 5 );
                      pos = 0;

            default:  sprintf( text, "  %6d - %6d", E1, E2 );
                      REPORT::rpt.Output( text, 5 );
                      break;
          }

          pos++;

          E1 = E3;
          E2 = -1;
        }

        else
        {
          E2 = E3;
        }
      }
    }
  }

  if( E1 > 0  &&  E2 > 0 )
  {
    switch( pos )
    {
      case  5:  REPORT::rpt.Output( "\n", 5 );

      default:  sprintf( text, "  %6d - %6d", E1, E2 );
                REPORT::rpt.Output( text, 5 );
                break;
    }
  }

  else if( E1 > 0 )
  {
    switch( pos )
    {
      case  5:  REPORT::rpt.Output( "\n", 5 );

      default:  sprintf( text, "           %6d", E1 );
                REPORT::rpt.Output( text, 5 );
                break;
    }
  }

  REPORT::rpt.Output( "\n", 5 );

# endif
}
