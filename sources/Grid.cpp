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
#include "Asciifile.h"
#include "Report.h"
#include "Shape.h"
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Project.h"
#include "Type.h"
#include "Times.h"

#include "Grid.h"


DRYREW GRID::dryRew;


GRID::GRID()
{
  np   = 0;
  node = NULL;

  ne   = 0;
  elem = NULL;

  firstDryRew = true;
}


GRID::~GRID()
{
}


void GRID::KillNode()
{
  delete[] node;
  np = 0;
}


void GRID::KillElem()
{
  delete[] elem;
  ne = 0;
}


int GRID::Getnp()  { return np; }
int GRID::Getne()  { return ne; }

NODE* GRID::Getnode( int n )  { return &node[n]; }
ELEM* GRID::Getelem( int n )  { return &elem[n]; }


void GRID::Setnode( int np, NODE* node )
{
  this->np   = np;
  this->node = node;
}


void GRID::Setelem( int ne, ELEM* elem )
{
  this->ne   = ne;
  this->elem = elem;
}


void GRID::Alloc( int np, int ne )
{
  this->np = np;

  if( np )
  {
    node = new NODE [np];
    if( !node )  REPORT::rpt.Error( kMemoryFault, "can not allocate memory (GRID::Alloc - 1)" );

    for( int i=0; i<this->np; i++ )  node[i].Setno(i);
  }

  this->ne = ne;

  if( ne )
  {
    elem = new ELEM [ne];
    if( !elem )  REPORT::rpt.Error( kMemoryFault, "can not allocate memory (GRID::Alloc - 2)" );

    for( int i=0; i<this->ne; i++ )  elem[i].Setno(i);
  }
}


void GRID::Free()
{
  if( np )  delete[] node;
  np = 0;

  if( ne )  delete[] elem;
  ne = 0;
}


void GRID::InputRegion( char* fileName, SUBDOM* subdom )
{
  char  text[200];
  char* textLine;

  ////////////////////////////////////////////////////////////////////////////////////////
  // open the region file

  ASCIIFILE* file = new ASCIIFILE( fileName, "r" );

  if( !file || !file->getid() )
    REPORT::rpt.Error( kOpenFileFault, "%s %s (GRID::InputRegion - 1)",
             "can not open region file", fileName );


  ////////////////////////////////////////////////////////////////////////////////////////
  // read number of nodes and number of elements

  int npreg = 0;
  int nereg = 0;
  int ndata = 0;

  textLine = file->nextLine();
  sscanf( textLine, "%d %d %d", &npreg, &nereg, &ndata );

  if( subdom->npr > 1 )
  {
    if( subdom->np != npreg )
      REPORT::rpt.Error( kParameterFault, "%s (GRID::InputRegion - 2)",
               "different number of nodes in subdomain and region file" );

    if( subdom->ne != nereg )
      REPORT::rpt.Error( kParameterFault, "%s (GRID::InputRegion - 3)",
               "different number of elements in subdomain and region file" );
  }

  else
  {
    // allocate memory for nodes and elements --------------------------------------------
    Alloc( npreg, nereg );
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // read nodes

  for( int i=0; i<npreg; i++ )
  {
    // read node name and coordinates
    int    name;
    double x, y, z;

    textLine = file->nextLine();
    sscanf( textLine, "%d %lf %lf %lf", &name, &x, &y, &z );

    if( name <= 0  ||  name > npreg )
    {
      sprintf( text, "%s %d - (GRID::InputRegion - 4)",
               "node numbers must be between 1 and", npreg );
      REPORT::rpt.Error( text );
    }

    // initializations
    NODE* nd = NULL;

    if( subdom->npr > 1 )
    {
      nd = subdom->node[name - 1];
    }
    else
    {
      nd = &this->node[name - 1];
    }

    if( nd )
    {
      nd->Setname( name );

      nd->x     = x;
      nd->y     = y;
      nd->z     = z;
      nd->zor   = z;
      nd->zero  = z - 1000.0;
      nd->noel  = 0;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // read elements

  for( int i=0; i<nereg; i++ )
  {
    int  mat;
    int  ncn, nnd;
    int  name, n[kMaxNodes2D];
    char shape[10];

    textLine = file->nextLine();
    sscanf( textLine, "%d %d %s ", &name, &mat, shape );

    // read only elements of type "tri" and "quad"
    if( strcmp(shape, "tri")  &&  strcmp(shape, "quad") )  continue;

    if( name <= 0  ||  name > nereg )
    {
      sprintf( text, "%s %d (GRID::InputRegion - 5)",
               "element numbers must be between 1 and", nereg );
      REPORT::rpt.Error( text );
    }

    ELEM* el = NULL;

    if( subdom->npr > 1 )
    {
      el = subdom->elem[name - 1];
    }
    else
    {
      el = &this->elem[name - 1];
    }

    if( el )
    {
      el->Setname( name );

      TYPE* elt = TYPE::Getno( mat );
      if( elt->no(0) > 0 )  el->type = elt->id();
      else                  el->type = -mat;

      ScanElement( el, textLine, shape, n, &ncn, &nnd );


      // set up pointers to nodes --------------------------------------------------------

      for( int j=0; j<nnd; j++ )
      {
        if( n[j] <= 0  ||  n[j] > npreg )
        {
          REPORT::rpt.Message( 0, "wrong node number %d in element %d %s\n",
                                  n[j], el->Getname(), "(GRID::InputRegion - 6)" );
        }

        else
        {
          if( subdom->npr > 1 )
          {
            el->nd[j] = subdom->node[n[j] - 1];
            if( !el->nd[j] )
            {
              REPORT::rpt.Error( kUnexpectedFault, "node %d of element %d not in domain",
                                                   n[j], el->Getname() );
            }
          }
          else
          {
            el->nd[j] = &this->node[n[j] - 1];
          }
        }
      }
    }
  }

  if( ndata > 0 )
  {
    // read labels and units -------------------------------------------------------------
    // ...

    int ncomp = 0;       // number of components

    textLine = file->nextLine();
    sscanf( textLine, "%d", &ncomp );

    for( int i=0; i<ncomp; i++ )
    {
      file->nextLine();
    }

    for( int i=0; i<npreg; i++ )
    {
      int    name = 0;
      double zero = 0.0;

      textLine = file->nextLine();
      sscanf( textLine, "%d %lf", &name, &zero );

      if( name <= 0  ||  name > npreg )
      {
        REPORT::rpt.Message( 0, "wrong node number %d in data %s\n",
                                name, "(GRID::InputRegion - 7)" );
      }


      NODE* nd = NULL;

      if( subdom->npr > 1 )
      {
        nd = subdom->node[name - 1];
      }
      else
      {
        nd = &this->node[name - 1];
      }

      if( nd )
      {
        nd->zero = zero;
      }
    }
  }

  delete file;


  // interpolate bottom elevation at midside nodes (linear interpolation) ----------------

  for( int e=0; e<Getne(); e++ )
  {
    ELEM* el = Getelem(e);

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    for( int i=ncn; i<nnd; i++ )
    {
      // get left and right corner node to midside node i ------------------------------

      int    il, ir;
      double left, rght;

      el->GetQShape()->getCornerNodes( i, &il, &ir );

      left = el->nd[il]->zor;
      rght = el->nd[ir]->zor;
      el->nd[i]->z = el->nd[i]->zor = 0.5 * (left + rght);
    }
  }
}



void GRID::InputControl( char* fileName, GRID* region, SUBDOM* subdom )
{
  char  text[200];
  char* textLine;

  ////////////////////////////////////////////////////////////////////////////////////////
  // open the region file

  ASCIIFILE* file = new ASCIIFILE( fileName, "r" );

  if( !file || !file->getid() )
    REPORT::rpt.Error( kOpenFileFault, "%s %s (GRID::InputControl - 1)",
             "can not open control file", fileName );


  ////////////////////////////////////////////////////////////////////////////////////////
  // read number of nodes and number of elements

  int npctr, nectr;

  textLine = file->nextLine();
  sscanf( textLine, "%d %d\n", &npctr, &nectr );


  ////////////////////////////////////////////////////////////////////////////////////////
  // not interested in nodes

  for( int i=0; i<npctr; i++ )  textLine = file->nextLine();


  ////////////////////////////////////////////////////////////////////////////////////////
  // read elements

  int* con[5];
  for( int i=0; i<5; i++ )
    con[i] = (int*) MEMORY::memo.Array_el( nectr );

  ne = 0;

  for( int i=0; i<nectr; i++ )
  {
    int  mat;
    int  name, n[3];
    char shape[10];

    textLine = file->nextLine();
    sscanf( textLine, "%d %d %s %d %d %d", &name, &mat, shape, &n[0], &n[1], &n[2] );

    // read only elements of type "line"
    if( strcmp(shape, "line") )  continue;

//    if( name <= 0  ||  name > nectr )
//    {
//      sprintf( text, "%s %d (GRID::InputControl - 2)\n",
//               "element numbers should be between 1 and", nectr );
//      REPORT::rpt.Message( 0, text );
//    }

    for( int j=0; j<3; j++ )
    {
      if( n[j] <= 0  ||  (subdom->npr == 1  &&  n[j] > region->Getnp())
                     ||  (subdom->npr  > 1  &&  n[j] > subdom->np) )
      {
        REPORT::rpt.Message( 0, "wrong node number %d in element %d %s\n",
                                n[j], name, "(GRID::InputControl - 3)" );
      }
    }

    if( subdom->npr > 1  &&
        (    !subdom->node[n[0] - 1]
          || !subdom->node[n[1] - 1]
          || !subdom->node[n[2] - 1] )  )  continue;

    con[0][ne] = n[0] - 1;
    con[1][ne] = n[1] - 1;
    con[2][ne] = n[2] - 1;

    con[3][ne] = name;
    con[4][ne] = mat;

    ne++;
  }


  // allocate memory for elements --------------------------------------------------------

  Alloc( 0, ne );

  for( int i=0; i<ne; i++ )
  {
    ELEM* el = &this->elem[i];

    el->Setshape( kLine );

    el->Setname( con[3][i] );
    el->type = con[4][i];


    // set up pointers to nodes ----------------------------------------------------------

    if( subdom->npr > 1 )
    {
      el->nd[0] = subdom->node[con[0][i]];
      el->nd[1] = subdom->node[con[1][i]];
      el->nd[2] = subdom->node[con[2][i]];
    }

    else
    {
      el->nd[0] = &region->node[con[0][i]];
      el->nd[1] = &region->node[con[1][i]];
      el->nd[2] = &region->node[con[2][i]];
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  // free allocated memory

  delete file;

  for( int i=0; i<5; i++ )  MEMORY::memo.Delete( con[i] );


  ////////////////////////////////////////////////////////////////////////////////////////
  // set number of nodes and a pointer to nodes in the region grid

  this->np   = region->np;
  this->node = region->node;
}


void GRID::InputInitial( int isAscii, char* name, TIME* time, SUBDOM* subdom, int zb_init )
{
  for( int n=0; n<np; n++ )
  {
    NODE* nd = &node[n];

    nd->v.U    = 0.0;
    nd->v.V    = 0.0;
    nd->v.K    = 0.0;
    nd->v.D    = 0.0;
    nd->v.C    = 0.0;
    nd->v.Qb   = 0.0;

    nd->v.dUdt = 0.0;
    nd->v.dVdt = 0.0;
    nd->v.dSdt = 0.0;

    nd->v.dKdt = 0.0;
    nd->v.dDdt = 0.0;

    nd->v.S    = nd->z;

    nd->qbo    = 0.0;
    nd->dz     = 0.0;
  }

  for( int e=0; e<ne; e++ )
  {
    ELEM* el = &elem[e];

    el->U  = 0.0;
    el->V  = 0.0;
    el->P  = 0.0;
    el->dz = 0.0;

    int ncn = el->Getncn();

    for( int i=0; i<ncn; i++ )  el->P += el->nd[i]->z;
    el->P /= ncn;
  }


  if( name[0] )
  {
    int npInit = 0;
    int neInit = 0;

    if( isAscii )
    {
      char* textLine;
      int   release = 0;

      ASCIIFILE* file = new ASCIIFILE( name, "r" );

      if( !file || !file->getid() )
        REPORT::rpt.Error( kOpenFileFault, "%s %s (GRID::InputInitial - 1)",
                           "can not open initial file", name );

      textLine = file->next();

      const char*  vars[] = { "U","V","W","S","dUdt","dVdt","dSdt","K","D","C","qb","Zb", "\0" };
      enum                  { kU, kV, kW, kS, kDUDT, kDVDT, kDSDT, kK, kD, kC, kQB, kZB  };
      double vals[]       = { 0.0,0.0,0.0,0.0,  0.0,   0.0,   0.0, 0.0,0.0,0.0, 0.0, 0.0 };

#     define kMaxtok  100
      int    present[kMaxtok];
      for( int i=0; i<kMaxtok; i++ )  present[i] = -1;

      char  list[500];
      char  seps[] = " ,\t\n\r";
      char* token;

      strcpy( list, textLine+1 );
      token = strtok( list, seps );

      for( int ntok=0; token!=NULL; ntok++ )
      {
        if( ntok >= kMaxtok )  break;

        if( ntok == 0 )
        {
          time->Set( token );                               // read time
        }
        else if( ntok == 2 )
        {
          sscanf( token, "%d", &release );                  // read release of data set
        }
        else if( ntok > 2  &&  release >= 40000 )
        {
          for( int i=0; vars[i][0]; i++ )
          {
            if( strcmp(token, vars[i]) == 0 )  present[ntok-3] = i;
          }
        }

        token = strtok( NULL, seps );
      }

      textLine = file->nextLine();
      sscanf( textLine, "%d %d", &npInit, &neInit );

      if(     (subdom->npr == 1  &&  npInit != np)
          ||  (subdom->npr  > 1  &&  npInit != subdom->np) )
        REPORT::rpt.Error( "wrong number of nodes in file (GRID::InputInitial - 2)" );

      for( int n=0; n<npInit; n++ )
      {
        // read name of node -------------------------------------------------------------
        int name;

        textLine = file->nextLine();
        sscanf( textLine, "%d", &name );

        NODE* nd = NULL;

        if( subdom->npr > 1 )  nd = subdom->node[name - 1];
        else                   nd = &node[name - 1];

        // initialize nodal values -------------------------------------------------------
        for( int i=0; vars[i][0]; i++ )  vals[i] = 0.0;

        if( nd )
        {
          vals[kZB] = nd->zor;
          vals[kS]  = nd->v.S;
        }

        // read node data ----------------------------------------------------------------
        if( release >= 40000 )
        {
          char  list[200];
          char  seps[] = " ,\t\n\r";
          char* token;

          strcpy( list, textLine );
          token = strtok( list, seps );

          for( int ntok=0; token!=NULL; ntok++ )
          {
            if( ntok == 0 )
            {
              sscanf( token, "%d", &name );
            }
            else if( present[ntok-1] >= 0 )
            {
              sscanf( token, "%lf", &vals[present[ntok-1]] );
            }

            token = strtok( NULL, seps );
          }
        }

        else if( release >= 280 )
        {
          double dummy;
          sscanf( textLine,
                  "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &name,
                  &vals[kU], &vals[kV], &vals[kS], &vals[kK], &vals[kD], &vals[kC],
                  &vals[kDUDT], &vals[kDVDT], &vals[kDSDT], &dummy );
        }

        else
        {
          double dummy;
          sscanf( textLine,
                  "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &name,
                  &vals[kU], &vals[kV], &vals[kW], &vals[kS], &vals[kK], &vals[kD], &vals[kC],
                  &vals[kDUDT], &vals[kDVDT], &vals[kDSDT], &dummy );
        }

        if( nd )
        {
          nd->v.U    = vals[kU];
          nd->v.V    = vals[kV];
          nd->v.S    = vals[kS];

          nd->v.dUdt = vals[kDUDT];
          nd->v.dVdt = vals[kDVDT];
          nd->v.dSdt = vals[kDSDT];

          nd->v.K    = vals[kK];
          nd->v.D    = vals[kD];

          nd->v.C    = vals[kC];

          nd->v.Qb   = vals[kQB];
          nd->qbo    = nd->v.Qb;

          nd->dz     = vals[kZB] - nd->zor;

          if( zb_init )
          {
            nd->z    =
            nd->zor  = vals[kZB];
          }
        }
      }

      for( int e=0; e<ne; e++ )
      {
        ELEM* el = el = &elem[e];

        el->U  = 0.0;
        el->V  = 0.0;
        el->P  = 0.0;
        el->dz = 0.0;

        int ncn = el->Getncn();

        for( int i=0; i<ncn; i++ )
        {
          el->U += el->nd[i]->v.U;
          el->V += el->nd[i]->v.V;
          el->P += el->nd[i]->v.S;
        }

        el->U /= ncn;
        el->V /= ncn;
        el->P /= ncn;
      }

      if( release >= 320  &&  neInit > 0 )
      {
        if(     (subdom->npr == 1  &&  neInit != ne)
            ||  (subdom->npr  > 1  &&  neInit != subdom->ne) )
          REPORT::rpt.Error( "wrong number of elements in file (GRID::InputInitial - 3)" );

        for( int e=0; e<neInit; e++ )
        {
          int   name;
          ELEM* el = NULL;

          textLine = file->nextLine();
          sscanf( textLine, "%d", &name );

          if( subdom->npr > 1 )  el = subdom->elem[name - 1];
          else                   el = &elem[name - 1];

          double U  = 0.0;
          double V  = 0.0;
          double P  = 0.0;
          double dz = 0.0;

          if( el )
          {
            U  = el->U;
            V  = el->V;
            P  = el->P;
            dz = el->dz;
          }

          if( release >= 40000 )  sscanf( textLine, "%d %lf %lf %lf %lf", &name, &U, &V, &P, &dz );
          else                    sscanf( textLine, "%d %lf", &name, &P );

          if( el )
          {
            el->U  = U;
            el->V  = V;
            el->P  = P;
            el->dz = dz;
          }
        }
      }

      delete file;
    }


    else
    {
      char   release[8];
      int    relno;
      size_t len;

      FILE* id = fopen( name, "rb" );

      if( !id )  REPORT::rpt.Error( "can not open initial file (GRID::InputInitial - 4)" );


      ////////////////////////////////////////////////////////////////////////////////////
      // read file header

      size_t ret;

      len = 7;
      ret = fread( release, sizeof(char), len, id );

      if( strcmp(release, "Release") == 0 )
      {
        len = 1;
        ret = fread( &relno, sizeof(int), len, id );

        if( relno >= 40000 )
        {
          char stime[23];
          len = 22;
          ret = fread( stime, sizeof(char), len, id );
          time->Set( stime );
        }
        else
        {
          double dtime;
          len = 1;
          ret = fread( &dtime, sizeof(double), len, id );
          time->Setsec( dtime );
        }

        len = 1;
        ret = fread( &npInit, sizeof(int), len, id );
        ret = fread( &neInit, sizeof(int), len, id );

        if(     (subdom->npr == 1  &&  npInit != np)
            ||  (subdom->npr  > 1  &&  npInit != subdom->np) )
          REPORT::rpt.Error( "wrong number of nodes in file (GRID::InputInitial - 5)" );

        if(     (subdom->npr == 1  &&  neInit != ne)
            ||  (subdom->npr  > 1  &&  neInit != subdom->ne) )
          REPORT::rpt.Error( "wrong number of elements in file (GRID::InputInitial - 6)" );
      }
      else
      {
        rewind( id );

        double dtime;
        len = 1;
        ret = fread( &dtime,   sizeof(double), len, id );
        time->Setsec( dtime );

        len = 1;
        ret = fread( &npInit, sizeof(int),    len, id );

        if(     (subdom->npr == 1  &&  npInit != np)
            ||  (subdom->npr  > 1  &&  npInit != subdom->np) )
          REPORT::rpt.Error( "wrong number of nodes in file (GRID::InputInitial - 7)" );
      }


      ////////////////////////////////////////////////////////////////////////////////////
      // read node data

      double* nd_data = (double*) MEMORY::memo.Array_nd( npInit );

      for( int i=0; i<11; i++ )
      {
        len = npInit;
        ret = fread( nd_data, sizeof(double), len, id );

        switch( i )
        {
          case  0:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.U = nd_data[n];
            }
            break;

          case  1:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.V = nd_data[n];
            }
            break;

          case  2:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.S = nd_data[n];
            }
            break;

          case  3:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.K = nd_data[n];
            }
            break;

          case  4:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.D = nd_data[n];
            }
            break;

          case  5:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.C = nd_data[n];
            }
            break;

          case  6:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.dUdt = nd_data[n];
            }
            break;

          case  7:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.dVdt = nd_data[n];
            }
            break;

          case  8:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )  nd->v.dSdt = nd_data[n];
            }
            break;

          case  9:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )
              {
                nd->v.Qb = nd_data[n];
                nd->qbo  = nd->v.Qb;
              }
            }
            break;

          case 10:
            for( int n=0; n<npInit; n++ )
            {
              NODE* nd = NULL;
              if( subdom->npr > 1 )  nd = subdom->node[n];
              else                   nd = &node[n];
              if( nd )
              {
                nd->dz     = nd_data[n] - nd->zor;
                nd->z      =
                nd->zor    = nd_data[n];
              }
            }
            break;
        }
      }

      MEMORY::memo.Detach( nd_data );


      ////////////////////////////////////////////////////////////////////////////////////
      // read element data

      for( int e=0; e<ne; e++ )
      {
        ELEM* el = &elem[e];

        el->U = 0.0;
        el->V = 0.0;
        el->P = 0.0;

        int ncn = el->Getncn();

        for( int i=0; i<ncn; i++ )
        {
          el->U += el->nd[i]->v.U;
          el->V += el->nd[i]->v.V;
          el->P += el->nd[i]->v.S;
        }

        el->U /= ncn;
        el->V /= ncn;
        el->P /= ncn;
      }

      double* ne_data = (double*) MEMORY::memo.Array_el( neInit );

      if( relno >= 40000  &&  neInit > 0  )
      {
        for( int i=0; i<4; i++ )
        {
          len = neInit;
          ret = fread( ne_data, sizeof(double), len, id );

          switch( i )
          {
            case 0:
              for( int e=0; e<neInit; e++ )
              {
                ELEM* el = NULL;
                if( subdom->npr > 1 )  el = subdom->elem[e];
                else                   el = &elem[e];
                if( el )  el->U = ne_data[e];
              }
              break;

            case 1:
              for( int e=0; e<neInit; e++ )
              {
                ELEM* el = NULL;
                if( subdom->npr > 1 )  el = subdom->elem[e];
                else                   el = &elem[e];
                if( el )  el->V = ne_data[e];
              }
              break;

            case 2:
              for( int e=0; e<neInit; e++ )
              {
                ELEM* el = NULL;
                if( subdom->npr > 1 )  el = subdom->elem[e];
                else                   el = &elem[e];
                if( el )  el->P = ne_data[e];
              }
              break;

            case 3:
              for( int e=0; e<neInit; e++ )
              {
                ELEM* el = NULL;
                if( subdom->npr > 1 )  el = subdom->elem[e];
                else                   el = &elem[e];
                if( el )  el->dz = ne_data[e];
              }
              break;
          }
        }
      }

      else if( relno >= 320  &&  neInit > 0  )
      {
        len = neInit;
        ret = fread( ne_data, sizeof(double), len, id );

        for( int e=0; e<neInit; e++ )
        {
          ELEM* el = NULL;

          if( subdom->npr > 1 )  el = subdom->elem[e];
          else                   el = &elem[e];

          if( el )
          {
            el->P = ne_data[e];
          }
        }
      }

      MEMORY::memo.Detach( ne_data );

      fclose( id );
    }


    for( int e=0; e<ne; e++ )
    {
      ELEM* el = &elem[e];

      int ncn = el->Getncn();

      double Z  = 0.0;
      double S = 0.0;

      for( int i=0; i<ncn; i++ )
      {
        Z += el->nd[i]->z;
        S += el->nd[i]->v.S;
      }

      Z /= ncn;
      S /= ncn;

      if( el->P < Z )  el->P = S;
    }
  }
}


void GRID::InitPrevious()
{
  for( int n=0; n<np; n++ )
  {
    NODE* nd = &node[n];

    nd->vo.U = nd->v.U;
    nd->vo.V = nd->v.V;
    nd->vo.S = nd->v.S;

    nd->vo.K = nd->v.K;
    nd->vo.D = nd->v.D;
    nd->vo.C = nd->v.C;

    nd->vo.dUdt = nd->v.dUdt;
    nd->vo.dVdt = nd->v.dVdt;
    nd->vo.dSdt = nd->v.dSdt;

    nd->vo.Qb   = nd->v.Qb;

    nd->vo.dKdt = 0.0;
    nd->vo.dDdt = 0.0;
  }
}


void GRID::ScanElement( ELEM* el,
                        char* textLine,
                        char* elemShape,
                        int*  n,
                        int*  ncn,
                        int*  nnd )
{
  int no, mat;

  if( strcmp(elemShape, "line") == 0 )
  {
    sscanf( textLine, "%d %d %s %d %d %d",
            &no, &mat, elemShape,
            &n[0], &n[1], &n[2] );
    *ncn = 2;
    *nnd = 3;

    el->Setshape( kLine );
  }

  else if( strcmp(elemShape, "tri") == 0 )
  {
    sscanf( textLine, "%d %d %s %d %d %d %d %d %d",
            &no, &mat, elemShape,
            &n[0], &n[1], &n[2], &n[3], &n[4], &n[5] );
    *ncn = 3;
    *nnd = 6;

    el->Setshape( kTriangle );
  }

  else if( strcmp(elemShape, "quad") == 0 )
  {
    sscanf( textLine, "%d %d %s %d %d %d %d %d %d %d %d",
            &no, &mat, elemShape,
            &n[0], &n[1], &n[2], &n[3], &n[4], &n[5], &n[6], &n[7] );
    *ncn = 4;
    *nnd = 8;

    el->Setshape( kSquare );
  }
}


void GRID::OutputData( PROJECT* project, int timeStep, char* time )
{
  if( !project->name.restartFile[0] )  return;

  char fileName[500];
  project->ReplaceMacro( project->name.restartFile, fileName, timeStep, 0 );

  // -------------------------------------------------------------------------------------
  // write data file
  if( project->name.ascii_restart )
  {
    FILE* id = fopen( fileName, "w" );

    if( !id )
    {
      REPORT::rpt.Error( kOpenFileFault,
                         "can not open restart file (GRID::outputData - 1)" );
      return;
    }

    fprintf( id, "# %22s   Release %d   U,V,S,dUdt,dVdt,dSdt,K,D,C,qb,Zb\n", time, (int)(project->release) );
    fprintf( id, "%d  %d\n", np, ne );

    // write nodal values ----------------------------------------------------------------
    for( int n=0; n<np; n++ )
    {
      double U    = node[n].v.U;
      double V    = node[n].v.V;
      double S    = node[n].v.S;
      double dUdt = node[n].v.dUdt;
      double dVdt = node[n].v.dVdt;
      double dSdt = node[n].v.dSdt;
      double K    = node[n].v.K;
      double D    = node[n].v.D;
      double C    = node[n].v.C;
      double Qb   = node[n].v.Qb;
      double Zb   = node[n].zor;

      fprintf( id, "%7d", node[n].Getname() );
      fprintf( id, " %14.6le %14.6le %14.6le", U, V, S );
      fprintf( id, " %14.6le %14.6le %14.6le", dUdt, dVdt, dSdt );
      fprintf( id, " %14.6le %14.6le %14.6le", K, D, C );
      fprintf( id, " %14.6le %14.6le\n", Qb, Zb );
    }

    // write element values --------------------------------------------------------------
    for( int e=0; e<ne; e++ )
    {
      double U  = elem[e].U;
      double V  = elem[e].V;
      double P  = elem[e].P;
      double dz = elem[e].dz;

      fprintf( id, "%7d", elem[e].Getname() );
      fprintf( id, " %14.6le %14.6le %14.6le %14.6le\n", U, V, P, dz );
    }

    fclose( id );
  }

  else
  {
    FILE* id = fopen( fileName, "wb" );

    if( !id )
    {
      REPORT::rpt.Error( kOpenFileFault,
                         "can not open restart file (GRID::outputData - 2)", 0 );
      return;
    }

    size_t len = 7;
    fwrite( "Release",                          sizeof(char),  len, id );

    len = 1;
    int rel = (int)(project->release);
    fwrite( &rel,                               sizeof(int),   len, id );

    len = 22;
    fwrite( time,                               sizeof(char),  len, id );

    len = 32;
    fwrite( "U,V,S,dUdt,dVdt,dSdt,K,D,C,qb,Zb", sizeof(char),  len, id );

    len = 1;
    fwrite( &np,                                sizeof(int),   len, id );
    fwrite( &ne,                                sizeof(int),   len, id );

    double* data;
    if( np > ne )  data = (double*) MEMORY::memo.Array_nd( np );
    else           data = (double*) MEMORY::memo.Array_el( ne );

    len = np;

    for( int n=0; n<np; n++ )  data[n] = node[n].v.U;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.V;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.S;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.K;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.D;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.C;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.dUdt;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.dVdt;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.dSdt;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].v.Qb;
    fwrite( data, sizeof(double), len, id );

    for( int n=0; n<np; n++ )  data[n] = node[n].zor;
    fwrite( data, sizeof(double), len, id );

    len = ne;

    for( int e=0; e<ne; e++ )  data[e] = elem[e].U;
    fwrite( data, sizeof(double), len, id );

    for( int e=0; e<ne; e++ )  data[e] = elem[e].V;
    fwrite( data, sizeof(double), len, id );

    for( int e=0; e<ne; e++ )  data[e] = elem[e].P;
    fwrite( data, sizeof(double), len, id );

    for( int e=0; e<ne; e++ )  data[e] = elem[e].dz;
    fwrite( data, sizeof(double), len, id );

    MEMORY::memo.Detach( data );

    fclose( id );
  }
}


void GRID::OutputGeom( PROJECT* project, int timeStep )
{
  char fileName[500];
  project->ReplaceMacro( project->name.geometryFile, fileName, timeStep, 0 );

  // -------------------------------------------------------------------------------------
  // write region file
  FILE* id = fopen( fileName, "w" );
  if( !id )
    REPORT::rpt.Error( "can not open geometry output file (GRID::outputGeom - 1)" );

  // write number of nodes and elements --------------------------------------------------
  fprintf( id, "%6d  %6d  0  0  0\n", np, ne );

  // write nodes -------------------------------------------------------------------------
  for( int n=0; n<np; n++ )
  {
    double X = node[n].x;
    double Y = node[n].y;
    double Z = node[n].zor;

    fprintf( id, "%6d  %12.6lf  %12.6lf  %12.6lf\n", node[n].Getname(), X, Y, Z );
  }

  // write element connectivity ----------------------------------------------------------
  for( int e=0; e<ne; e++ )
  {
    ELEM* el = &elem[e];

    char elemShape[6];

    if( isFS(el->flag, ELEM::kRegion) )
    {
      switch( el->GetshapeID() )
      {
        //case kLine:      strcpy (elemShape, "line");   break;
        case kTriangle:  strcpy (elemShape, "tri");    break;
        case kSquare:    strcpy (elemShape, "quad");   break;
      }

      int mat = TYPE::Getid(el->type)->no(el->type);

      fprintf( id, "%6d  %5d  %5s", el->Getname(), mat, elemShape );

      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )  fprintf( id, "  %6d", el->nd[i]->Getname() );

      fprintf( id, "\n" );
    }
  }

  fclose( id );
}
