// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class SUBDOM
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
#include "Memory.h"
#include "Node.h"
#include "Elem.h"
#include "Grid.h"
#include "Model.h"
#include "Eqs.h"
#include "Project.h"

#include "Subdom.h"


//#define kDebug


//////////////////////////////////////////////////////////////////////////////////////////
// global variables


//////////////////////////////////////////////////////////////////////////////////////////
// constructor & destructor

SUBDOM::SUBDOM()
{
  pid = 0;
  npr = 1;

  inface = NULL;
  subbuf = NULL;
}

SUBDOM::~SUBDOM()
{
  if( inface )  delete[] inface;
  if( subbuf )  delete[] subbuf;
}

void SUBDOM::Input( char* subdomFileName, char* regionFileName, GRID* region )
{
  char  text[200];
  char* textLine;

  // -------------------------------------------------------------------------------------
  // read region file and allocate memory for arrays SD_NODE sdnd[] and SD_ELEM sdel[]

  ASCIIFILE* regionFile = new ASCIIFILE( regionFileName, "r" );
  if( !regionFile || !regionFile->getid() )
    REPORT::rpt.Error( kOpenFileFault, "%s %s (SUBDOM::Input - 1)",
                       "can not open region file", regionFileName );

  textLine = regionFile->nextLine();
  sscanf( textLine, " %d %d", &np, &ne );

  sdnd = new SD_NODE[np];
  sdel = new SD_ELEM[ne];

  node = new NODE* [np];
  elem = new ELEM* [ne];

  if( !sdnd || !sdel || !node || !elem )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (SUBDOM::Input - 2)" );

  // temp array to hold element-node-connectivity
  int* con[kMaxNodes2D+1];
  for( int i=0; i<=kMaxNodes2D; i++ )
    con[i] = (int*) MEMORY::memo.Array_el( ne );

  // not interested in node co-ordinates
  for( int i=0; i<np; i++ )
    regionFile->nextLine();

  // read element connectivity
  for( int i=0; i<ne; i++ )
  {
    int  name;
    int  mat;
    char shape[20];

    regionFile->nextLine();
    sscanf( textLine, "%d %d %s", &name, &mat, shape );

    if( strcmp(shape, "tri") == 0 )
    {
      int n[6];
      sscanf( textLine, "%d %d %s %d %d %d %d %d %d",
              &name, &mat, shape,
              &n[0], &n[1], &n[2], &n[3], &n[4], &n[5] );

      int no = name - 1;

      con[0][no] = 6;
      con[1][no] = n[0] - 1;    con[4][no] = n[3] - 1;
      con[2][no] = n[1] - 1;    con[5][no] = n[4] - 1;
      con[3][no] = n[2] - 1;    con[6][no] = n[5] - 1;
    }
    else if( strcmp(shape, "quad") == 0 )
    {
      int n[8];
      sscanf( textLine, "%d %d %s %d %d %d %d %d %d %d %d",
              &name, &mat, shape,
              &n[0], &n[1], &n[2], &n[3], &n[4], &n[5], &n[6], &n[7] );

      int no = name - 1;

      con[0][no] = 8;
      con[1][no] = n[0] - 1;    con[5][no] = n[4] - 1;
      con[2][no] = n[1] - 1;    con[6][no] = n[5] - 1;
      con[3][no] = n[2] - 1;    con[7][no] = n[6] - 1;
      con[4][no] = n[3] - 1;    con[8][no] = n[7] - 1;
    }
  }

  delete regionFile;


  // -------------------------------------------------------------------------------------
  // read subdomain file and set elements

  ASCIIFILE* subdomFile = new ASCIIFILE( subdomFileName, "r" );
  if( !subdomFile || !subdomFile->getid() )
    REPORT::rpt.Error( kOpenFileFault, "%s %s (SUBDOM::Input - 2)",
                       "can not open subdomain file", subdomFileName );

  int sub_ne = 0;
  textLine = subdomFile->nextLine();
  sscanf( textLine, " %d", &sub_ne );

  if( sub_ne != ne )
    REPORT::rpt.Error( kOpenFileFault,
             "wrong number of elements in subdomain file (SUBDOM::Input - 2)" );


  int nsub = 0;

  for( int e=0; e<ne; e++ )
  {
    int name;
    int sub;

    textLine = subdomFile->nextLine();
    sscanf( textLine, " %d %d", &name, &sub );

    int no = name - 1;
    sdel[no].sub = sub - 1;

    if( sub > nsub ) nsub = sub;
  }

  delete subdomFile;


  // -------------------------------------------------------------------------------------
  // output info on domain decomposition

  if( nsub != npr )
    REPORT::rpt.Error( kParameterFault,
             "differing number of processes and subdomains (SUBDOM::Input - 3)" );

  sprintf( text, "\n (SUBDOM::Input)         %s %d\n",
                 "performing domain decomposition:", nsub );
  REPORT::rpt.Output( text, 3 );


  // -------------------------------------------------------------------------------------
  // definition of interface nodes in subdomains
  //
  //            O--o--O--o--O--o--O--o--O--o--O--o--O
  //            |     |     |    pid    |     |     |
  //            o     o     o     4     o     o     o
  //            |     |     |     |     |     |     |
  //  ========= H==h==H==h==H==h==H==h==H==h==H==h==H === interface (3-4)
  //            |     |     |     |     |     |     |
  //            o     o     o     o     o     o     o
  //            |     |     |           |     |     |
  //            O--o--O--o--O--  pid  --O--o--O--o--O
  //            |     |     |     3     |     |     |
  //            o     o     o           o     o     o
  //            |     |     |     |     |     |     |
  //  ========= Q==q==Q==q==Q==q==X==q==Q==q==Q==q==Q === interface (1-3), (2-3)
  //            |     |     |     +     |     |     |
  //            o     o pid o     o     o pid o     o     X = crosspoint
  //            |     |  1  |     +     |  2  |     |
  //            O--o--O--o--O--o--O--o--O--o--O--o--O
  //
  // upstream interface nodes (interface 3-1, 3-2)
  // downstream interface nodes (interface 3-4)

  // initialisation ----------------------------------------------------------------------

  for( int n=0; n<np; n++ )
  {
    sdnd[n].sub  = -1;
    sdnd[n].mark = false;

    CF( sdnd[n].flag, SD_NODE::kInface );      // interface
    CF( sdnd[n].flag, SD_NODE::kInface_UP );   // upstream interface
    CF( sdnd[n].flag, SD_NODE::kInface_DN );   // downstream interface
  }


  // allocate memory for interface class INFACE ------------------------------------------

  inface = new INFACE[nsub];
  if( !inface )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (SUBDOM::Input - 4)" );


  // mark nodes and elements belonging to the subdomain pid ------------------------------

  for( int e=0; e<ne; e++ )
  {
    sdel[e].mark = false;

    if( sdel[e].sub == pid )
    {
      sdel[e].mark = true;

      for( int i=1; i<=con[0][e]; i++ )
      {
        int n = con[i][e];

        sdnd[n].sub  = pid;
        sdnd[n].mark = true;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // set up a list of links from nodes to interfaces
  // link[    0][<node->no>]  holds the number of interfaces
  // link[1,...][<node->no>]  holds the interface id

  // check for existence of interfaces�---------------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    if( sdel[e].sub != pid )
    {
      for( int i=1; i<=con[0][e]; i++ )
      {
        int n = con[i][e];

        if( sdnd[n].sub == pid )  inface[sdel[e].sub].exist = true;
      }
    }
  }


  // count number of interfaces and allocate memory for list of links --------------------

  int nlink = 0;

  for( int s=0; s<nsub; s++ )
  {
    if( inface[s].exist )  nlink++;
  }

  int** link = new int*[nlink+1];
  if( !link )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (SUBDOM::Input - 5)" );

  link[0] = new int[(nlink+1)*np];
  if( !link[0] )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (SUBDOM::Input - 6)" );

  for( int i=1; i<=nlink; i++ )  link[i]    = link[0] + i * np;
  for( int n=0; n<np; n++ )      link[0][n] = 0;


  // set up list of links ----------------------------------------------------------------

  for( int e=0; e<ne; e++ )
  {
    if( sdel[e].sub != pid )
    {
      for( int i=1; i<=con[0][e]; i++ )
      {
        int n = con[i][e];

        if( sdnd[n].sub == pid )        // interface node to current subdomain pid
        {
          for( int l=1; l<=link[0][n]; l++ )
          {
            if( link[l][n] == sdel[e].sub )
            {
              n = -1;
              break;
            }
          }

          if( n >= 0 )
          {
            link[0][n]++;
            link[link[0][n]][n] = sdel[e].sub;
          }
        }
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // mark interface nodes: set "sdnd[n].sub" to the largest pid of attached subdomains

  for( int e=0; e<ne; e++ )
  {
    if( sdel[e].sub != pid )
    {
      for( int i=1; i<=con[0][e]; i++ )
      {
        int n = con[i][e];

        if( sdnd[n].sub == pid )
        {
          sdnd[n].sub = sdel[e].sub;
        }
        else if( sdnd[n].mark && sdnd[n].sub != sdel[e].sub )
        {
          if( sdnd[n].sub < sdel[e].sub )  sdnd[n].sub = sdel[e].sub;
        }
      }
    }
  }


  // set "sdnd[n].flag" ------------------------------------------------------------------

  for( int n=0; n<np; n++ )
  {
    if( sdnd[n].mark  &&  sdnd[n].sub != pid )
    {
      SF( sdnd[n].flag, SD_NODE::kInface );

           if( sdnd[n].sub < pid )  SF( sdnd[n].flag, SD_NODE::kInface_UP );
      else if( sdnd[n].sub > pid )  SF( sdnd[n].flag, SD_NODE::kInface_DN );
    }
  }


  // -------------------------------------------------------------------------------------
  // count nodes and elements in subdomain pid

  npdom = 0;
  nedom = 0;

  for( int n=0; n<np; n++ )
  {
    if( sdnd[n].mark )  npdom++;
  }

  for( int e=0; e<ne; e++ )
  {
    if( sdel[e].mark )  nedom++;
  }


  // print information on size of subdomains
# ifdef _MPI_
  int nxdom[2];

  if( pid == 0 )
  {
    MPI_Status status;

    sprintf( text, "\n %s\n %s  %5d    %7d  %7d\n",
                   "(SUBDOM::Input)         subdomain  elements   nodes",
                   "                        ", pid+1, nedom, npdom );
    REPORT::rpt.Output( text, 3 );

    for( int s=1; s<npr; s++ )
    {
      MPI_Recv( nxdom, 2, MPI_INT, s, 1, MPI_COMM_WORLD, &status );

      sprintf( text, " %s  %5d    %7d  %7d\n",
                     "                        ", s+1, nxdom[0], nxdom[1] );
      REPORT::rpt.Output( text, 3 );
    }
  }
  else
  {
    nxdom[0] = nedom;
    nxdom[1] = npdom;
    MPI_Send( nxdom, 2, MPI_INT, 0, 1, MPI_COMM_WORLD );
  }
# endif

  // -------------------------------------------------------------------------------------
  // allocate memory for nodes and elements

  region->Alloc( npdom, nedom );


  // -------------------------------------------------------------------------------------
  // ordering: 1. interior nodes             (nd->sub  == pid)
  //           2. upstream interface nodes   (nd->sub  <  pid)
  //           3. downstream interface nodes (nd->sub  >  pid)
  //
  // NOTE: A reordering of equations will be also performed later in EQS::ResetEqOrder().

  npdom = 0;

  for( int n=0; n<np; n++ )
  {
    node[n] = NULL;

    if( sdnd[n].mark  &&  sdnd[n].sub == pid )
    {
      node[n] = region->Getnode( npdom );
      npdom++;
    }
  }

  for( int n=0; n<np; n++ )
  {
    if( sdnd[n].mark  &&  isFS(sdnd[n].flag, SD_NODE::kInface_UP) )
    {
      node[n] = region->Getnode( npdom );
      npdom++;
    }
  }

  for( int n=0; n<np; n++ )
  {
    if( sdnd[n].mark  &&  isFS(sdnd[n].flag, SD_NODE::kInface_DN) )
    {
      node[n] = region->Getnode( npdom );
      npdom++;
    }
  }


  for( int n=0; n<np; n++ )
  {
    NODE* nd = node[n];

    if( nd )
    {
      if( isFS(sdnd[n].flag, SD_NODE::kInface) )     SF( nd->flag, NODE::kInface );
      if( isFS(sdnd[n].flag, SD_NODE::kInface_DN) )  SF( nd->flag, NODE::kInface_DN );
      if( isFS(sdnd[n].flag, SD_NODE::kInface_UP) )  SF( nd->flag, NODE::kInface_UP );
    }
  }


  // set pointer to elements in region file ----------------------------------------------

  nedom = 0;

  for( int e=0; e<ne; e++ )
  {
    elem[e]= NULL;

    if( sdel[e].mark )
    {
      elem[e] = region->Getelem( nedom );
      nedom++;
    }
  }


  // -------------------------------------------------------------------------------------
  // create list of interface nodes

  sprintf( text, "\n (SUBDOM::Input)         number of interface nodes...\n");
  REPORT::rpt.Output( text, 3 );


  // count number of nodes on interfaces -------------------------------------------------

  for( int n=0; n<np; n++ )
  {
    for( int l=1; l<=link[0][n]; l++ )
    {
      inface[link[l][n]].np++;
    }
  }


  // allocate memory for node pointers ---------------------------------------------------

  for( int s=0; s<nsub; s++ )
  {
    if( inface[s].np > 0 )
    {
      int npinf = inface[s].np;
      int nvinf = kSimDF * npinf;

      inface[s].node = new NODE* [npinf];

      inface[s].recv = new double[nvinf];
      inface[s].send = new double[nvinf];
      inface[s].ria1 = new char  [nvinf];
      inface[s].sia1 = new char  [nvinf];
      inface[s].ria2 = new int   [nvinf];
      inface[s].sia2 = new int   [nvinf];

      if( !inface[s].node  || !inface[s].recv || !inface[s].send
                           || !inface[s].ria1 || !inface[s].ria2
                           || !inface[s].sia1 || !inface[s].sia2 )
        REPORT::rpt.Error( kMemoryFault, "can not allocate memory (SUBDOM::Input - 11)" );

      sprintf( text, "                         %03d <- %03d;    np = %5d\n",
                     pid+1, s+1, inface[s].np );
      REPORT::rpt.Output( text, 3 );

      inface[s].np = 0;
    }
  }


  // set list of interface nodes ---------------------------------------------------------

  for( int n=0; n<np; n++ )
  {
    NODE* nd = node[n];

    if( nd )
    {
      for( int l=1; l<=link[0][n]; l++ )
      {
        int s = link[l][n];

        inface[s].node[inface[s].np] = nd;
        inface[s].np++;
      }
    }
  }

# ifdef kDebug
  {
    for( int s=0; s<nsub; s++ )
    {
      if( inface[s].np > 0 )
      {
        char fname[20];
        sprintf( fname, "subdom_%03d_%03d.geo", pid+1, s+1 );
        FILE* dbg = fopen( fname, "w" );

        fprintf( dbg, "1\nPOLYLINE %d\n", inface[s].np );

        for( int n=0; n<inface[s].np; n++ )
        {
          fprintf( dbg, "%15.8le %15.8le %15.8le\n",
                        inface[s].node[n]->x,
                        inface[s].node[n]->y,
                        inface[s].node[n]->z );
        }

        fclose( dbg );
      }
    }
  }
# endif

  // -------------------------------------------------------------------------------------
  // free allocated memory (1)

  for( int i=0; i<=kMaxNodes2D; i++ )  MEMORY::memo.Detach( con[i] );

  delete[] sdnd;
  delete[] sdel;


  // -------------------------------------------------------------------------------------
  // copy information of interfaces from link[][] to subdomain pointer NODE::sub

  int nl = 0;

  for( int i=0; i<np; i++ )
  {
    if( node[i] )  nl += link[0][i];
  }

  subbuf = new SUB[nl];
  if( !subbuf )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (SUBDOM::Input - 12)" );

  nl = 0;

  for( int i=0; i<np; i++ )
  {
    if( node[i] )
    {
      SUB* sub = NULL;

      for( int l=1; l<=link[0][i]; l++ )
      {
        if( l == 1 )
        {
          sub = node[i]->sub = &subbuf[nl];
          nl++;
        }
        else
        {
          sub = sub->next = &subbuf[nl];
          nl++;
        }

        sub->no = link[l][i];
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // free allocated memory (2)

  delete[] link[0];
  delete[] link;
}


//////////////////////////////////////////////////////////////////////////////////////////

void SUBDOM::SetInface( GRID* region )
{
  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    CF( nd->flag, NODE::kInface );
    CF( nd->flag, NODE::kInface_DN );
    CF( nd->flag, NODE::kInface_UP );

    int inf = false;
    int idn = false;

    SUB* sub = nd->sub;
    while( sub )
    {
      // check if at least one domain adjacent to "nd" is wet
      if( !sub->dry )
      {
        inf = true;

        // if the adjacent domain has a higher pid than me, it's a downstream interface
        if( sub->no > pid )  idn = true;
      }

      sub = sub->next;
    }

    if( inf )
    {
      SF( nd->flag, NODE::kInface );

      if( idn )  SF( nd->flag, NODE::kInface_DN );
      else       SF( nd->flag, NODE::kInface_UP );
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

int SUBDOM::Mpi_max( int num )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    int local = num;
    MPI_Allreduce( &local, &num, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  }
# endif

  return num;
}

double SUBDOM::Mpi_max( double num )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    double local = num;
    MPI_Allreduce( &local, &num, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  }
# endif

  return num;
}

double SUBDOM::Mpi_maxabs( double num )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    if( pid == 0 )
    {
      MPI_Status status;
      double mpi_num;

      for( int s=1; s<npr; s++ )
      {
        MPI_Recv( &mpi_num, 1, MPI_DOUBLE, s, 1, MPI_COMM_WORLD, &status );
        if( fabs(mpi_num) > fabs(num) ) num = mpi_num;
      }
    }
    else
    {
      MPI_Send( &num, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
    }

    MPI_Bcast( &num, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  }
# endif

  return num;
}

int SUBDOM::Mpi_min( int num )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    int local = num;
    MPI_Allreduce( &local, &num, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
  }
# endif

  return num;
}

double SUBDOM::Mpi_min( double num )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    double local = num;
    MPI_Allreduce( &local, &num, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  }
# endif

  return num;
}


int SUBDOM::Mpi_sum( int data )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    int local = data;
    MPI_Allreduce( &local, &data, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  }
# endif

  return data;
}

/*
int SUBDOM::Mpi_sum( int data )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    if( pid == 0 )
    {
      MPI_Status status;
      int mpi_data;

      for( int s=1; s<npr; s++ )
      {
        MPI_Recv( &mpi_data, 1, MPI_INT, s, 1, MPI_COMM_WORLD, &status );
        data += mpi_data;
      }
    }
    else
    {
      MPI_Send( &data, 1, MPI_INT, 0, 1, MPI_COMM_WORLD );
    }

    MPI_Bcast( &data, 1, MPI_INT, 0, MPI_COMM_WORLD );
  }
# endif

  return data;
}
*/

double SUBDOM::Mpi_sum( double data )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    double local = data;
    MPI_Allreduce( &local, &data, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  }
# endif

  return data;
}

/*
double SUBDOM::Mpi_sum( double data )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    if( pid == 0 )
    {
      MPI_Status status;
      double mpi_data;

      for( int s=1; s<npr; s++ )
      {
        MPI_Recv( &mpi_data, 1, MPI_DOUBLE, s, 1, MPI_COMM_WORLD, &status );
        data += mpi_data;
      }
    }
    else
    {
      MPI_Send( &data, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
    }

    MPI_Bcast( &data, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  }
# endif

  return data;
}
*/

//////////////////////////////////////////////////////////////////////////////////////////
// Assemble the nodal vector "vec[]" of length "npdom" across subdomains.

void SUBDOM::Mpi_assemble( double* vec )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    // copy vector "vec[]" to a temporary array ------------------------------------------
    double* tmp = (double*) MEMORY::memo.Array_nd( npdom );
    memcpy( tmp, vec, npdom*sizeof(double) );

    // loop on all interfaces 0 <= s < npr: exchange vector data -------------------------
    for( int s=0; s<npr; s++ )
    {
      MPI_Status status;

      int npinf = inface[s].np;             // number of nodes on interface to domain s

      if( npinf > 0 )
      {
        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];     // pointer to the interface node n
          int   no = nd->Getno();           // node number local in this domain

          inface[s].send[n] = tmp[no];      // copy the vector element to the send array
        }

        MPI_Sendrecv( inface[s].send, npinf, MPI_DOUBLE, s, 1,
                      inface[s].recv, npinf, MPI_DOUBLE, s, 1,
                      MPI_COMM_WORLD, &status );

        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];
          int   no = nd->Getno();

          vec[no] += inface[s].recv[n];
        }
      }
    }

    MEMORY::memo.Detach( tmp );
  }
# endif
}


void SUBDOM::Mpi_assemble( int* cnt )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    // copy vector "vec[]" to a temporary array ------------------------------------------
    int* tmp = (int*) MEMORY::memo.Array_nd( npdom );
    memcpy( tmp, cnt, npdom*sizeof(int) );

    // loop on all interfaces 0 <= s < npr: exchange vector data -------------------------
    for( int s=0; s<npr; s++ )
    {
      MPI_Status status;

      int npinf = inface[s].np;             // number of nodes on interface to domain s

      if( npinf > 0 )
      {
        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];     // pointer to the interface node n
          int   no = nd->Getno();           // node number local in this domain

          inface[s].sia2[n] = tmp[no];      // copy the vector element to the send array
        }

        MPI_Sendrecv( inface[s].sia2, npinf, MPI_INT, s, 1,
                      inface[s].ria2, npinf, MPI_INT, s, 1,
                      MPI_COMM_WORLD, &status );

        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];
          int   no = nd->Getno();

          cnt[no] += inface[s].ria2[n];
        }
      }
    }

    MEMORY::memo.Detach( tmp );
  }
# endif
}


//////////////////////////////////////////////////////////////////////////////////////////
// Average the nodal vector "vec[]" of length "npdom" across subdomains.

void SUBDOM::Mpi_average( double* vec )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    // copy vector "vec[]" to a temporary array ------------------------------------------
    int*    cnt = (int*)    MEMORY::memo.Array_nd( npdom );
    double* tmp = (double*) MEMORY::memo.Array_nd( npdom );

    memcpy( tmp, vec, npdom*sizeof(double) );

    for( int n=0; n<npdom; n++ ) cnt[n] = 1;

    // loop on all interfaces: exchange vector data --------------------------------------
    for( int s=0; s<npr; s++ )
    {
      MPI_Status status;

      int np = inface[s].np;

      if( np > 0 )
      {
        for( int n=0; n<np; n++ )
        {
          NODE* nd = inface[s].node[n];
          inface[s].send[n] = tmp[nd->Getno()];
        }

        MPI_Sendrecv( inface[s].send, np, MPI_DOUBLE, s, 1,
                      inface[s].recv, np, MPI_DOUBLE, s, 1,
                      MPI_COMM_WORLD, &status );

        for( int n=0; n<np; n++ )
        {
          NODE* nd = inface[s].node[n];
          vec[nd->Getno()] += inface[s].recv[n];
          cnt[nd->Getno()]++;
        }
      }
    }

    // average vector --------------------------------------------------------------------
    for( int n=0; n<npdom; n++ )
    {
      if( cnt[n] > 1 )  vec[n] /= cnt[n];
    }

    // detach memory ---------------------------------------------------------------------
    MEMORY::memo.Detach( cnt );
    MEMORY::memo.Detach( tmp );
  }
# endif
}


//////////////////////////////////////////////////////////////////////////////////////////
// Reduce the nodal vector "vec[]" of length "npdom" to its maximum across subdomains.

void SUBDOM::Mpi_max( double* vec )
{
# ifdef _MPI_
  if( npr > 1 )
  {
    // copy vector "vec[]" to a temporary array ------------------------------------------
    double* tmp = (double*) MEMORY::memo.Array_nd( npdom );
    memcpy( tmp, vec, npdom*sizeof(double) );

    // loop on all interfaces 0 <= s < npr: exchange vector data -------------------------
    for( int s=0; s<npr; s++ )
    {
      MPI_Status status;

      int npinf = inface[s].np;             // number of nodes on interface to domain s

      if( npinf > 0 )
      {
        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];     // pointer to the interface node n
          int   no = nd->Getno();           // node number local in this domain

          inface[s].send[n] = tmp[no];      // copy the vector element to the send array
        }

        MPI_Sendrecv( inface[s].send, npinf, MPI_DOUBLE, s, 1,
                      inface[s].recv, npinf, MPI_DOUBLE, s, 1,
                      MPI_COMM_WORLD, &status );

        for( int n=0; n<npinf; n++ )
        {
          NODE* nd = inface[s].node[n];
          int   no = nd->Getno();

          if( inface[s].recv[n] > vec[no] )  vec[no] = inface[s].recv[n];
        }
      }
    }

    MEMORY::memo.Detach( tmp );
  }
# endif
}
