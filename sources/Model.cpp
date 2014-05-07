// ======================================================================================
//
// Copyright (C) 1992-2012  by  P.M. SCHROEDER
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
#include "Shape.h"
#include "Time.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Memory.h"
#include "Project.h"
#include "Statist.h"
#include "Tmserhead.h"

#include "Model.h"

//////////////////////////////////////////////////////////////////////////////////////////

MODEL::MODEL()
{
  np   = 0;
  node = NULL;

  ne   = 0;
  elem = NULL;

  list = NULL;

  init = 0;

  boundList = NULL;

  // allocate memory for grids -----------------------------------------------------------
  region  = new GRID();
  control = new GRID();
  bound   = new GRID();

  if( !region || !control || !bound )
    REPORT::rpt.Error( kMemoryFault, "can not allocate memory (MODEL::MODEL - 1)" );

  subdom = NULL;

  // initialize variables for output
  phi   = NULL;
  rot   = NULL;
  curv  = NULL;
  Lx    = NULL;
  Ly    = NULL;
  rc    = NULL;
  d90   = NULL;
  d50   = NULL;
  Hr    = NULL;
  Hd    = NULL;
  kd    = NULL;
  hp    = NULL;
  dp    = NULL;
  sp    = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////

MODEL::~MODEL()
{
  delete region;
  delete control;
  delete bound;
}

//////////////////////////////////////////////////////////////////////////////////////////

int MODEL::Getinit()
{
  return init;
}

//////////////////////////////////////////////////////////////////////////////////////////

void MODEL::Incinit()
{
  init++;
}


//////////////////////////////////////////////////////////////////////////////////////////
// initialize model list of elements and nodes:
// used methods:  - connectivity()
//                - boundary()
//                - setLocation()

void MODEL::Initialize()
{
  init++;

  // set up connection of nodes to elements ----------------------------------------------
  region->Connection( ELEM::kDry );

  // compute boundary elements -----------------------------------------------------------
  Boundary();

  // delete existing list of nodes and elements ------------------------------------------
  if( node )  delete[] node;
  if( elem )  delete[] elem;

  node = NULL;
  elem = NULL;

  // create list of nodes ----------------------------------------------------------------
  np = 0;

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);
    if( !isFS(nd->flag, NODE::kDry) )  np++;
  }

  if( np )
  {
    node = new NODE* [np];

    if( !node )
      REPORT::rpt.Error( kMemoryFault, "can not allocate memory (MODEL::initialize - 1)" );

    np = 0;

    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE* nd = region->Getnode(n);

      if( !isFS(nd->flag,NODE::kDry) )
      {
        node[np] = nd;
        np++;
      }
    }
  }

  // create list of elements (including boundary elements) -------------------------------
  ne = 0;

  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* el = region->Getelem(e);

    if( !isFS(el->flag, ELEM::kDry) )  ne++;
  }

  ne += bound->Getne();


  if( ne )
  {
    elem = new ELEM* [ne];

    if( !elem )
      REPORT::rpt.Error( kMemoryFault, "can not allocate memory (MODEL::initialize - 2)" );

    ne = 0;

    for( int re=0; re<region->Getne(); re++ )
    {
      ELEM* el = region->Getelem(re);

      if( !isFS(el->flag, ELEM::kDry) )
      {
        elem[ne] = el;                                  // set pointer to region element
        ne++;

        int ncn = el->Getncn();
        int nnd = el->Getnnd();

        for( int i=ncn; i<nnd; i++ )
        {
          if( isFS(el->nd[i]->flag, NODE::kBound) )
          {
            elem[ne] = boundList[el->nd[i]->Getno()];   // set pointer to boundary element
            ne++;
          }
        }
      }
    }
  }

  // set location of boundary elements and nodes -----------------------------------------
  SetLocation();
}

//////////////////////////////////////////////////////////////////////////////////////////

int MODEL::Input( PROJECT* project, NAME* name, TIME* prevTime, SED* sed )
{
  char text[200];


  // MPI: read subdomain file and determine sub domain nodes and elements ----------------
  if( subdom->npr > 1 )
  {
    subdom->Input( name->subdomFile, name->regionFile, region );
  }

  // read region elements ----------------------------------------------------------------
  region->InputRegion( name->regionFile, subdom );

  sprintf( text, "\n %s %d\n %s %d\n",
           "(MODEL::input)          number of nodes            :", region->Getnp(),
           "                        number of region elements  :", region->Getne() );
  REPORT::rpt.Output( text, 2 );

  // set region-flag ---------------------------------------------------------------------
  int nType = TYPE::Get();

  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* el = region->Getelem(e);

    SF( el->flag, ELEM::kRegion );
    //if( el->type < 0  ||  el->type >= nType )
    //{
    //  REPORT::rpt.Error( kInputFault, "bad element type number %d %s\n",
    //                                  el->type, "(MODEL::Input #1)" );
    //}
  }

  // set flag for corner and midside nodes -----------------------------------------------
  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* el = region->Getelem(e);

    int ncn = el->Getncn();
    int nnd = el->Getnnd();

    int i;

    for( i=0;   i<ncn; i++ )  SF( el->nd[i]->flag, NODE::kCornNode );
    for( i=ncn; i<nnd; i++ )  SF( el->nd[i]->flag, NODE::kMidsNode );
  }

  int ncp = 0;

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    if( isFS(nd->flag, NODE::kCornNode) ) ncp++;
  }

  sprintf( text, "\n %s %d\n",
           "(MODEL::input)          number of corner nodes     :", ncp );
  REPORT::rpt.Output( text, 2 );

  // read control elements ---------------------------------------------------------------
  if( name->controlFile[0] )
  {
    control->InputControl( name->controlFile, region, subdom );

    sprintf( text, "\n %s %d\n",
             "(MODEL::input)          number of control elements :", control->Getne() );
    REPORT::rpt.Output( text, 2 );

    for( int e=0; e<control->Getne(); e++ )  SF( control->Getelem(e)->flag, ELEM::kControl );
  }

  // initialize node data ----------------------------------------------------------------
  region->InputInitial( name->ascii_initial, name->initialFile, prevTime, subdom, sed->zb_init );

  if( project->timeint.set )  *prevTime = project->timeint.startTime;

  // initialize list of nodes and elements -----------------------------------------------
  Initialize();

  // mark boundary nodes of the whole grid with the flag NODE::kGridBnd ------------------
  // this flag will never change
  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);
    if( isFS(nd->flag, NODE::kBound) ) SF( nd->flag, NODE::kGridBnd );
  }

  if( name->initialFile[0] )  return 1;
  else                        return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////

void MODEL::Output( PROJECT* project, int timeStep )
{
  // count number of wet corner nodes and elements ---------------------------------------
  int npwet = 0;
  int npc   = 0;

  for( int n=0; n<region->Getnp(); n++ )
  {
    NODE* nd = region->Getnode(n);

    if( isFS(nd->flag, NODE::kCornNode) )
    {
      npc++;

      if( !isFS(nd->flag, NODE::kDry) ) npwet++;
    }
  }

  int newet = 0;

  for( int e=0; e<region->Getne(); e++ )
  {
    ELEM* el = region->Getelem(e);

    if( !isFS(el->flag, ELEM::kDry) ) newet++;
  }

  // -------------------------------------------------------------------------------------
  // compute some special parameter like rotation of flow, element size, ...
  for( int i=0; i<project->vpcomp; i++ )
  {
    AttachOutput( project, project->vpoutlist[i] );
  }

  // -------------------------------------------------------------------------------------
  // write region file
  if( project->name.rgAvsFile[0] )
  {
    char fileName[500];

    project->ReplaceMacro( project->name.rgAvsFile, fileName, timeStep, 0 );

    FILE* id = fopen( fileName, "w" );
    if( !id )
      REPORT::rpt.Error( kOpenFileFault, "%s %s (MODEL::output - 1)",
                         "can not open file", fileName );

    // write number of nodes, elements and node data -------------------------------------
    fprintf( id, "%d  %d  %d  %d  %d\n", region->Getnp(), region->Getne(),
             project->vpdata, project->vedata, 0 );

    // write nodes -----------------------------------------------------------------------
    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE* nd = region->Getnode(n);

      double X = nd->x;
      double Y = nd->y;
      double Z = nd->zor;

      fprintf( id, "%7d  %12.6lf  %12.6lf  %12.6lf\n", nd->Getname(), X, Y, Z );

      nd->mark = true;
    }

    // write element connectivity --------------------------------------------------------
    for( int e=0; e<region->Getne(); e++ )
    {
      ELEM* el = region->Getelem(e);

      char elemShape[6];

      switch( el->GetshapeID() )
      {
        case kLine:      strcpy( elemShape, "line" );   break;
        case kTriangle:  strcpy( elemShape, "tri" );    break;
        case kSquare:    strcpy( elemShape, "quad" );   break;
      }

      int mat = TYPE::Getid(el->type)->no(el->type);
      if( mat <= 0 )  mat = abs(el->type);

      fprintf( id, "%6d  %5d  %5s", el->Getname(), mat, elemShape );

      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )  fprintf( id, "  %6d", el->nd[i]->Getname() );

      fprintf( id, "\n" );

      el->mark = true;
    }

    // write nodal values ----------------------------------------------------------------
    UCDOutput( id, project, region, phi, rot, curv, Lx, Ly,
               rc, d90, d50, Hr, Hd, kd, hp, dp, sp,
               project->sed.qbc, project->sed.PLs, project->sed.sx, project->sed.sy,
               project->sed.dzds, project->sed.dzmx, project->sed.dhds );

    fclose( id );
  }

  // -------------------------------------------------------------------------------------
  // write reduced region file (only corner nodes)
  if( project->name.redAvsFile[0] )
  {
    char fileName[500];

    project->ReplaceMacro( project->name.redAvsFile, fileName, timeStep, 0 );

    FILE* id = fopen( fileName, "w" );
    if( !id )
      REPORT::rpt.Error( kOpenFileFault, "%s %s (MODEL::output - 1)",
                         "can not open file", fileName );

    // write number of nodes, elements and node data -------------------------------------
    fprintf( id, "%d  %d  %d  %d  %d\n", npc, region->Getne(),
             project->vpdata, project->vedata, 0 );

    // write nodes -----------------------------------------------------------------------
    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE* nd = region->Getnode(n);

      if( isFS(nd->flag, NODE::kCornNode) )
      {
        double X = nd->x;
        double Y = nd->y;
        double Z = nd->zor;

        fprintf( id, "%6d  %12.6lf  %12.6lf  %12.6lf\n", nd->Getname(), X, Y, Z );

        nd->mark = true;
      }
      else
      {
        nd->mark = false;
      }
    }

    // write element connectivity --------------------------------------------------------
    for( int e=0; e<region->Getne(); e++ )
    {
      ELEM* el = region->Getelem(e);

      char elemShape[6];

      switch( el->GetshapeID() )
      {
        case kLine:      strcpy( elemShape, "line" );   break;
        case kTriangle:  strcpy( elemShape, "tri" );    break;
        case kSquare:    strcpy( elemShape, "quad" );   break;
      }

      int mat = TYPE::Getid(el->type)->no(el->type);
      if( mat <= 0 )  mat = abs(el->type);

      fprintf( id, "%6d  %5d  %5s", el->Getname(), mat, elemShape );

      int ncn = el->Getncn();

      for( int i=0; i<ncn; i++ )  fprintf( id, "  %6d", el->nd[i]->Getname() );

      fprintf( id, "\n" );

      el->mark = true;
    }

    // write nodal values ----------------------------------------------------------------
    UCDOutput( id, project, region, phi, rot, curv, Lx, Ly,
               rc, d90, d50, Hr, Hd, kd, hp, dp, sp,
               project->sed.qbc, project->sed.PLs, project->sed.sx, project->sed.sy,
               project->sed.dzds, project->sed.dzmx, project->sed.dhds );

    fclose( id );
  }

  // -------------------------------------------------------------------------------------
  // write region file (only dry corner nodes)
  if( project->name.wetAvsFile[0] )
  {
    char fileName[500];

    project->ReplaceMacro( project->name.wetAvsFile, fileName, timeStep, 0 );

    FILE* id = fopen( fileName, "w" );
    if( !id )
      REPORT::rpt.Error( kOpenFileFault, "can not open file %s (MODEL::output - 2)", fileName );

    // write number of nodes, elements and node data -------------------------------------
    fprintf( id, "%d  %d  %d  %d  %d\n", npwet, newet,
             project->vpdata, project->vedata, 0 );


    // write nodes -----------------------------------------------------------------------
    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE* nd = region->Getnode(n);

      if(    !isFS(nd->flag, NODE::kDry)
             &&  isFS(nd->flag, NODE::kCornNode) )
      {
        double X = nd->x;
        double Y = nd->y;
        double Z = nd->z;

        fprintf( id, "%6d  %12.6lf  %12.6lf  %12.6lf\n", nd->Getname(), X, Y, Z );

        nd->mark = true;
      }
      else
      {
        nd->mark = false;
      }
    }

    // write element connectivity --------------------------------------------------------
    for( int e=0; e<region->Getne(); e++ )
    {
      ELEM* el = region->Getelem(e);

      char elemShape[6];

      if( !isFS(el->flag, ELEM::kDry) )
      {
        switch( el->GetshapeID() )
        {
          case kLine:      strcpy( elemShape, "line" );   break;
          case kTriangle:  strcpy( elemShape, "tri" );    break;
          case kSquare:    strcpy( elemShape, "quad" );   break;
        }

        int mat = TYPE::Getid(el->type)->no(el->type);

        fprintf( id, "%6d  %5d  %5s", el->Getname(), mat, elemShape );

        int ncn = el->Getncn();

        for( int i=0; i<ncn; i++ )
        {
          fprintf( id, "  %6d", el->nd[i]->Getname() );
        }

        fprintf( id, "\n" );

        el->mark = true;
      }
      else
      {
        el->mark = false;
      }
    }

    // write nodal values ----------------------------------------------------------------
    UCDOutput( id, project, region, phi, rot, curv, Lx, Ly,
               rc, d90, d50, Hr, Hd, kd, hp, dp, sp,
               project->sed.qbc, project->sed.PLs, project->sed.sx, project->sed.sy,
               project->sed.dzds, project->sed.dzmx, project->sed.dhds );

    fclose( id );
  }

  // -------------------------------------------------------------------------------------
  // write control file
  if( project->name.ctAvsFile[0] )
  {
    char fileName[500];

    project->ReplaceMacro( project->name.ctAvsFile, fileName, timeStep, 0 );

    FILE* id = fopen( fileName, "w" );
    if( !id )
      REPORT::rpt.Error( kOpenFileFault, "%s %s (MODEL::output - 4)",
                         "can not open file", fileName );

    // mark and count nodes at control elements ------------------------------------------
    for( int n=0; n<control->Getnp(); n++ )
    {
      NODE* nd = control->Getnode(n);
      nd->mark = false;
    }

    for( int e=0; e<control->Getne(); e++ )
    {
      ELEM* el = control->Getelem(e);

      int nnd = el->Getnnd();
      for( int i=0; i<nnd; i++ )  el->nd[i]->mark = true;
    }

    int nc = 0;

    for( int n=0; n<control->Getnp(); n++ )
    {
      NODE* nd = control->Getnode(n);
      if( nd->mark )  nc++;
    }

    // write number of nodes, elements and node data -------------------------------------
    fprintf( id, "%d  %d  %d  %d  %d\n", nc, control->Getne(),
             project->vpdata, 0, 0 );

    // write nodes -----------------------------------------------------------------------
    for( int n=0; n<control->Getnp(); n++ )
    {
      NODE* nd = control->Getnode(n);

      if( nd->mark )
      {
        double X = nd->x;
        double Y = nd->y;
        double Z = nd->zor;

        fprintf( id, "%6d  %12.6lf  %12.6lf  %12.6lf\n",  nd->Getname(), X, Y, Z );
      }
    }

    // write element connectivity --------------------------------------------------------
    for( int e=0; e<control->Getne(); e++ )
    {
      ELEM* el = control->Getelem(e);

      char elemShape[6];

      switch( el->GetshapeID() )
      {
        case kLine:      strcpy(elemShape, "line");   break;
        case kTriangle:  strcpy(elemShape, "tri");    break;
        case kSquare:    strcpy(elemShape, "quad");   break;
      }

      fprintf( id, "%6d  %5d  %5s", el->Getname(), el->type, elemShape );

      int nnd = el->Getnnd();

      for( int i=0; i<nnd; i++ )
      {
        fprintf( id, "  %6d", el->nd[i]->Getname() );
      }

      fprintf( id, "\n" );

      el->mark = true;
    }

    // write nodal values ----------------------------------------------------------------
    UCDOutput( id, project, control, phi, rot, curv, Lx, Ly,
               rc, d90, d50, Hr, Hd, kd, hp, dp, sp,
               project->sed.qbc, project->sed.PLs, project->sed.sx, project->sed.sy,
               project->sed.dzds, project->sed.dzmx, project->sed.dhds );

    fclose( id );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void MODEL::OutputSeriesHeader( PROJECT *project, TMSER *tmser )
{
  // count number of corner nodes and mark them ------------------------------------------
  int npc = 0;

  for( int n=0; n<region->Getnp(); n++ )
  {
    if( isFS(region->Getnode(n)->flag,NODE::kCornNode) ) npc++;
  }

  // write the header and the grid geometry to the file (overwrite existing files)--------
  TMSERHEAD header;
  TMSERHEAD header_old;

  FILE *id = fopen( tmser->filename, "rb+" );
  if( !id )
  {
    id = fopen( tmser->filename, "wb+" );
    if( !id ) REPORT::rpt.Error( kOpenFileFault, "%s %s (MODEL::OutputSeriesHeader #1)",
                                 "cannot open file", tmser->filename );

    header.np = npc;
    header.ne = region->Getne();

    header.first = tmser->first;
    header.last  = 0;
    header.step  = tmser->step;

    project->timeint.startTime.Get( header.startTime );
    project->timeint.deltaTime.Get( header.deltaTime );

    header.vcomp = tmser->vcomp;
    header.vdata = tmser->vdata;

    fwrite( header.buffer, sizeof(char), kTMSERHEADSize, id );

    // write idents of variables and their vector length
    for( int i=0; i<header.vcomp; i++ )
    {
      int  dim;
      char ident[10];

      int j = tmser->voutlist[i];

      dim = project->valist[j].dim;
      strcpy( ident, project->valist[j].name );

      fwrite( &dim,  sizeof(int),   1, id );
      fwrite( ident, sizeof(char), 10, id );
    }

    // write node names and coordinates (x,y,z) --------------------------------------------
    int    *ndid = (int*)    MEMORY::memo.Array_nd( region->Getnp() );
    double *xcor = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    double *ycor = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    double *zcor = (double*) MEMORY::memo.Array_nd( region->Getnp() );

    npc = 0;
    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE *nd = region->Getnode(n);

      if( isFS(nd->flag,NODE::kCornNode) )
      {
        ndid[npc] = nd->Getname();
        xcor[npc] = nd->x;
        ycor[npc] = nd->y;
        zcor[npc] = nd->zor;
        npc++;
      }
    }

    fwrite( ndid, sizeof(int),    header.np, id );
    fwrite( xcor, sizeof(double), header.np, id );
    fwrite( ycor, sizeof(double), header.np, id );
    fwrite( zcor, sizeof(double), header.np, id );

    // detach the temporarily used arrays --------------------------------------------------
    MEMORY::memo.Detach( ndid );
    MEMORY::memo.Detach( xcor );
    MEMORY::memo.Detach( ycor );
    MEMORY::memo.Detach( zcor );

    // write element names, materials and connectivity -------------------------------------
    int *elid = (int*) MEMORY::memo.Array_el( region->Getne() );
    int *type = (int*) MEMORY::memo.Array_el( region->Getne() );
    int *nd1  = (int*) MEMORY::memo.Array_el( region->Getne() );
    int *nd2  = (int*) MEMORY::memo.Array_el( region->Getne() );
    int *nd3  = (int*) MEMORY::memo.Array_el( region->Getne() );
    int *nd4  = (int*) MEMORY::memo.Array_el( region->Getne() );

    for( int e=0; e<region->Getne(); e++ )
    {
      ELEM *el = region->Getelem(e);

      elid[e] = el->Getname();
      type[e] = el->type;

      switch( el->shape )
      {
        case kTri:
          nd1[e] = el->Getnode(0)->Getname();
          nd2[e] = el->Getnode(1)->Getname();
          nd3[e] = el->Getnode(2)->Getname();
          nd4[e] = 0;
          break;

        case kQuad:
          nd1[e] = el->Getnode(0)->Getname();
          nd2[e] = el->Getnode(1)->Getname();
          nd3[e] = el->Getnode(2)->Getname();
          nd4[e] = el->Getnode(3)->Getname();
          break;
      }
    }

    fwrite( elid, sizeof(int), header.ne, id );
    fwrite( type, sizeof(int), header.ne, id );
    fwrite( nd1,  sizeof(int), header.ne, id );
    fwrite( nd2,  sizeof(int), header.ne, id );
    fwrite( nd3,  sizeof(int), header.ne, id );
    fwrite( nd4,  sizeof(int), header.ne, id );

    // detach the temporarily used arrays --------------------------------------------------
    MEMORY::memo.Detach( elid );
    MEMORY::memo.Detach( type );
    MEMORY::memo.Detach( nd1 );
    MEMORY::memo.Detach( nd2 );
    MEMORY::memo.Detach( nd3 );
    MEMORY::memo.Detach( nd4 );

  }

  else if( fread( header_old.buffer, sizeof(char), kTMSERHEADSize, id ) < kTMSERHEADSize
          || header_old.vcomp != tmser->vcomp || header_old.vdata != tmser->vdata )
//          || header_old.step  != tmser->step )
  {
    REPORT::rpt.Error( kReadFileFault, "%s %s %s (MODEL::OutputSeries #2)",
                       "can not read file header", tmser->filename, "oder rts-Ausgaben nicht kompatibel" );
  }




  fclose( id );


}

//////////////////////////////////////////////////////////////////////////////////////////

void MODEL::OutputSeries( PROJECT *project, int timeStep, TMSER *tmser, bool tmser_first_call )
{
  // write the header, if this is the first time step
  //if( timeStep == tmser->first ) OutputSeriesHeader( project, tmser );

  //  static int firstCall = true;
  //  if( firstCalltmser[its])
  //  {
  //    tmser->first = timeStep;
  //    OutputSeriesHeader( project, tmser );
  //    firstCalltmser[its] = false;
  //  }

  if( !tmser_first_call )
  {
    tmser->first = timeStep;
    OutputSeriesHeader( project, tmser );
  }


  TMSERHEAD header;

  for( int i=0; i<tmser->vcomp; i++ )
  {
    AttachOutput( project, project->valist[tmser->voutlist[i]].id );
  }

  FILE *id = fopen( tmser->filename, "rb+" );
  if( !id ) REPORT::rpt.Error( kOpenFileFault, "%s %s (MODEL::OutputSeries #1)",
                               "can not open file", tmser->filename );

  double *x = (double*) MEMORY::memo.Array_nd( region->Getnp() );
  double *y = (double*) MEMORY::memo.Array_nd( region->Getnp() );

  SED *sed = &project->sed;

  // read the file header
  fseek( id, 0L, SEEK_SET );

  if( fread( header.buffer, sizeof(char), kTMSERHEADSize, id ) < kTMSERHEADSize )
  {
    REPORT::rpt.Error( kReadFileFault, "%s %s (MODEL::OutputSeries #2)",
                       "can not read file header", tmser->filename );
  }

  if( header.last >= timeStep )
  {
//  header.last = timeStep;

//  fseek( id, 0L, SEEK_SET );
//  fwrite( header.buffer, sizeof(char), kTMSERHEADSize, id );

  // calculate actual right position in rts
  int ntimesteps = ( timeStep - header.first ) / tmser->step;
  unsigned long pos = kTMSERHEADSize
      + header.vcomp * (sizeof(int) + 10*sizeof(char) )
      + header.np * ( sizeof(int) + 3*sizeof(double) )
      + header.ne * ( 6 * sizeof(int) )
      + header.np * ( header.vdata * ntimesteps * sizeof(double) );

  fseek( id, pos, SEEK_SET);
  }
  else
  {
    header.last = timeStep;
    fseek( id, 0L, SEEK_SET );
    fwrite( header.buffer, sizeof(char), kTMSERHEADSize, id );
    fseek( id, 0L, SEEK_END);
  }

  for( int i=0; i<tmser->vcomp; i++ )
  {
    int npc = 0;

    for( int n=0; n<region->Getnp(); n++ )
    {
      NODE *nd = region->Getnode(n);

      if( isFS(nd->flag,NODE::kCornNode) )
      {
        switch( project->valist[tmser->voutlist[i]].id )
        {
          case PROJECT::kUV:
            x[npc] = nd->v.U;
            y[npc] = nd->v.V;
            break;

          case PROJECT::kS:
            x[npc] = region->Getnode(n)->v.S;
            break;

          case PROJECT::kK:
            x[npc] = region->Getnode(n)->v.K;
            break;

          case PROJECT::kD:
            x[npc] = region->Getnode(n)->v.D;
            break;

          case PROJECT::kC:
            x[npc] = region->Getnode(n)->v.C;
            break;

          case PROJECT::kCB:
          {
            double U  = nd->v.U;
            double V  = nd->v.V;
            double Us = sqrt( U*U + V*V );
            if( Us > 0.0 )  x[npc] = nd->v.Qb / Us;
            else            x[npc] = 0.0;
          }
            break;

          case PROJECT::kQB:
            x[npc] = region->Getnode(n)->v.Qb;
            break;

          case PROJECT::kQBE:
            x[npc] = sed->qbc[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kLS:
            x[npc] = sed->PLs[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kSEDUV:
          {
            int no = nd->Getno();
            x[npc] = nd->v.Qb*sed->sx[no];
            y[npc] = nd->v.Qb*sed->sy[no];
          }
            break;

          case PROJECT::kDUVDT:
            x[npc] = nd->v.dUdt;
            y[npc] = nd->v.dVdt;
            break;

          case PROJECT::kDSDT:
            x[npc] = region->Getnode(n)->v.dSdt;
            break;

          case PROJECT::kZ:
            x[npc] = region->Getnode(n)->z;
            break;

          case PROJECT::kH:
            x[npc] = nd->v.S - nd->z;
            break;

          case PROJECT::kUS:
          {
            double U  = nd->v.U;
            double V  = nd->v.V;
            x[npc] = sqrt( U*U + V*V );
          }
            break;

          case PROJECT::kCF:
            x[npc] = region->Getnode(n)->cf;
            break;

          case PROJECT::kRC:
            x[npc] = rc[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kD90:
            x[npc] = d90[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kD50:
            x[npc] = d50[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kHR:
            x[npc] = Hr[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kHD:
            x[npc] = Hd[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kKD:
            x[npc] = kd[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kHP:
            x[npc] = hp[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kDP:
            x[npc] = dp[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kSP:
            x[npc] = sp[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kTAU:
          {
            double H    = nd->v.S - nd->z;
            double U    = nd->v.U;
            double V    = nd->v.V;
            double Us   = sqrt( U*U + V*V );
            double Utau = sed->GetUtau( Us, H, 0.0, project );
            x[npc] = project->rho * Utau * Utau;
          }
            break;

          case PROJECT::kMAN:
          {
            double H  = nd->v.S - nd->z;
            if( H > 0.0 ) x[npc] = sqrt( nd->cf * pow(H,1.0/3.0) / project->g );
            else          x[npc] =  0.0;
          }
            break;

          case PROJECT::kVT:
            x[npc] = region->Getnode(n)->vt;
            break;

          case PROJECT::kEst:
          {
            double K   = nd->v.K;
            double D   = nd->v.D;
            double H   = nd->v.S - nd->z;
            double U   = nd->v.U;
            double V   = nd->v.V;
            double Ust = sqrt( nd->cf * (U*U + V*V) );
            double cm  = project->KD.cm;
            if( nd->v.D > 0.0  &&  Ust > 0.0  &&  H > 0.0 ) x[npc] = cm*K*K/D / Ust / H;
            else                                            x[npc] = 0.0;
          }
            break;

          case PROJECT::kEXX:
          {
            double r = sqrt( 4.0*nd->exy*nd->exy + (nd->exx-nd->eyy)*(nd->exx-nd->eyy) );
            x[npc] = 0.5*(nd->exx + nd->eyy + r) * nd->vt;
          }
            break;

          case PROJECT::kEYY:
          {
            double r = sqrt( 4.0*nd->exy*nd->exy + (nd->exx-nd->eyy)*(nd->exx-nd->eyy) );
            x[npc] = 0.5*(nd->exx + nd->eyy - r) * nd->vt;
          }
            break;

          case PROJECT::kDUU:
          {
            double H  = nd->v.S - nd->z;
            double U  = nd->v.U;
            double V  = nd->v.V;
            double Us = sqrt( U*U + V*V );
            x[npc] = H * Us * Us * nd->Dxx;
          }
            break;

          case PROJECT::kDUV:
          {
            double H  = nd->v.S - nd->z;
            double U  = nd->v.U;
            double V  = nd->v.V;
            double Us = sqrt( U*U + V*V );
            x[npc] = H * Us * Us * nd->Dxy;
          }
            break;

          case PROJECT::kDVV:
          {
            double H  = nd->v.S - nd->z;
            double U  = nd->v.U;
            double V  = nd->v.V;
            double Us = sqrt( U*U + V*V );
            x[npc] = H * Us * Us * nd->Dyy;
          }
            break;

          case PROJECT::kVSEC:
            x[npc] = region->Getnode(n)->Vsec;
            break;

          case PROJECT::kUVBOT:
          {
            double H  = nd->v.S - nd->z;
            double U  = nd->v.U;
            double V  = nd->v.V;
            double Us = sqrt( U*U + V*V );
            if( Us > 1.0e-9 )
            {
              x[npc] = nd->v.U - nd->Vsec * V/Us;
              y[npc] = nd->v.V + nd->Vsec * U/Us;
            }
            else
            {
              x[npc] = 0.0;
              y[npc] = 0.0;
            }
          }
            break;

          case PROJECT::kDZ:
            x[npc] = region->Getnode(n)->dz;
            break;

          case PROJECT::kRE:
          {
            double H  = nd->v.S - nd->z;
            double U  = nd->v.U;
            double V  = nd->v.V;
            double Us = sqrt( U*U + V*V );
            x[npc] = 4.0 * Us * H / project->vk;
          }
            break;

          case PROJECT::kFR:
          {
            double H  = nd->v.S - nd->z;
            double U  = nd->v.U;
            double V  = nd->v.V;
            double Us = sqrt( U*U + V*V );
            if( H <= 0.0 ) x[npc] = 0.0;
            else           x[npc] = Us / sqrt( project->g * H );
          }
            break;

          case PROJECT::kPE:
          {
            int    no = nd->Getno();
            double H  = nd->v.S - nd->z;
            double U  = nd->v.U;
            double V  = nd->v.V;
            if( nd->vt > 0.0 ) x[npc] = sqrt( U*Lx[no]*U*Lx[no] + V*Ly[no]*V*Ly[no] ) / nd->vt;
            else               x[npc] = 0.0;
          }
            break;

          case PROJECT::kCU:
          {
            int    no = nd->Getno();
            double H  = nd->v.S - nd->z;
            double U  = nd->v.U;
            double V  = nd->v.V;
            double dt = project->timeint.deltaTime.Getsec();
            x[npc] = sqrt(U/Lx[no]*U/Lx[no]+V/Ly[no]*V/Ly[no])*dt;
          }
            break;

          case PROJECT::kPHI:
            x[npc] = phi[region->Getnode(n)->Getno()];
            break;

          case PROJECT::kROT:
          {
            double H  = nd->v.S - nd->z;
            if( H < 0.0 ) H = 0.0;
            double U  = nd->v.U;
            double V  = nd->v.V;
            double Us = sqrt( U*U + V*V );
            if( Us > 0.0 ) x[npc] = rot[region->Getnode(n)->Getno()] * H / Us;
            else           x[npc] = 0.0;
          }
            break;

          case PROJECT::kCURV:
            x[npc] = curv[region->Getnode(n)->Getno()];
            break;
        }

        npc++;
      }
    }

    //    fseek( id, 0L, SEEK_END );
    //    fseek( id, 0L, SEEK_CUR );

    if( project->valist[tmser->voutlist[i]].dim == 1 )
    {
      fwrite( x, sizeof(double), header.np, id );
    }
    else if( project->valist[tmser->voutlist[i]].dim == 2 )
    {
      fwrite( x, sizeof(double), header.np, id );
      fwrite( y, sizeof(double), header.np, id );
    }
  }

  fclose( id );

  // detach temporary used arrays
  MEMORY::memo.Detach( x );
  MEMORY::memo.Detach( y );
}

//////////////////////////////////////////////////////////////////////////////////////////

void MODEL::AttachOutput( PROJECT* project, int val )
{

  if( !project->sed.Getinit()  &&  (   val == PROJECT::kSEDUV
                                       || val == PROJECT::kLS
                                       || val == PROJECT::kQBE
                                       || val == PROJECT::kDZDS
                                       || val == PROJECT::kDZMX
                                       || val == PROJECT::kDHDS ) )
  {
    project->sed.Initialize( project );
  }

  // -----------------------------------------------------------------------------------

  if( !phi  &&  val == PROJECT::kPHI )
  {
    phi = project->M2D->Phi2D();
  }

  // -----------------------------------------------------------------------------------

  if( !rot  &&  val == PROJECT::kROT )
  {
    rot = project->M2D->Rot2D();
  }

  // -----------------------------------------------------------------------------------

  if( !curv  &&  val == PROJECT::kCURV )
  {
    curv = project->M2D->Curv2D();
  }

  // -----------------------------------------------------------------------------------

  if( !Lx  &&
      (    val == PROJECT::kCU
           || val == PROJECT::kPE
           || val == PROJECT::kKINER ) )
  {
    Lx = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    Ly = (double*) MEMORY::memo.Array_nd( region->Getnp() );

    int* cnt = (int*) MEMORY::memo.Array_nd( region->Getnp() );

    for( int n=0; n<region->Getnp(); n++ )
    {
      Lx[n]  = 0.0;
      Ly[n]  = 0.0;
      cnt[n] = 0;
    }

    for( int e=0; e<region->Getne(); e++ )
    {
      ELEM* el = region->Getelem(e);

      if( isFS( el->flag, ELEM::kRegion ) )
      {
        int ncn = el->Getncn();
        int nnd = el->Getnnd();

        double xmax = el->nd[0]->x;
        double xmin = el->nd[0]->x;
        double ymax = el->nd[0]->y;
        double ymin = el->nd[0]->y;

        for( int i=1; i<ncn; i++ )
        {
          double x = el->nd[i]->x;
          double y = el->nd[i]->y;

          if( x > xmax )  xmax = x;
          if( x < xmin )  xmin = x;
          if( y > ymax )  ymax = y;
          if( y < ymin )  ymin = y;
        }

        for( int i=0; i<nnd; i++ )
        {
          int no = el->nd[i]->Getno();
          cnt[no]++;

          Lx[no] += xmax - xmin;
          Ly[no] += ymax - ymin;
        }
      }
    }

    for( int n=0; n<region->Getnp(); n++ )
    {
      if( cnt[n] > 0 )
      {
        Lx[n] /= cnt[n];
        Ly[n] /= cnt[n];
      }
    }

    MEMORY::memo.Detach( cnt );
  }

  // -----------------------------------------------------------------------------------

  if( !rc  &&  (   val == PROJECT::kRC
                   || val == PROJECT::kD90
                   || val == PROJECT::kD50) )
  {
    rc  = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    d90 = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    d50 = (double*) MEMORY::memo.Array_nd( region->Getnp() );

    int* cnt = (int*) MEMORY::memo.Array_nd( region->Getnp() );

    for( int n=0; n<region->Getnp(); n++ )
    {
      rc[n]  = 0.0;
      d90[n] = 0.0;
      d50[n] = 0.0;
      cnt[n] = 0;
    }

    for( int e=0; e<region->Getne(); e++ )
    {
      ELEM* el = region->Getelem(e);

      if( isFS( el->flag, ELEM::kRegion ) )
      {
        TYPE* type = TYPE::Getid( el->type );

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          int no = el->nd[i]->Getno();

          rc[no]  += type->rcoef;
          d90[no] += type->d90;
          d50[no] += type->d50;

          cnt[no]++;
        }
      }
    }

    for( int n=0; n<region->Getnp(); n++ )
    {
      if( cnt[n] > 0 )
      {
        rc[n]  /= cnt[n];
        d90[n] /= cnt[n];
        d50[n] /= cnt[n];
      }
    }

    MEMORY::memo.Detach( cnt );
  }

  if( !Hr  &&  (   val == PROJECT::kHR
                   || val == PROJECT::kHD
                   || val == PROJECT::kKD ) )
  {
    Hr = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    Hd = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    kd = (double*) MEMORY::memo.Array_nd( region->Getnp() );

    int* cnt = (int*) MEMORY::memo.Array_nd( region->Getnp() );

    for( int n=0; n<region->Getnp(); n++ )
    {
      Hr[n]  = 0.0;
      Hd[n]  = 0.0;
      kd[n]  = 0.0;
      cnt[n] = 0;
    }

    for( int e=0; e<region->Getne(); e++ )
    {
      ELEM* el = region->Getelem(e);

      if( isFS( el->flag, ELEM::kRegion ) )
      {
        TYPE* type = TYPE::Getid( el->type );

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          int no = el->nd[i]->Getno();

          double H = el->nd[i]->v.S - el->nd[i]->z;
          double U = el->nd[i]->v.U;
          double V = el->nd[i]->v.V;

          double Us = sqrt( U*U + V*V );

          if( type->form > 0  &&  H > 0.0 )
          {
            double Hde, Lde, Hre, Lre;

            if( type->kslaw == 2 )
            {
              type->Dune( Us, H, project->kappa, project->vk, project->g, project->rho,
                          project->sed.rhob, project->sed.d50, project->sed.d90,
                          &Hde, &Lde, &Hre, &Lre );
            }
            else
            {
              type->Dune( Us, H, project->kappa, project->vk, project->g, project->rho,
                          project->sed.rhob, type->d50, type->d90,
                          &Hde, &Lde, &Hre, &Lre );
            }

            Hr[no] += Hre;
            Hd[no] += Hde;

            if( Lde > 0.01 )
            {
              kd[no] += type->duneCoef * Hde * ( 1.0 - exp(-25.0*Hde/Lde) );

              //              if( Lre > 0.01 )  kd[no] += 20.0 * 0.7 * Hre * Hre / Lre;
              if( Lre > 0.01 )  kd[no] += type->duneCoef * Hre * ( 1.0 - exp(-25.0*Hre/Lre) );
            }
            else if( Lre > 0.01 )
            {
              //              kd[no] += 20.0 * Hre * Hre / Lre;
              kd[no] += type->duneCoef * Hre * ( 1.0 - exp(-25.0*Hre/Lre) );
            }

            cnt[no]++;
          }
        }
      }
    }

    for( int n=0; n<region->Getnp(); n++ )
    {
      if( cnt[n] > 0 )
      {
        Hr[n] /= cnt[n];
        Hd[n] /= cnt[n];
        kd[n] /= cnt[n];
      }
    }

    MEMORY::memo.Detach( cnt );
  }

  if( !hp  &&  (   val == PROJECT::kHP
                   || val == PROJECT::kDP
                   || val == PROJECT::kSP ) )
  {
    hp = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    dp = (double*) MEMORY::memo.Array_nd( region->Getnp() );
    sp = (double*) MEMORY::memo.Array_nd( region->Getnp() );

    int* cnt = (int*) MEMORY::memo.Array_nd( region->Getnp() );

    for( int n=0; n<region->Getnp(); n++ )
    {
      hp[n]  = 0.0;
      dp[n]  = 0.0;
      sp[n]  = 0.0;
      cnt[n] = 0;
    }

    for( int e=0; e<region->Getne(); e++ )
    {
      ELEM* el = region->Getelem(e);

      if( isFS( el->flag, ELEM::kRegion ) )
      {
        TYPE* type = TYPE::Getid( el->type );

        int nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
        {
          int no = el->nd[i]->Getno();
          cnt[no]++;

          hp[no] += type->hp;

          dp[no] += type->dp[0];
          sp[no] += type->sp[0];
        }
      }
    }

    for( int n=0; n<region->Getnp(); n++ )
    {
      if( cnt[n] > 0 )
      {
        hp[n] /= cnt[n];
        dp[n] /= cnt[n];
        sp[n] /= cnt[n];
      }
    }

    MEMORY::memo.Detach( cnt );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void MODEL::DetachOutput( PROJECT *project )
{
  if( phi )  MEMORY::memo.Detach( phi );
  if( rot )  MEMORY::memo.Detach( rot );
  if( curv ) MEMORY::memo.Detach( curv );
  if( Lx )   MEMORY::memo.Detach( Lx );
  if( Ly )   MEMORY::memo.Detach( Ly );
  if( rc )   MEMORY::memo.Detach( rc );
  if( d90 )  MEMORY::memo.Detach( d90 );
  if( d50 )  MEMORY::memo.Detach( d50 );
  if( Hr )   MEMORY::memo.Detach( Hr );
  if( Hd )   MEMORY::memo.Detach( Hd );
  if( kd )   MEMORY::memo.Detach( kd );
  if( hp )   MEMORY::memo.Detach( hp );
  if( dp )   MEMORY::memo.Detach( dp );
  if( sp )   MEMORY::memo.Detach( sp );

  project->sed.Detach();

  // initialize variables
  phi   = NULL;
  rot   = NULL;
  curv  = NULL;
  Lx    = NULL;
  Ly    = NULL;
  rc    = NULL;
  d90   = NULL;
  d50   = NULL;
  Hr    = NULL;
  Hd    = NULL;
  kd    = NULL;
  hp    = NULL;
  dp    = NULL;
  sp    = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////

void MODEL::UCDOutput( FILE*    id,
                       PROJECT* project,
                       GRID*    grid,
                       double*  phi,
                       double*  rot,
                       double*  curv,
                       double*  Lx,
                       double*  Ly,
                       double*  rc,
                       double*  d90,
                       double*  d50,
                       double*  Hr,
                       double*  Hd,
                       double*  kd,
                       double*  hp,
                       double*  dp,
                       double*  sp,
                       double*  qbc,
                       double*  PLs,
                       double*  sx,
                       double*  sy,
                       double*  dzds,
                       double*  dzmx,
                       double*  dhds )
{
  // write information on values and dimension -------------------------------------------

  fprintf( id, "%d", project->vpcomp );

  for( int i=0; i<project->vpcomp; i++ )
  {
    int vl = project->vpoutlist[i];
    fprintf( id, " %d", project->valist[vl].vec );
  }

  fprintf( id, "\n" );

  for( int i=0; i<project->vpcomp; i++ )
  {
    int vl = project->vpoutlist[i];
    fprintf( id, "%s, %s\n", project->valist[vl].name, project->valist[vl].unit );
  }


  // write nodal values ------------------------------------------------------------------

  for( int n=0; n<grid->Getnp(); n++ )
  {
    NODE* nd = grid->Getnode(n);
    int   no = nd->Getno();

    double U  = nd->v.U;
    double V  = nd->v.V;
    double H  = nd->v.S - nd->z;

    double Us = sqrt( U*U + V*V );

    if( H < 0.0 )  H = 0.0;

    if( nd->mark )
    {
      double r;

      fprintf( id, "%6d", nd->Getname() );

      for( int i=0; i<project->vpcomp; i++ )
      {
        switch( project->vpoutlist[i] )
        {
          case PROJECT::kUV:
            fprintf( id, " %14.6le %14.6le %14.6le", U, V, 0.0 );
            break;

          case PROJECT::kS:
            fprintf( id, " %14.6le", nd->v.S );
            break;

          case PROJECT::kK:
            fprintf( id, " %14.6le", nd->v.K );
            break;

          case PROJECT::kD:
            fprintf( id, " %14.6le", nd->v.D );
            break;

          case PROJECT::kC:
            fprintf( id, " %14.6le", nd->v.C );
            break;

          case PROJECT::kCB:
            if( Us > 0.0 )  fprintf( id, " %14.6le", nd->v.Qb / Us );
            else            fprintf( id, " %14.6le", 0.0 );
            break;

          case PROJECT::kQB:
            fprintf( id, " %14.6le", nd->v.Qb );
            break;

          case PROJECT::kQBE:
            fprintf( id, " %14.6le", qbc[no] );
            break;

          case PROJECT::kLS:
            if( PLs )  fprintf( id, " %14.6le", PLs[no] );
            else       fprintf( id, " %14.6le", 0.0 );
            break;

          case PROJECT::kSEDUV:
            fprintf( id, " %14.6le %14.6le %14.6le", nd->v.Qb*sx[no], nd->v.Qb*sy[no], 0.0 );
            break;

          case PROJECT::kDUVDT:
            fprintf( id, " %14.6le %14.6le %14.6le", nd->v.dUdt, nd->v.dVdt, 0.0 );
            break;

          case PROJECT::kDSDT:
            fprintf( id, " %14.6le", nd->v.dSdt );
            break;

          case PROJECT::kZ:
            fprintf( id, " %14.6le", nd->z );
            break;

          case PROJECT::kH:
            fprintf( id, " %14.6le", H );
            break;

          case PROJECT::kUS:
            fprintf( id, " %14.6le", Us );
            break;

          case PROJECT::kCF:
            fprintf( id, " %14.6le", nd->cf );
            break;

          case PROJECT::kRC:
            fprintf( id, " %14.6le", rc[no] );
            break;

          case PROJECT::kD90:
            fprintf( id, " %14.6le", d90[no] );
            break;

          case PROJECT::kD50:
            fprintf( id, " %14.6le", d50[no] );
            break;

          case PROJECT::kHR:
            fprintf( id, " %14.6le", Hr[no] );
            break;

          case PROJECT::kHD:
            fprintf( id, " %14.6le", Hd[no] );
            break;

          case PROJECT::kKD:
            fprintf( id, " %14.6le", kd[no] );
            break;

          case PROJECT::kHP:
            fprintf( id, " %14.6le", hp[no] );
            break;

          case PROJECT::kDP:
            fprintf( id, " %14.6le", dp[no] );
            break;

          case PROJECT::kSP:
            fprintf( id, " %14.6le", sp[no] );
            break;

          case PROJECT::kTAU:
          {
            double Utau = project->sed.GetUtau( Us, H, 0.0, project );
            fprintf( id, " %14.6le", project->rho * Utau * Utau );
          }
            break;

          case PROJECT::kMAN:
            if( H > 0.0 ) fprintf( id, " %14.6le", sqrt(nd->cf*pow(H,1.0/3.0)/project->g) );
            else          fprintf( id, " %14.6le", 0.0 );
            break;

          case PROJECT::kVT:
            fprintf( id, " %14.6le", nd->vt );
            break;

          case PROJECT::kEst:
            if( nd->v.D > 0.0  &&  Us > 0.0 )
            {
              double Ust = sqrt( nd->cf ) * Us;
              double vt  = project->KD.cm * nd->v.K * nd->v.K / nd->v.D;
              double est = vt / Ust / H;

              fprintf( id, " %14.6le", est );
            }
            else
            {
              fprintf( id, " %14.6le", 0.0 );
            }
            break;

          case PROJECT::kEXX:
            r = sqrt( 4.0*nd->exy*nd->exy + (nd->exx-nd->eyy)*(nd->exx-nd->eyy) );
            fprintf( id, " %14.6le", 0.5*(nd->exx + nd->eyy + r) * nd->vt );
            break;

          case PROJECT::kEYY:
            r = sqrt( 4.0*nd->exy*nd->exy + (nd->exx-nd->eyy)*(nd->exx-nd->eyy) );
            fprintf( id, " %14.6le", 0.5*(nd->exx + nd->eyy - r) * nd->vt );
            break;

          case PROJECT::kDUU:
            fprintf( id, " %14.6le", H * Us * Us * nd->Dxx );
            break;

          case PROJECT::kDUV:
            fprintf( id, " %14.6le", H * Us * Us * nd->Dxy );
            break;

          case PROJECT::kDVV:
            fprintf( id, " %14.6le", H * Us * Us * nd->Dyy );
            break;

          case PROJECT::kVSEC:
            fprintf( id, " %14.6le", nd->Vsec );
            break;

          case PROJECT::kUVBOT:
            if( Us > 1.0e-9 )
            {
              double Ub = nd->v.U - nd->Vsec * V/Us;
              double Vb = nd->v.V + nd->Vsec * U/Us;
              fprintf( id, " %14.6le %14.6le %14.6le", Ub, Vb, 0.0 );
            }
            else
            {
              fprintf( id, " %14.6le %14.6le %14.6le", 0.0, 0.0, 0.0 );
            }

            break;

          case PROJECT::kDZ:
            fprintf( id, " %14.6le", nd->dz );
            break;

          case PROJECT::kRE:
            fprintf( id, " %14.6le", 4.0 * Us * H / project->vk );
            break;

          case PROJECT::kFR:
            if( H <= 0.0 )  fprintf( id, " %14.6le", 0.0 );
            else            fprintf( id, " %14.6le", Us/sqrt(project->g*H) );
            break;

          case PROJECT::kPE:
            if( nd->vt > 0.0 )
            {
              double Pe = sqrt( U*Lx[no]*U*Lx[no] + V*Ly[no]*V*Ly[no] ) / nd->vt;
              fprintf( id, " %14.6le", Pe );
            }
            else
            {
              fprintf( id, " %14.6le", 0.0 );
            }
            break;

          case PROJECT::kCU:
          {
            double dt = project->timeint.deltaTime.Getsec();
            fprintf( id, " %14.6le", sqrt(U/Lx[no]*U/Lx[no]+V/Ly[no]*V/Ly[no])*dt );
          }
            break;

          case PROJECT::kPHI:
            fprintf( id, " %14.6le", phi[no] );
            break;

          case PROJECT::kROT:
            if( Us > 0.0 )  fprintf( id, " %14.6le", rot[no]*H/Us );
            else            fprintf( id, " %14.6le", 0.0 );
            break;

          case PROJECT::kCURV:
            fprintf( id, " %14.6le", curv[no] );
            break;


            // for internal usage (debugging) ...
          case PROJECT::kDZDS:
            if( dzds )  fprintf( id, " %14.6le", dzds[no] );
            break;

          case PROJECT::kDZMX:
            if( dzmx )  fprintf( id, " %14.6le", dzmx[no] );
            break;

          case PROJECT::kDHDS:
            if( dhds )  fprintf( id, " %14.6le", dhds[no] );
            break;

            // output of statistics
          case PROJECT::kMEANUV:
            fprintf( id, " %14.6le", project->statist->GetMeanU(no) );
            fprintf( id, " %14.6le", project->statist->GetMeanV(no) );
            fprintf( id, " %14.6le", 0.0 );
            break;

          case PROJECT::kMEANS:
            fprintf( id, " %14.6le", project->statist->GetMeanS(no) );
            break;

          case PROJECT::kMEANUS:
            fprintf( id, " %14.6le", project->statist->GetMeanUs(no) );
            break;

          case PROJECT::kMEANH:
            fprintf( id, " %14.6le", project->statist->GetMeanH(no) );
            break;

          case PROJECT::kMEANVT:
            fprintf( id, " %14.6le", project->statist->GetMeanVt(no) );
            break;

          case PROJECT::kVARU:
            fprintf( id, " %14.6le", project->statist->GetVarU(no) );
            break;

          case PROJECT::kVARV:
            fprintf( id, " %14.6le", project->statist->GetVarV(no) );
            break;

          case PROJECT::kVARUV:
            fprintf( id, " %14.6le", project->statist->GetVarUV(no) );
            break;

          case PROJECT::kKINE:
            fprintf( id, " %14.6le", project->statist->GetKinE(no) );
            break;

          case PROJECT::kSDEVH:
            fprintf( id, " %14.6le", project->statist->GetSdevH(no) );
            break;

          case PROJECT::kVARVT:
            fprintf( id, " %14.6le", project->statist->GetVarVt(no) );
            break;

          case PROJECT::kKINER:
          {
            //            //double L2 = 0.094 * 0.094 * Lx[no] * Ly[no];
            //            double L2 = sqrt( project->KD.cm*project->KD.cd ) * TYPE::deflt.lm * TYPE::deflt.lm * Lx[no] * Ly[no]; // deflt.lm ggf. anpassen
            //            double kr = 0.0;
            //            if( L2 > 1.0e-9 )  kr = project->statist->GetVtVt(no) / L2;

            double kr = 0.0;
            kr = -1.5 * nd->uu;

            fprintf( id, " %14.6le", kr );
            break;
          }

          case PROJECT::kFLDRATE:
            fprintf( id, " %14.6le", project->statist->GetFldRate(no) );
            break;
          case PROJECT::kKINRATIO:
          {
            // kinEr
            //            double L2 = 0.094 * 0.094 * Lx[no] * Ly[no];
            double L2 = sqrt( project->KD.cm*project->KD.cd ) * TYPE::deflt.lm * TYPE::deflt.lm * Lx[no] * Ly[no]; // deflt.lm ggf. anpassen
            double kr = 0.0;
            if( L2 > 1.0e-9 )  kr = project->statist->GetVtVt(no) / L2;
            // kinE
            double kE = project->statist->GetKinE(no);

            double kinRatio = 0.0;
            if ( (kE + kr ) > 1.0e-9 ) kinRatio = kE / (kE + kr );

            fprintf( id, " %14.6le", kinRatio );
            break;
          }

          case PROJECT::kMAXUV:
            fprintf( id, " %14.6le", project->statist->GetMaxU(no) );
            fprintf( id, " %14.6le", project->statist->GetMaxV(no) );
            fprintf( id, " %14.6le", 0.0);
            break;

          case PROJECT::kMINUV:
            fprintf( id, " %14.6le", project->statist->GetMinU(no) );
            fprintf( id, " %14.6le", project->statist->GetMinV(no) );
            fprintf( id, " %14.6le", 0.0);
            break;

          case PROJECT::kMAXUS:
            fprintf( id, " %14.6le", project->statist->GetMaxUs(no) );
            break;

          case PROJECT::kMINUS:
            fprintf( id, " %14.6le", project->statist->GetMinUs(no) );
            break;

          case PROJECT::kMAXTAU:
            fprintf( id, " %14.6le", project->statist->GetMaxTau(no) );
            break;

          case PROJECT::kMAXU:
            fprintf( id, " %14.6le", project->statist->GetMaxU_scalar(no) );
            break;
          case PROJECT::kMINU:
            fprintf( id, " %14.6le", project->statist->GetMinU_scalar(no) );
            break;
          case PROJECT::kMAXV:
            fprintf( id, " %14.6le", project->statist->GetMaxV_scalar(no) );
            break;
          case PROJECT::kMINV:
            fprintf( id, " %14.6le", project->statist->GetMinV_scalar(no) );
            break;
        }
      }

      fprintf( id, "\n" );
    }
  }

  // write information on values and dimension -------------------------------------------

  if( project->vecomp >= 0 )
  {
    fprintf( id, "%d", project->vecomp );

    for( int i=0; i<project->vecomp; i++ )
    {
      int vl = project->veoutlist[i];
      fprintf( id, " %d", project->valist[vl].vec );
    }

    fprintf( id, "\n" );

    for( int i=0; i<project->vecomp; i++ )
    {
      int vl = project->veoutlist[i];
      fprintf( id, "%s, %s\n", project->valist[vl].name, project->valist[vl].unit );
    }


    // write element values --------------------------------------------------------------

    for( int e=0; e<grid->Getne(); e++ )
    {
      ELEM* el = grid->Getelem(e);

      if( el->mark )
      {
        TYPE* type = NULL;
        if( el->shape != kLine ) type = TYPE::Getid( el->type );

        fprintf( id, "%6d", el->Getname() );

        for( int i=0; i<project->vecomp; i++ )
        {
          switch( project->veoutlist[i] )
          {
            case PROJECT::kS:
              fprintf( id, " %14.6le", el->P );
              break;

            case PROJECT::kUV:
              fprintf( id, " %14.6le %14.6le", el->U, el->V );
              break;

            case PROJECT::kDZ:
              fprintf( id, " %14.6le", el->dz );
              break;

            case PROJECT::kRC:
              if( type ) fprintf( id, " %14.6le", type->rcoef );
              else       fprintf( id, " %14.6le", -1.0 );
              break;

            case PROJECT::kD90:
              if( type ) fprintf( id, " %14.6le", type->d90 );
              else       fprintf( id, " %14.6le", -1.0 );
              break;

            case PROJECT::kD50:
              if( type ) fprintf( id, " %14.6le", type->d50 );
              else       fprintf( id, " %14.6le", -1.0 );
              break;

            case PROJECT::kHP:
              if( type ) fprintf( id, " %14.6le", type->hp );
              else       fprintf( id, " %14.6le", -1.0 );
              break;

            case PROJECT::kDP:
              if( type ) fprintf( id, " %14.6le", type->dp[0] );
              else       fprintf( id, " %14.6le", -1.0 );
              break;

            case PROJECT::kSP:
              if( type ) fprintf( id, " %14.6le", type->sp[0] );
              else       fprintf( id, " %14.6le", -1.0 );
              break;
          }
        }

        fprintf( id, "\n" );
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////
// assemble vector "vec[]" across subdomains
//
//void MODEL::Mpi_assemble( double* vec, PROJECT* project )
//{
//# ifdef _MPI_
//  if( project->subdom.npr > 1 )
//  {
//    SUBDOM* subdom = &project->subdom;
//    INFACE* inface = subdom->inface;
//
//    int rgnp = region->Getnp();
//
//    // copy vector "vec[]" to a temporary array ------------------------------------------
//    double* tmp = (double*) MEMORY::memo.Array_nd( rgnp );
//    memcpy( tmp, vec, rgnp*sizeof(double) );
//
//    // loop on all interfaces 0 <= s < npr: exchange vector data -------------------------
//    for( int s=0; s<subdom->npr; s++ )
//    {
//      MPI_Status status;
//
//      int npinf = inface[s].np;             // number of nodes on interface to domain s
//
//      if( npinf > 0 )
//      {
//        for( int n=0; n<npinf; n++ )
//        {
//          NODE* nd = inface[s].node[n];     // pointer to the interface node n
//          int   no = nd->Getno();           // node number local in this domain
//
//          inface[s].send[n] = tmp[no];      // copy the vector element to the send array
//        }
//
//        MPI_Sendrecv( inface[s].send, npinf, MPI_DOUBLE, s, 1,
//                      inface[s].recv, npinf, MPI_DOUBLE, s, 1,
//                      MPI_COMM_WORLD, &status );
//
//        for( int n=0; n<npinf; n++ )
//        {
//          NODE* nd = inface[s].node[n];
//          int   no = nd->Getno();
//
//          vec[no] += inface[s].recv[n];
//        }
//      }
//    }
//
//    MEMORY::memo.Detach( tmp );
//  }
//# endif
//}
//
//
////////////////////////////////////////////////////////////////////////////////////////////
//// average vector "vec[]" across subdomains
//
//void MODEL::Mpi_average( double* vec, PROJECT* project )
//{
//# ifdef _MPI_
//  if( project->subdom.npr > 1 )
//  {
//    SUBDOM* subdom = &project->subdom;
//    INFACE* inface = subdom->inface;
//
//    int rgnp = region->Getnp();
//
//    // copy vector "vec[]" to a temporary array ------------------------------------------
//    int*    cnt = (int*)    MEMORY::memo.Array_nd( rgnp );
//    double* tmp = (double*) MEMORY::memo.Array_nd( rgnp );
//
//    memcpy( tmp, vec, rgnp*sizeof(double) );
//
//    for( int n=0; n<rgnp; n++ ) cnt[n] = 1;
//
//    // loop on all interfaces: exchange vector data --------------------------------------
//    for( int s=0; s<subdom->npr; s++ )
//    {
//      MPI_Status status;
//
//      int np = inface[s].np;
//
//      if( np > 0 )
//      {
//        for( int n=0; n<np; n++ )
//        {
//          NODE* nd = inface[s].node[n];
//          inface[s].send[n] = tmp[nd->Getno()];
//        }
//
//        MPI_Sendrecv( inface[s].send, np, MPI_DOUBLE, s, 1,
//                      inface[s].recv, np, MPI_DOUBLE, s, 1,
//                      MPI_COMM_WORLD, &status );
//
//        for( int n=0; n<np; n++ )
//        {
//          NODE* nd = inface[s].node[n];
//          vec[nd->Getno()] += inface[s].recv[n];
//          cnt[nd->Getno()]++;
//        }
//      }
//    }
//
//    // average vector --------------------------------------------------------------------
//    for( int n=0; n<rgnp; n++ )
//    {
//      if( cnt[n] > 1 )  vec[n] /= cnt[n];
//    }
//
//    // detach memory ---------------------------------------------------------------------
//    MEMORY::memo.Detach( cnt );
//    MEMORY::memo.Detach( tmp );
//  }
//# endif
//}
