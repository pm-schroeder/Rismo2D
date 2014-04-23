// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class STATIST
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
#include "Node.h"
#include "Elem.h"
#include "Shape.h"
#include "Type.h"
#include "Grid.h"
#include "Model.h"
#include "Project.h"
#include "Subdom.h"

#include "Statist.h"


STATIST::STATIST()
{
}


void STATIST::Init( int np )
{
  this->n    = 0;

  this->np   = np;

  this->U    = new double [np];
  this->V    = new double [np];
  this->S    = new double [np];
  this->H    = new double [np];
  this->Vt   = new double [np];

  this->UU   = new double [np];
  this->UV   = new double [np];
  this->VV   = new double [np];
  this->HH   = new double [np];

  this->VtVt = new double [np];

  this->nwet = new int [np];

  if(     !this->U   ||  !this->V    ||  !this->S   ||  !this->H
      ||  !this->Vt  ||  !this->UU   ||  !this->VV  ||  !this->UV
      ||  !this->HH  ||  !this->VtVt ||  !this->nwet )
    REPORT::rpt.Error( "can not allocate memory - STATIST::STATIST (1)" );


  for( int i=0; i<np; i++ )
  {
    this->nwet[i] = 0;
    this->U[i]    = 0.0;
    this->V[i]    = 0.0;
    this->S[i]    = 0.0;
    this->H[i]    = 0.0;
    this->Vt[i]   = 0.0;

    this->UU[i]   = 0.0;
    this->VV[i]   = 0.0;
    this->UV[i]   = 0.0;
    this->HH[i]   = 0.0;
    this->VtVt[i] = 0.0;
  }
}


STATIST::~STATIST()
{
  delete[] U;
  delete[] V;
  delete[] S;
  delete[] H;
  delete[] Vt;

  delete[] UU;
  delete[] UV;
  delete[] VV;
  delete[] HH;
  delete[] VtVt;

  delete[] nwet;
}


double STATIST::GetMeanU( int no )
{
  if( nwet[no] > 0 ) return U[no] / nwet[no];
  else               return 0.0;
}

double STATIST::GetMeanV( int no )
{
  if( nwet[no] > 0 ) return V[no] / nwet[no];
  else               return 0.0;
}

double STATIST::GetMeanS( int no )
{
  if( nwet[no] > 0 ) return S[no] / nwet[no];
  else               return 0.0;
}

double STATIST::GetMeanH( int no )
{
  if( nwet[no] > 0 ) return H[no] / nwet[no];
  else               return 0.0;
}

double STATIST::GetMeanUs( int no )
{
  if( nwet[no] > 0 )
  {
    double mU = U[no] / nwet[no];
    double mV = V[no] / nwet[no];

    return sqrt( mU*mU + mV*mV );
  }
  else
  {
    return 0.0;
  }
}

double STATIST::GetMeanVt( int no )
{
  if( nwet[no] > 0 ) return Vt[no] / nwet[no];
  else               return 0.0;
}

double STATIST::GetVarU( int no )
{
  if( nwet[no] > 1 )  return (UU[no] - U[no]*U[no]/nwet[no]) / (nwet[no]-1);
  else                return 0.0;
}

double STATIST::GetVarUV( int no )
{
  if( nwet[no] > 1 )  return (UV[no] - U[no]*V[no]/nwet[no]) / (nwet[no]-1);
  else                return 0.0;
}

double STATIST::GetVarV( int no )
{
  if( nwet[no] > 1 )  return (VV[no] - V[no]*V[no]/nwet[no]) / (nwet[no]-1);
  else                return 0.0;
}

double STATIST::GetKinE( int no )
{
  if( nwet[no] > 1 )  return 0.5 * ( GetVarU(no) + GetVarV(no) );
  else                return 0.0;
}

double STATIST::GetSdevH( int no )
{
  if( nwet[no] > 1 )
  {
    double hh = (HH[no] - H[no]*H[no]/nwet[no]) / (nwet[no]-1);
    
    if( hh <= 0.0 )  return 0.0;
    else             return sqrt( hh );
  }
  else
  {
    return 0.0;
  }
}

double STATIST::GetVarVt( int no )
{
  if( nwet[no] > 1 )  return (VtVt[no] - Vt[no]*Vt[no]/nwet[no]) / (nwet[no]-1);
  else                return 0.0;
}

double STATIST::GetVtVt( int no )
{
  return VtVt[no] / nwet[no];
}

double STATIST::GetFldRate( int no )
{
  if( nwet[no] > 1 )  return double(nwet[no]) / n;
  else                return 0.0;
}


void STATIST::Read( int np, char *statisticFile, SUBDOM *subdom )
{
  if( np != this->np  ||  !statisticFile[0] ) return;

  ////////////////////////////////////////////////////////////////////////////////////////
  // open the file and read header lines

  char  text[200];
  char* textLine;

  ASCIIFILE* file = new ASCIIFILE( statisticFile, "r" );

  if( !file || !file->getid() )
    REPORT::rpt.Error( kOpenFileFault, "%s %s (STATIST::Read - 1)",
                       "can not open statistic file", statisticFile );

  textLine = file->next();
  sscanf( textLine, "# %d", &this->n );

  int npsta;
  textLine = file->nextLine();
  sscanf( textLine, "%d", &npsta );

  if(     (subdom->npr == 1  &&  npsta != np)
      ||  (subdom->npr  > 1  &&  npsta != subdom->np) )
    REPORT::rpt.Error( kOpenFileFault, "%s %s (STATIST::Read - 2)",
                       "wrong number of nodes in statistic file", statisticFile );


  ////////////////////////////////////////////////////////////////////////////////////////
  // read the file

  for( int i=0; i<npsta; i++ )
  {
    int name;

    textLine = file->nextLine();
    sscanf( textLine, "%d", &name );

    if( name <= 0  ||  name > npsta )
    {
      sprintf( text, "%s %d - (STATIST::Read - 3)",
               "node numbers must be between 1 and", npsta );
      REPORT::rpt.Error( text );
    }

    int no = -1;

    if( subdom->npr > 1 )
    {
      NODE* nd = subdom->node[name - 1];
      if( nd )  no = nd->Getno();
    }
    else
    {
      no = name - 1;
    }

    if( no >= 0 )
    {
      sscanf( textLine, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                        &name, &this->nwet[no],
                               &this->U[no],
                               &this->V[no],
                               &this->H[no],
                               &this->Vt[no],
                               &this->UU[no],
                               &this->VV[no],
                               &this->UV[no],
                               &this->HH[no],
                               &this->VtVt[no] );
    }
  }

  delete file;
}


void STATIST::Write( MODEL   *model,
                     int      release,
                     char    *statisticFile,
                     char    *regionFile,
                     int      timeStep,
                     PROJECT *project )
{
  GRID *rg = model->region;

  if( rg->Getnp() != this->np )  return;


  // ---------------------------------------------------------------------------------------
  // write statistic file

  if( statisticFile[0] )
  {
    char filename[500];
    project->ReplaceMacro( statisticFile, filename, timeStep, 0 );

    FILE *id = fopen( filename, "w" );

    if( !id )
      REPORT::rpt.Error( kOpenFileFault, "%s %s (STATIST::Write #1)",
                                         "can not open statistic file", filename );

    fprintf( id, "# %d   Release %d   nwet,U,V,H,Vt,uu,vv,uv,hh,vtvt\n", this->n, release );
    fprintf( id, "%d\n", np );

    for( int n=0; n<np; n++ )
    {
      NODE* nd = rg->Getnode(n);

      fprintf( id, "%8d %8d  %14.6le %14.6le %14.6le %14.6le",
                   nd->Getname(), this->nwet[n],
                                  this->U[n],
                                  this->V[n],
                                  this->H[n],
                                  this->Vt[n] );

      fprintf( id, " %14.6le %14.6le %14.6le %14.6le %14.6le\n",
                        this->UU[n],
                        this->VV[n],
                        this->UV[n],
                        this->HH[n],
                        this->VtVt[n] );
    }

    fclose( id );
  }


  // ---------------------------------------------------------------------------------------
  // write region file

  if( regionFile[0] )
  {
    char filename[500];
    project->ReplaceMacro( regionFile, filename, timeStep, 0 );

    FILE* id = fopen( filename, "w" );

    if( !id )
      REPORT::rpt.Error( kOpenFileFault, "%s %s (STATIST::Write #2)",
                        "can not open statistic file", filename );


    // write number of nodes, elements and node data ---------------------------------------

    fprintf( id, "%6d  %6d  11  0  0\n", rg->Getnp(), rg->Getne() );


    // write nodes -------------------------------------------------------------------------

    for( int i=0; i<rg->Getnp(); i++ )
    {
      NODE* nd = rg->Getnode(i);

      fprintf( id, "%6d  %12.6lf  %12.6lf  %12.6lf\n",
                   nd->Getname(), nd->x, nd->y, nd->zor );
    }


    // write element connectivity ----------------------------------------------------------

    for( int e=0; e<rg->Getne(); e++ )
    {
      ELEM* el = rg->Getelem(e);

      int  nnd;
      char elemShape[6];

      if( isFS(el->flag, ELEM::kRegion) )
      {
        switch( el->GetshapeID() )
        {
          case kLine:      strcpy( elemShape, "line" );   break;
          case kTriangle:  strcpy( elemShape, "tri" );    break;
          case kSquare:    strcpy( elemShape, "quad" );   break;
        }

        int mat = TYPE::Getid(el->type)->no(el->type);

        fprintf( id, "%6d  %5d  %5s", el->Getname(), mat, elemShape );

        nnd = el->Getnnd();

        for( int i=0; i<nnd; i++ )
          fprintf( id, "  %6d", el->nd[i]->Getname() );

        fprintf( id, "\n" );
      }
    }


    // write information on dimension ------------------------------------------------------

    fprintf( id, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
                 "11 3 1 1 1 1 1 1 1 1 1 1",
                 "UV, m/s",
                 "S, mNN",
                 "K, m2/s2",
                 "D, m2/s3",
                 "Vt, m2/s",
                 "uu, m2/s2",
                 "uv, m2/s2",
                 "vv, m2/s2",
                 "k, m2/s2",
                 "vt, m4/s2",
                 "h, m" );


    // write nodal values ------------------------------------------------------------------

    int n = this->n;

    if( n > 1 )
    {
      for( int i=0; i<rg->Getnp(); i++ )
      {
        NODE* nd = rg->Getnode(i);

        double U, V, H, S, Vt, K, D;
        double uu, uv, vv, h, k, vt;

        U  = GetMeanU(i);
        V  = GetMeanV(i);
        S  = GetMeanS(i);
        H  = GetMeanH(i);
        Vt = GetMeanVt(i);

        uu = GetVarU(i);
        uv = GetVarUV(i);
        vv = GetVarV(i);
        h  = GetSdevH(i);
        k  = 0.75*(uu + vv);
        vt = GetVarVt(i);

        K = nd->v.K;
        D = nd->v.D;

        fprintf( id, "%6d  %14.6le %14.6le %14.6le %14.6le %14.6le %14.6le %14.6le",
                 nd->Getname(), U, V, 0.0, S, K, D, Vt );

        fprintf( id, " %14.6le %14.6le %14.6le %14.6le %14.6le %14.6le\n",
                 uu, uv, vv, h, k, vt );
      }
    }

    else
    {
      for( int i=0; i<rg->Getnp(); i++ )
      {
        NODE* nd = rg->Getnode(i);

        double U, V, H, S, K, D, Vt;

        U  = this->U[i];
        V  = this->V[i];
        S  = this->S[i];
        H  = this->H[i];
        Vt = this->Vt[i];

        K = nd->v.K;
        D = nd->v.D;

        fprintf( id, "%6d  %14.6le %14.6le %14.6le %14.6le %14.6le %14.6le %14.6le",
                 nd->Getname(), U, V, 0.0, S, K, D, Vt );

        fprintf( id, " %14.6le %14.6le %14.6le %14.6le %14.6le %14.6le\n",
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      }
    }


    fclose( id );
  }
}


void STATIST::Sum( MODEL* model )
{
  this->n++;

  GRID* rg = model->region;

  for( int i=0; i<rg->Getnp(); i++ )
  {
    NODE* nd = rg->Getnode(i);

    if( !isFS(nd->flag,NODE::kDry) )
    {
      this->nwet[i]++;

      double H = nd->v.S - nd->z;

      this->U[i]    += nd->v.U;
      this->V[i]    += nd->v.V;
      this->S[i]    += nd->v.S;
      this->H[i]    += H;

      this->Vt[i]   += nd->vt;

      this->UU[i]   += nd->v.U * nd->v.U;
      this->VV[i]   += nd->v.V * nd->v.V;
      this->UV[i]   += nd->v.U * nd->v.V;
      this->HH[i]   += H * H;
      this->VtVt[i] += nd->vt * nd->vt;
    }
  }
}

void STATIST::Reset( MODEL* model )
{
  char text[500];
  sprintf( text, "\n (STATIST::Reset)        reset sum of values in statistics");
  REPORT::rpt.Output( text, 3 );

  GRID* rg = model->region;

  this->n = 0;

  for( int i=0; i<rg->Getnp(); i++ )
  {
    this->nwet[i] = 0;
    this->U[i]    = 0.0;
    this->V[i]    = 0.0;
    this->S[i]    = 0.0;
    this->H[i]    = 0.0;
    this->Vt[i]   = 0.0;

    this->UU[i]   = 0.0;
    this->VV[i]   = 0.0;
    this->UV[i]   = 0.0;
    this->HH[i]   = 0.0;
    this->VtVt[i] = 0.0;

  }
}
