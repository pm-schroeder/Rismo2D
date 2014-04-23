// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class TYPE
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

#include "Datkey.h"
#include "Defs.h"
#include "Asciifile.h"
#include "Report.h"

#include "Type.h"

int   TYPE::nType = 0;
TYPE* TYPE::pType = NULL;
TYPE  TYPE::deflt;
int   TYPE::nUndef = 0;
int   TYPE::list_of_undefined[100];

int   TYPE::release = 0;


TYPE::TYPE()
{
  _id        =  0;
  _no        =  0;

  rtype      =  3;
  rcoef      =  0.0;
  kslaw      =  0;
  ksfact     =  1.0;

  d90        =  0.01;
  d50        =  0.002;
  form       =  0;
  duneCoef   =  0.7;

  dp[0]      =  0.0;
  sp[0]      =  0.0;

  hp         =  0.0;
  dp[1]      =  0.0;
  sp[1]      =  0.0;

  vt         =  0.0001;
  st         =  1.0;
  KDest      =  0.6;
  lm         =  0.1;

  estx       =  6.0;
  esty       =  0.6;

  betaSf     =  1.0;

  for( int i=0; i<kMaxUndef; i++ )  list_of_undefined[i] = 0;
}


TYPE::~TYPE()
{
  if( nType )
  {
    nType = 0;
    delete[] pType;
  }
}


TYPE& TYPE::operator =( const TYPE& t )
{
  rtype    = t.rtype;
  rcoef    = t.rcoef;
  kslaw    = t.kslaw;
  ksfact   = t.ksfact;

  d90      = t.d90;
  d50      = t.d50;
  form     = t.form;
  duneCoef = t.duneCoef;

  dp[0]    = t.dp[0];
  sp[0]    = t.sp[0];

  hp       = t.hp;
  dp[1]    = t.dp[1];
  sp[1]    = t.sp[1];

  vt       = t.vt;
  st       = t.st;
  estx     = t.estx;
  esty     = t.esty;
  KDest    = t.KDest;
  lm       = t.lm;

  betaSf   = t.betaSf;

  return *this;
}


//////////////////////////////////////////////////////////////////////////////////////////

void TYPE::Import( const char* filename )
{
  ASCIIFILE* tbl = new ASCIIFILE( filename, "r" );

  if( !tbl  ||  !tbl->getid() )
  {
    REPORT::rpt.Error( kOpenFileFault, "%s %s (TYPE::import - 1)",
                       "can not open roughness table", filename );
  }

  // check version of input file ---------------------------------------------------------

  char* textLine = tbl->next();
  sscanf( textLine, " $RISMO2D %d", &release );
  REPORT::rpt.Screen( 2, "\n ### release of roughness table:         %d ###\n\n", release );

       if( release >= 40000 )   Import_40000( tbl );
  else if( release >= 30900 )   Import_30900( tbl );
  else                          Import_00000( tbl );

  delete tbl;
}


void TYPE::Import_00000( ASCIIFILE* tbl )
{
  tbl->rewind();

  char text[400];

  // read element type specifications ---------------------------------------------------

  int nTypes;

  char* textLine = tbl->nextLine();
  sscanf( textLine, " %d", &nTypes );

  TYPE::Init( nTypes );


  sprintf( text, "  %d element type specifications:\n\n %s\n", TYPE::Get(),
                 "  no   bottom        dp       sp" );
  REPORT::rpt.Output( text, 5 );

  for( int i=0; i<nTypes; i++ )
  {
    int    no;
    int    idummy;
    double ddummy;

    TYPE* type = TYPE::Getid( i+1 );

    textLine = tbl->nextLine();
    sscanf( textLine, " %d %d %lf %d %lf %lf %lf",
                      &no,
                      &type->rtype, &type->rcoef,
                      &idummy,      &ddummy,
                      &type->dp[0], &type->sp[0] );
    type->Setno( no );

    sprintf( text, "  %2d   %1d %8.2le",
                   type->no(0),
                   type->rtype,  type->rcoef );
    REPORT::rpt.Output( text, 5 );

    sprintf( text, "   %7.4lf  %7.4lf\n", type->dp[0], type->sp[0] );
    REPORT::rpt.Output( text, 5 );
  }


  REPORT::rpt.Output( "\n\n  no       vt       st     KDest    lm      lat\n", 5 );

  for( int i=0; i<nTypes; i++ )
  {
    int no;

    textLine = tbl->nextLine();

    sscanf( textLine, " %d", &no );
    TYPE* type = TYPE::Getno( no );

    sscanf( textLine, " %d %lf %lf %lf %lf",
                      &no, &(type->vt), &(type->st),
                           &(type->KDest), &(type->lm) );

    sprintf( text, "  %2d   %9.6lf  %6.3lf  %6.3lf  %6.3lf\n",
                   type->no(0), type->vt, type->st, type->KDest, type->lm );
    REPORT::rpt.Output( text, 5 );
  }


  REPORT::rpt.Output( "\n\n  no    estx    esty\n", 5 );

  for( int i=0; i<nTypes; i++ )
  {
    int no;

    textLine = tbl->nextLine();

    sscanf( textLine, " %d", &no );
    TYPE* type = TYPE::Getno( no );

    sscanf( textLine, " %d %lf %lf",
                      &no, &(type->estx), &(type->esty) );

    sprintf( text, "  %2d   %6.3lf  %6.3lf", type->no(0), type->estx, type->esty );
    REPORT::rpt.Output( text, 5 );
  }
}


void TYPE::Import_30900( ASCIIFILE* tbl )
{
  enum
  {
    kROUGHNESS,  kTURBULENCE,  kSEDIMENT
  };

  DATKEY datkey[] =
  {
    kROUGHNESS,   "ROUGHNESS",     // 00
    kTURBULENCE,  "TURBULENCE",    // 01
    kSEDIMENT,    "SEDIMENT"       // 02
  };

  int nkey = 3;

  // -------------------------------------------------------------------------------------

  int   nRoughness  = 0;
  int   nTurbulence = 0;
  int   nSediment   = 0;

  char* textLine;


  // -------------------------------------------------------------------------------------
  // the file is read three times

  for( int rd=0; rd<3; rd++ )
  {
    while( !feof(tbl->getid()) )
    {
      if( !(textLine = tbl->nextLine()) )  break;

      char key[50];
      int  ikey = -99;

      // remove leading blanks
      while( *textLine == ' ' )  textLine++;

      if( textLine[0] == '$' )
      {
        sscanf( textLine, "$%s", key );

        for( int i=0; i<nkey; i++ )
        {
          if( strcmp(key,datkey[i].name) == 0 )
          {
            ikey = datkey[i].id;
            break;
          }
        }
      }

      switch( ikey )
      {
        case kROUGHNESS:
          if( rd == 0 )
          {
            int no;
            sscanf( textLine, "$ROUGHNESS %d", &no );
            if( no > 0 )  nRoughness++;
          }
          else if( rd == 1 )
          {
            int    no;
            int    idummy;
            double ddummy;

            sscanf( textLine, "$ROUGHNESS %d", &no );

            if( no > 0 )
            {
              TYPE* type = TYPE::Getid( nRoughness+1 );

              sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %lf %lf",
                                &no,
                                &type->rtype, &type->rcoef,
                                &idummy,      &ddummy,
                                &type->dp[0], &type->sp[0] );
              type->Setno( no );

              nRoughness++;
            }
            else
            {
              sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %lf %lf",
                                &no,
                                &deflt.rtype, &deflt.rcoef,
                                &idummy,      &ddummy,
                                &deflt.dp[0], &deflt.sp[0] );
            }
          }
          break;

        case kTURBULENCE:
          if( rd == 2 )
          {
            int no;
            sscanf( textLine, "$TURBULENCE %d", &no );

            if( no > 0 )
            {
              TYPE* type = TYPE::Getno( no );

              sscanf( textLine, "$TURBULENCE %d %lf %lf %lf %lf",
                                &no, &type->vt, &type->st, &type->KDest, &type->lm );

              nTurbulence++;
            }
            else
            {
              sscanf( textLine, "$TURBULENCE %d %lf %lf %lf %lf",
                                &no, &deflt.vt, &deflt.st, &deflt.KDest, &deflt.lm );
            }
          }
          break;

        case kSEDIMENT:
          if( rd == 2 )
          {
            int no;
            sscanf( textLine, "$SEDIMENT %d", &no );

            if( no > 0 )
            {
              TYPE* type = TYPE::Getno( no );

              sscanf( textLine, "$SEDIMENT %d %lf %lf", // %lf %lf %lf %lf %lf",
                                &no, &type->estx, &type->esty );
                              //&type->M, &type->tauc, &type->taus,
                              //&type->us, &type->rhob );

              nSediment++;
            }
            else
            {
              sscanf( textLine, "$SEDIMENT %d %lf %lf", // %lf %lf %lf %lf %lf",
                                &no, &deflt.estx, &deflt.esty );
                              //&deflt.M, &deflt.tauc, &deflt.taus,
                              //&deflt.us, &deflt.rhob );
            }
          }
          break;
      }
    }

    if( rd == 0 )
    {
      TYPE::Init( nRoughness );
      nRoughness = 0;
    }

    tbl->rewind();
  }


  // report roughness specifications -----------------------------------------------------

  char  text[400];

  sprintf( text, "\n   %d element type specifications:\n", TYPE::Get() );
  REPORT::rpt.Output( text, 5 );

  REPORT::rpt.Output( "\n   no   bottom        dp       sp\n", 5 );

  for( int i=0; i<TYPE::Get(); i++ )
  {
    TYPE* type = TYPE::Getid( i+1 );

    sprintf( text, "   %2d   %1d %8.2le",
                   type->no(0),
                   type->rtype,  type->rcoef );
    REPORT::rpt.Output( text, 5 );

    sprintf( text, "   %7.4lf  %7.4lf\n", type->dp[0], type->sp[0] );
    REPORT::rpt.Output( text, 5 );
  }


  REPORT::rpt.Output( "\n\n   no       vt       st    estx    esty     KDest    lm\n", 5 );

  for( int i=0; i<TYPE::Get(); i++ )
  {
    TYPE* type = TYPE::Getid( i+1 );

    sprintf( text, "   %2d   %9.6lf  %6.3lf  %6.3lf  %6.3lf  %6.3lf  %6.3lf\n",
                   type->no(0), type->vt, type->st, type->estx, type->esty, type->KDest, type->lm );
    REPORT::rpt.Output( text, 5 );
  }
}


void TYPE::Import_40000( ASCIIFILE* tbl )
{
  enum
  {
    kBOTTOM=1, kFORM, kVEGETATION, kTURBULENCE, kROUGHNESS
  };

  DATKEY datkey[] =
  {
    kBOTTOM,      "BOTTOM",        // 1
    kFORM,        "FORM",          // 2
    kVEGETATION,  "VEGETATION",    // 3
    kTURBULENCE,  "TURBULENCE",    // 4
    kROUGHNESS,   "ROUGHNESS"      // 5
  };

  int nkey = 5;

  // -------------------------------------------------------------------------------------

  int   nBottom     = 0;
  int   nTurbulence = 0;
  int   nRoughness  = 0;

  char* textLine;


  // -------------------------------------------------------------------------------------
  // the file is read three times
  // rd == 0: obtain the number of different roughness classes "nRoughness"
  //       1: read roughness parameter and set number of material zone
  //       2: read parameters for TURBULENCE

  for( int rd=0; rd<3; rd++ )
  {
    while( !feof(tbl->getid()) )
    {
      if( !(textLine = tbl->nextLine()) )  break;

      char key[50];
      int  ikey = -99;

      // remove leading blanks
      while( *textLine == ' ' )  textLine++;

      if( textLine[0] == '$' )
      {
        sscanf( textLine, "$%s", key );

        for( int i=0; i<nkey; i++ )
        {
          if( strcmp(key,datkey[i].name) == 0 )
          {
            ikey = datkey[i].id;
            break;
          }
        }
      }

      switch( ikey )
      {
        case kBOTTOM:
          if( rd == 0 )
          {
            int no;
            sscanf( textLine, "$BOTTOM %d", &no );
            if( no > 0 ) nBottom++;
          }
          else if( rd == 1 )
          {
            int    no;
            char   rlaw[10];
            sscanf( textLine, "$BOTTOM %d", &no );

            if( no > 0 )
            {
              TYPE* type = TYPE::Getid( nBottom+1 );

              sscanf( textLine, "$BOTTOM %d %s %lf %d %lf",
                      &no, rlaw, &type->rcoef, &type->kslaw, &type->ksfact );

              if(      strcmp(rlaw,"COWI") == 0 )  type->rtype = kCOWI;
              else if( strcmp(rlaw,"NIKU") == 0  ) type->rtype = kNIKU;
              else if( strcmp(rlaw,"CHEZ") == 0  ) type->rtype = kCHEZ;
              else if( strcmp(rlaw,"MANN") == 0  ) type->rtype = kMANN;

              type->Setno( no );

              nBottom++;
            }
          }
          break;

        case kFORM:
          if( rd == 2 )
          {
            int    no;
            sscanf( textLine, "$FORM %d", &no );

            if( no > 0 )
            {
              TYPE* type = TYPE::Getno( no );

              sscanf( textLine, "$FORM %d %lf %lf %d %lf",
                      &no, &type->d90, &type->d50, &type->form, &type->duneCoef );
            }
          }
          break;

        case kVEGETATION:
          if( rd == 2 )
          {
            int    no;
            sscanf( textLine, "$VEGETATION %d", &no );

            if( no > 0 )
            {
              TYPE* type = TYPE::Getno( no );

              sscanf( textLine, "$VEGETATION %d %lf %lf %lf %lf %lf",
                      &no, &type->dp[0], &type->sp[0], &type->hp,
                           &type->dp[1], &type->sp[1] );
            }
          }
          break;

        case kROUGHNESS:
          if( rd == 0 )
          {
            int no;
            sscanf( textLine, "$ROUGHNESS %d", &no );
            if( no > 0 )  nRoughness++;
          }
          else if( rd == 1 )
          {
            int    no;
            int    idummy;
            double ddummy;
            sscanf( textLine, "$ROUGHNESS %d", &no );

            if( no > 0 )
            {
              TYPE* type = TYPE::Getid( nRoughness+1 );

              if( release >= 40075 )
              {
                sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %d %lf %lf %d %lf %lf %lf %lf %lf %lf",
                                  &no,
                                  &type->rtype, &type->rcoef,
                                  &idummy,      &ddummy,
                                  &type->kslaw, &type->d90,   &type->d50,
                                  &type->form,  &type->duneCoef,
                                  &type->dp[0], &type->sp[0],
                                  &type->hp,    &type->dp[1], &type->sp[1] );
              }
              else if( release >= 40070 )
              {
                sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %d %lf %lf %lf %lf %lf %lf",
                                  &no,
                                  &type->rtype, &type->rcoef,
                                  &idummy,      &ddummy,
                                  &type->form,  &type->duneCoef,
                                  &type->dp[0], &type->sp[0],
                                  &type->hp,    &type->dp[1],    &type->sp[1] );
                type->d90 = type->rcoef;
              }
              else if( release >= 40060 )
              {
                sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %lf %lf %lf %lf %lf %d",
                                  &no,
                                  &type->rtype,    &type->rcoef,
                                  &idummy,         &ddummy,
                                  &type->dp[0], &type->sp[0],
                                  &type->hp, &type->dp[1], &type->sp[1], &type->form );
                type->d90 = type->rcoef;
              }
              else
              {
                sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %lf %lf",
                                  &no,
                                  &type->rtype,    &type->rcoef,
                                  &idummy,         &ddummy,
                                  &type->dp[0], &type->sp[0] );
                type->d90 = type->rcoef;
              }
              type->Setno( no );

              nRoughness++;
            }
            else
            {
              if( release >= 40075 )
              {
                sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %d %lf %lf %d %lf %lf %lf %lf %lf %lf",
                                  &no,
                                  &deflt.rtype, &deflt.rcoef,
                                  &idummy,      &ddummy,
                                  &deflt.kslaw, &deflt.d90, &deflt.d50,
                                  &deflt.form,  &deflt.duneCoef,
                                  &deflt.dp[0], &deflt.sp[0],
                                  &deflt.hp, &deflt.dp[1], &deflt.sp[1] );
              }
              else if( release >= 40070 )
              {
                sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %d %lf %lf %lf %lf %lf %lf",
                                  &no,
                                  &deflt.rtype,    &deflt.rcoef,
                                  &idummy,         &ddummy,
                                  &deflt.form, &deflt.duneCoef,
                                  &deflt.dp[0], &deflt.sp[0],
                                  &deflt.hp, &deflt.dp[1], &deflt.sp[1] );
                deflt.d90 = deflt.rcoef;
              }
              else if( release >= 40060 )
              {
                sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %lf %lf %lf %lf %lf",
                                  &no,
                                  &deflt.rtype,    &deflt.rcoef,
                                  &idummy,         &ddummy,
                                  &deflt.dp[0], &deflt.sp[0],
                                  &deflt.hp, &deflt.dp[1], &deflt.sp[1] );
                deflt.d90 = deflt.rcoef;
              }
              else
              {
                sscanf( textLine, "$ROUGHNESS %d %d %lf %d %lf %lf %lf",
                                  &no,
                                  &deflt.rtype,    &deflt.rcoef,
                                  &idummy,         &ddummy,
                                  &deflt.dp[0], &deflt.sp[0] );
                deflt.d90 = deflt.rcoef;
              }
            }
          }
          break;

        case kTURBULENCE:
          if( rd == 2 )
          {
            int no;
            sscanf( textLine, "$TURBULENCE %d", &no );

            if( no > 0 )
            {
              TYPE* type = TYPE::Getno( no );

              sscanf( textLine, "$TURBULENCE %d %lf %lf %lf %lf %lf %lf %lf",
                                &no, &type->vt, &type->st,
                                     &type->estx, &type->esty,
                                     &type->KDest, &type->lm, &type->betaSf );

              nTurbulence++;
            }
            else
            {
              sscanf( textLine, "$TURBULENCE %d %lf %lf %lf %lf %lf %lf %lf",
                                &no, &deflt.vt, &deflt.st,
                                     &deflt.estx, &deflt.esty,
                                     &deflt.KDest, &deflt.lm, &deflt.betaSf );
            }
          }
          break;
      }
    }

    if( rd == 0 )
    {
      if( nBottom > nRoughness ) nRoughness = nBottom;

      TYPE::Init( nRoughness );

      nRoughness = 0;
      nBottom    = 0;
    }

    tbl->rewind();
  }


  // report roughness specifications -----------------------------------------------------

  char  text[400];

  sprintf( text, "\n   %d element type specifications:\n\n", TYPE::Get() );
  REPORT::rpt.Output( text, 5 );

  REPORT::rpt.Output( " grain roughness\n", 5 );
  REPORT::rpt.Output( "   no          type  coef     kslaw ksfact\n", 5 );

  for( int i=0; i<TYPE::Get(); i++ )
  {
    TYPE* type = TYPE::Getid( i+1 );

    char rlaw[5];
    switch( type->rtype )
    {
      case kCOWI:  strcpy( rlaw, "COWI" ); break;
      case kNIKU:  strcpy( rlaw, "NIKU" ); break;
      case kMANN:  strcpy( rlaw, "MANN" ); break;
      case kCHEZ:  strcpy( rlaw, "CHEZ" ); break;
    }
    sprintf( text, "   %9d  %4s  %8.2le   %1d   %6.3lf\n",
                   type->no(0), rlaw,  type->rcoef, type->kslaw, type->ksfact );
    REPORT::rpt.Output( text, 5 );
  }

  REPORT::rpt.Output( "\n" );
  REPORT::rpt.Output( " form roughness\n", 5 );
  REPORT::rpt.Output( "   no          d90       d50      form  duneCf\n", 5 );
  for( int i=0; i<TYPE::Get(); i++ )
  {
    TYPE* type = TYPE::Getid( i+1 );

    sprintf( text, "   %9d  %8.6lf  %8.6lf   %1d   %6.4lf\n",
                   type->no(0), type->d90, type->d50, type->form, type->duneCoef );
    REPORT::rpt.Output( text, 5 );
  }

  REPORT::rpt.Output( "\n" );
  REPORT::rpt.Output( " vegetation parameter\n", 5 );
  REPORT::rpt.Output( "   no         dp1       sp1       hp2     dp2       sp2 \n", 5 );

  for( int i=0; i<TYPE::Get(); i++ )
  {
    TYPE* type = TYPE::Getid( i+1 );

    sprintf( text, "   %9d    %8.4lf %8.4lf", type->no(0), type->dp[0], type->sp[0] );
    REPORT::rpt.Output( text, 5 );

    sprintf( text, "  %6.4lf %8.4lf %8.4lf\n", type->hp,  type->dp[1], type->sp[1] );
    REPORT::rpt.Output( text, 5 );
  }

  REPORT::rpt.Output( "\n" );
  REPORT::rpt.Output( " turbulence parameter\n", 5 );
  REPORT::rpt.Output( "   no         vt         st      estx    esty    KDest   lm      betaSf\n", 5 );
  for( int i=0; i<TYPE::Get(); i++ )
  {
    TYPE* type = TYPE::Getid( i+1 );

    sprintf( text, "   %9d   %9.6lf  %6.3lf  %6.3lf  %6.3lf  %6.3lf  %6.3lf  %6.3lf\n",
             type->no(0), type->vt, type->st, type->estx, type->esty, type->KDest,
             type->lm, type->betaSf );
    REPORT::rpt.Output( text, 5 );
  }
}
