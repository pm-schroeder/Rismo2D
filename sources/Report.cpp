// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class REPORT
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "Report.h"


REPORT REPORT::rpt;


//////////////////////////////////////////////////////////////////////////////////////////

REPORT::REPORT()
{
  pid   = 0;
  level = 0;
  out   = NULL;
  ntm   = 0;
  first = 1;
  last  = 1000000;
  step  = 1;
}

//////////////////////////////////////////////////////////////////////////////////////////

REPORT::~REPORT()
{
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Open( char* text, int level, int pid )
{
  this->pid   = pid;
  this->level = level;

  out = fopen( text, "w" );
  setvbuf( out, NULL, _IONBF, 0 );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Close()
{
  fprintf( out, "\n" );
  fclose( out );
  out = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Error( const char* format, ... )
{
  char text[200];

  va_list argPtr;
  va_start( argPtr, format );

  sprintf( text, "\n=== ERROR: " );  Text( text );
  vsprintf( text, format, argPtr );  Text( text );
  Text( "\n" );

  va_end( argPtr );

# ifdef _MPI_
  MPI_Finalize();
# endif

  exit( 0 );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Error( const int theError, const char* format, ... )
{
  char text[200];

  va_list argPtr;
  va_start( argPtr, format );

  sprintf( text, "\n=== ERROR %d: ", theError );  Text( text );
  vsprintf( text, format, argPtr );               Text( text );
  Text( "\n" );

  va_end( argPtr );

# ifdef _MPI_
  MPI_Finalize();
# endif

  exit( 0 );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Warning( const int theError, const char* format, ... )
{
  char text[200];

  va_list argPtr;
  va_start( argPtr, format );

  if( theError > 0 )
  {
    Text( "\n" );
    sprintf( text, "=== WARNING %3d: ", theError );  Text( text );
    vsprintf( text, format, argPtr );                Text( text );
    Text( "\n" );
  }
  else
  {
    sprintf( text, "                 " );  Text( text );
    vsprintf( text, format, argPtr );      Text( text );
    Text( "\n" );
  }

  va_end( argPtr );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Message( int level, const char* format, ... )
{
  char text[1000];

  va_list argPtr;
  va_start( argPtr, format );

  vsprintf( text, format, argPtr );

  va_end( argPtr );

  if( level <= this->level ) Text( text );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Text( const char* text )
{
  if( out )  fprintf( out, "%s", text );
  if( pid == 0 )  printf( "%s", text );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Screen( int level, const char* format, ... )
{
  if( pid == 0  &&  level <= this->level )
  {
    char text[1000];

    va_list argPtr;
    va_start( argPtr, format );

    vsprintf( text, format, argPtr );

    va_end( argPtr );

    printf( "%s", text );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Output( const char *text, int level, int pid )
{
  if( out && level <= this->level )  fprintf( out, "%s", text );
}


//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::ScreenLine1( int level )
{
  Screen( level, "--------------------------------------------------------------------------------\n" );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::OutputLine1( int level )
{
  Output( "--------------------------------------------------------------------------------\n", level );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Line1( int level )
{
  Message( level, "--------------------------------------------------------------------------------\n" );
}


//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::ScreenLine2( int level )
{
  Screen( level, "================================================================================\n" );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::OutputLine2( int level )
{
  Output( "================================================================================\n", level );
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Line2( int level )
{
  Message( level, "================================================================================\n" );
}

//////////////////////////////////////////////////////////////////////////////////////////

char* REPORT::License()
{
  static char license[] =
  "\
  All rights reserved.\n\n\
  This program is free software; you can redistribute it\n\
  under the terms of the Rismo Public License / Version 1.0.\n\n\
  This program is distributed in the hope that it will be useful,\n\
  but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n\
  Rismo Public License for more details.\n\n\
  You should have received a copy of the Rismo Public License\n\
  Version 1.0 along with this program; if not, write to:\n\
  Dr. Schroeder; Sperberweg 22/1; 75045 Walzbachtal; Germany\n\
  Email: rismo@hnware.com\n";

  return license;
}

//////////////////////////////////////////////////////////////////////////////////////////

void REPORT::Copyright( double release, int screen )
{
  char text[100];
  char chrel[20];
  sprintf( chrel, "%-9.3lf", release );

# ifdef _MPI_
  sprintf( text, "[RISMO_2D # Release %c.%c%c.%c%c-%c%c%c (MPI) # Copyright (C) 1992-2014, P.M. Schroeder]\n\n",
           chrel[0], chrel[1], chrel[2], chrel[3], chrel[4], chrel[6], chrel[7], chrel[8] );
# else
  sprintf( text, "[RISMO_2D # Release %c.%c%c.%c%c-%c%c%c # Copyright (C) 1992-2014, P.M. Schroeder]\n\n",
           chrel[0], chrel[1], chrel[2], chrel[3], chrel[4], chrel[6], chrel[7], chrel[8] );
# endif

  if( screen )
  {
    Screen( 0, text );

    ScreenLine2( 0 );
    Screen( 0, License() );
    ScreenLine2( 0 );

    Screen( 0, "\n" );
  }
  else
  {
    Output( text );

    OutputLine2( 0 );
    Output( License() );
    OutputLine2( 0 );

    Output( "\n", 0 );
  }
}
