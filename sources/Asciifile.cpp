// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class ASCIIFILE
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
#include "Report.h"
#include "Asciifile.h"


ASCIIFILE::ASCIIFILE( const char* fileName, const char* mode, int maxLineLength )
{
  line = new char[maxLineLength];
  name = new char[strlen(fileName) + 1];

  if( !line  ||  !name )
    REPORT::rpt.Error( kMemoryFault, "cannot allocate memory (ASCIIFILE - 1)" );

  strcpy( name, fileName );


  id = fopen( name, mode );

  length = maxLineLength;
  lineNo = 0;
}


ASCIIFILE::~ASCIIFILE()
{
  delete[] line;
  delete[] name;

  if( id )  fclose( id );
}


FILE* ASCIIFILE::getid()      { return id; }
int   ASCIIFILE::getlineNo()  { return lineNo; }
char* ASCIIFILE::getname()    { return name; }
void  ASCIIFILE::rewind()     { fseek(id, 0L, SEEK_SET); lineNo=0; }


char* ASCIIFILE::next()
{
  if( feof( id ) )  return (char*) 0;

  if( !fgets(line, length, id) ) return (char*) 0;

  lineNo++;

  return line;
}


char* ASCIIFILE::nextLine()
{
  int  comment;

  comment = true;

  while( comment )
  {
    if( feof( id ) )  return (char*) 0;

    if( !fgets(line, length, id) )  return (char*) 0;

    if( line[0] != '*'  &&  line[0] != '#')  comment = false;

    lineNo++;
  }

  return line;
}


void ASCIIFILE::print( const char* format, ... )
{
  va_list argPtr;
  va_start( argPtr, format );

  vfprintf( id, format, argPtr );

  va_end( argPtr );
}


void ASCIIFILE::echoPrint( const char* format, ... )
{
  va_list argPtr;
  va_start( argPtr, format );

  vfprintf( id, format, argPtr );
  vprintf( format, argPtr );

  va_end( argPtr );
}
