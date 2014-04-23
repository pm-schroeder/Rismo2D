// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software.
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// ======================================================================================
//
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
