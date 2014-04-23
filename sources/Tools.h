// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// T O O L S
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Tools.h   : definition file of the class.
// Tools.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// Definition of some global functions.
//
// -------------------------------------------------------------------------------------------------
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
// -------------------------------------------------------------------------------------------------
//
// HISTORY
//
//    date              changes
// ------------  ----  -----------------------------------------------------------------------------
//  01.01.1992    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////// ======================================================================================

#include "Defs.h"

void    printTime( int );
void    initTheClock();
char*   convertTime( clock_t );

void    doubleToTime( double, char [18] );
double  timeToDouble( char [18] );
