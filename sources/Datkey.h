// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// D A T K E Y
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Datkey.h   : definition file of the class.
// Datkey.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// A simple data wrapper.
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
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.200x     sc     first implementation
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef DATKEY_INCL
#define DATKEY_INCL

struct DATKEY
{
  int         id;
  const char* name;
};

#endif
