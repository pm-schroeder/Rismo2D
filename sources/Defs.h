// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// D E F S
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Defs.h  : file containing some global used defines.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// A simple header file containing some basic defines.
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
//  01.01.1992     sc     first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef DEFS_INCL
#define DEFS_INCL

//#define _MPI_DBG
//#define _MPI_

#ifdef _MPI_
#include <mpi.h>
#endif

#define LINUX
//#define MS_WIN

//#define kRangeCheck
//#define kIteratCount

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>

#define kRelease     40521.002          // release

#define PI           3.141592653589793

#define false        0
#define true         1

#define kZero        1.0e-40
#define kEpsilon     1.0e-20

#define kMaxDF       6                  // maximum degree of freedom at nodes
#define kSimDF       5                  // simultaneous degree of freedom (maximum)

#define kMaxCycles  36                  // maximum number of cycles

#define REALPR float

// makros for bit manipulation -----------------------------------------------------------
#define SF(F,b)   ((F) |=  ((unsigned int) b))         // set flag
#define CF(F,b)   ((F) &= ~((unsigned int) b))         // clear flag
#define isFS(F,b) ((F)  &  ((unsigned int) b))         // is flag set?

#endif
