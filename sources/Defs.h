// ======================================================================================
//
// Copyright (C) 1992-2013  by  P.M. SCHROEDER
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

#define PI           3.141592653589793

#define false        0
#define true         1

#define kZero        1.0e-40
#define kEpsilon     1.0e-20

#define kMaxDF       6         // maximum degree of freedom at nodes
#define kSimDF       5         // simultaneous degree of freedom (maximum)

#define kMaxCycles  36         // maximum number of cycles

#define REALPR float

// makros for bit manipulation -----------------------------------------------------------
#define SF(F,b)   ((F) |=  ((unsigned int) b))         // set flag
#define CF(F,b)   ((F) &= ~((unsigned int) b))         // clear flag
#define isFS(F,b) ((F)  &  ((unsigned int) b))         // is flag set?

#endif
