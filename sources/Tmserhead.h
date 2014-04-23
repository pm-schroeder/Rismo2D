// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// T M S E R H E A D
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Tmserhead.h   : definition file of the class.
// Tmserhead.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the interface to the Rismo time series file (*.rts).
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
//  29.03.2010    sc    Rismo-Version 4.01.00, new keywords: kTM_NODE, kTM_LINE
//                      class RELOC to ensure backward compatibility
//  13.10.2012    sc    Rismo-Version 4.03.00, new keyword: kTM_STATIONARY
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TMSERHEAD_H
#define	TMSERHEAD_H

enum { kTMSERHEADSize = 17*sizeof(int) };

union TMSERHEAD
{
  char buffer[kTMSERHEADSize];
  struct
  {
    int    np;
    int    ne;
    int    first;
    int    last;
    int    step;
    int    startTime[5];
    int    deltaTime[5];
    int    vcomp;                  // number of data components in time series file
    int    vdata;
  };
};

#endif	// TMSERHEAD_H

