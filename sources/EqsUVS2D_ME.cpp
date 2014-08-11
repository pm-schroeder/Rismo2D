// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class EQS_UVS2D_ME
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

#include "EqsUVS2D_ME.h"


EQS_UVS2D_ME::EQS_UVS2D_ME() : EQS_UVS2D()
{
  neq = 0;

  tri_id[0]  = 0;
  tri_id[1]  = 1;
  tri_id[2]  = 2;
  tri_id[3]  = 6;

  quad_id[0] = 0;
  quad_id[1] = 1;
  quad_id[2] = 2;
  quad_id[3] = 3;
  quad_id[4] = 8;
}


EQS_UVS2D_ME::~EQS_UVS2D_ME()
{
}
