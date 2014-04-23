// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class NODE
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

#include "Node.h"


NODE::NODE()
{
  no   = 0;
  name = 0;

  x    = 0.0;
  y    = 0.0;
  z    = 0.0;
  zor  = 0.0;

  dz   = 0.0;
  zero = 0.0;
  qbo  = 0.0;

  cf   = 0.0;
  cfw  = 0.0;

  vt   = 0.0;
  exx  = 0.0;
  exy  = 0.0;
  eyy  = 0.0;

  uu   = 0.0;
  uv   = 0.0;
  vv   = 0.0;

  Dxx  = 0.0;
  Dxy  = 0.0;
  Dyy  = 0.0;

  Vsec = 0.0;

  flag = 0;
  noel = 0;
  mark = 0;

  countDown = 0;

  sub  = NULL;
}


NODE::~NODE()
{
  if( sub )  delete sub;
}


NODE NODE::operator =( const NODE& n )
{
  x   = n.x;
  y   = n.y;
  z   = n.z;
  zor = n.zor;

  no        = n.no;
  name      = n.name;

  flag      = n.flag;
  sub       = n.sub;
  countDown = n.countDown;
  dz        = n.dz;
  cf        = n.cf;
  cfw       = n.cfw;
  vt        = n.vt;
  exx       = n.exx;
  exy       = n.exy;
  eyy       = n.eyy;
  uu        = n.uu;
  uv        = n.uv;
  vv        = n.vv;
  v         = n.v;
  vo        = n.vo;
  bc        = n.bc;

  noel      = 0;
  el        = NULL;

  return* this;
}
