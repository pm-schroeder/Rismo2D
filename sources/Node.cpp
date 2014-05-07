// ======================================================================================
//
// Copyright (C) 1992-2008  by  P.M. SCHROEDER
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
