// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class SOLVER
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
#include "Elem.h"
#include "Eqs.h"

#include "Solver.h"


int      SOLVER::m_neqs   = 0;
SOLVER** SOLVER::m_solver = NULL;


SOLVER::SOLVER()
{
  solverType  = 0;
  preconType  = 0;

  path[0]     = '\0';

  no          = 0;
  solverType  = 1;
  eqs         = NULL;

  mfw         = 500;
  size        = 0;

  preconType  = 0;
  mceq        = 100;
  proceed     = 0;
  maxIter     = 1000;
  maxDiff     = 0.001;

  mkyrl       = 10;

  accuracy    = 1.0;
  iterCountCG = 0;
};


SOLVER::~SOLVER()
{
};
