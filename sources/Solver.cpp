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
