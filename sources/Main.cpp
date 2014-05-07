////////////////////////////////////////////////////////////////////////////////////////////////////
//                                  R I S M O  2 D
//
// =================================================================================================
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
// =================================================================================================
//
//      Dr. Schroeder
//      Sperberweg 22/1
//      75045 Walzbachtal
//      Germany
//
//      Email: rismo@hnware.com
//
// =================================================================================================
// History
// -------
//                                    R I S M O
//
// The Finite Element program Rismo is designed to solve the time averaged
// Navier-Stokes-equations (Reynolds equations). For the first time it was
// developed at the
//
//                  Institute of Hydraulic Engineering
//                                 and
//                      Water Resources Development
//
//              University of Technology in Aachen, Germany
//
// for the purpose of natural river simulation. The work was financed
// by the German Research Foundation (DFG Ro 365/31-6).
//
// I would like to convey my thanks to all those people who have supported
// the development of this computer code. I'm very gratefule to Ian P. King
// who has been so kind to make the source code of his FE programs (RMAxx)
// available. This code has been a very valuable piece of illustration during
// the process of development.
//
//                                      Paul Michael Schroeder, November 1992
//
//
// -------------- 10.11.1992 -----------------------------------------------------------------------
// Version 0.0    is a two-dimensional depth averaged implementation,
//                which does use a simple eddy viscosity model to take
//                account for turbulence mechanisms in natural rivers.
//                A special algorithm was integrated to treat friction
//                losses from flood plain flow through non-submerged
//                vegetation (function lindner()).
//
// -------------- 14.11.1992 -----------------------------------------------------------------------
// Version 1.0    includes the k-epsilon turbulence model to determine
//                a variable eddy viscosity over the flow field. The
//                implementation of the k-epsilon model was previously
//                tested by Dr. C.J. Stein in the original RMA2 code.
//                Changes have been made to the treatment of Reynolds
//                stresses in depth averaged Navier-Stokes equations.
//                The k-epsilon model was stabilized by a pre-initia-
//                lization with a linearized version.
//
// -------------- 12.12.1992 -----------------------------------------------------------------------
// Version 1.1    includes the choice between linear and quadratic shape
//                functions (it was found that the linear integration
//                does not work in common cases)
//
//         1.1.1  enables the incremental transition from a constant to
//                a variable eddy viscosity.
//
// -------------- 23.03.1993 -----------------------------------------------------------------------
// Version 1.2    includes a dry and rewet algorithm. All elements with
//                one or more dry nodes (flow depth less than a specified
//                limit for drying) are rejected from the computation.
//
// -------------- 29.03.1993 -----------------------------------------------------------------------
// Version 1.3    includes the computation of instationary flow.
//
//         1.3.x  errors in dry/rewet algorithm fixed
//
// -------------- 01.07.1993 -----------------------------------------------------------------------
// Version 1.4    includes the Lagrange multiplier equation to determine
//                a divergence free flow field.
//
// -------------- 01.10.1993 -----------------------------------------------------------------------
// Version 1.5    has changed the computation of the friction coefficient
//                for bottom roughness and flow through non-submerged
//                vegetation from element to node level.
//
//
// =================================================================================================
//
//                                    R I S M O  2 D
//
// Starting with Version 2.0 the computer code has been developed to a
// threedimensional Finite-Element-Program. Because of the increasing
// complexity of the code I decided to split the twodimensional and
// threedimensional parts. And I named the twodimensional part Rismo2D
// which got the first version number 0.3.
//
//                                      Paul Michael Schroeder, Januar 1996
//
// ---------------- 10.11.1996 ---------------------------------------------------------------------
// Version 1.2      has implemented some new features:
//                  - a new dry/rewet algorithm, which keeps partially
//                    dry elements within the system matrix
//                  - a new very stable relaxation for steady and unsteady
//                    flow analysis
//                  - different data files for geometry and numeric results
//
// ---------------- 21.10.1997 ---------------------------------------------------------------------
// Version 2.0      fixed errors (due to the upgrade to Version 1.2)
//                  - a geometry output file was added to the list of file
//                    names, which holds the changes, made to the ordering
//                    of elements (reorder cycle)
//
// ---------------- 20.11.1997 ---------------------------------------------------------------------
// Version 2.1      changes in the generation of boundary conditions
//
// ----------------  8.12.1997 ---------------------------------------------------------------------
// Version 2.2      statistic files reestablished
//
// ---------------- 18.02.1998 ---------------------------------------------------------------------
// Version 2.2.1    minor changes; output of changed bottom elevation due
//                  to erosion/sedimentation; use of node->zor instead of
//                  node->z in all output files
//
// ---------------- 19.02.1998 ---------------------------------------------------------------------
// Version 2.2.2    minor changes in initialization of boundary conditions
//                  for nodes (initBC)
//
// ---------------- 27.02.1998 ---------------------------------------------------------------------
// Version 2.3      advection-diffusion computation adapted to the dry-
//                  rewet-method 2 (marsh-nodes)
//
// ---------------- 20.03.1998 ---------------------------------------------------------------------
// Version 2.4      dynamic memory allocation changed from C-functions
//                  malloc() and free() to C++ allocators new and delete
//
// ---------------- 18.03.1999 ---------------------------------------------------------------------
// Version 2.5.1    possible error in function InitS() (initialization of
//                  water surface from sections) removed
//
// ---------------- 20.05.1999 ---------------------------------------------------------------------
// Version 2.5.2    functions assemble..() and setIndexMat() adapted to
//                  solve also for element equations in doUVSL()
//
// ---------------- 30.08.1999 ---------------------------------------------------------------------
// Version 2.5.5    scaling factors read from input data file
//
// ---------------- 12.10.1999 ---------------------------------------------------------------------
// Version 2.6.1    last version before C++ classes: PROJECT, MODEL, GRID
//
// ---------------- 01.12.1999 ---------------------------------------------------------------------
// Version 2.7.1    version adapted for supercomputer Origin 2000
//
// =================================================================================================
//
// ---------------- 27.10.2000 ---------------------------------------------------------------------
// Version 3.00.00  further C++ classes
//
// ---------------- 19.07.2002 ---------------------------------------------------------------------
// Version 3.02.01  number for element material (in file *.dat) can be chosen freely
//                  internally the material number is no longer used as array index
//                  changes in: TYPE and all methods using TYPE
//
// ---------------- 27.11.2002 ---------------------------------------------------------------------
// Version 3.02.02  outflow boundary condition according to a specified relation
//                  of flowdepth and discharge may be specified
//
// ---------------- 03.12.2002 ---------------------------------------------------------------------
// Version 3.2.3    changes made to inflow condition in method setInflow()
//
// ---------------- 16.12.2002 ---------------------------------------------------------------------
// Version 3.02.04  instability in new outflow condition (see 3.2.2) relaxed
//
// ---------------- 12.09.2003 ---------------------------------------------------------------------
// Version 3.03.04  cycle kSurfaceToVol = 35; to compute surface in elements (ELEM::S)
//
// ---------------- 15.12.2003 ---------------------------------------------------------------------
// Version 3.03.05  Errors in initialisation of ELDERs model: GRID::eddyViso()
//                  initialisation of nd[i]->Exx, nd[i]->Exy und nd[i]->Eyy
//
// ---------------- 30.06.2004 ---------------------------------------------------------------------
// Version 3.04.05  some smaller changes
//
// ---------------- 03.07.2004 ---------------------------------------------------------------------
// Version 3.04.06  first release facing MPI-Version
//                  method PROJECT::Compute() developed from main()
//
// ---------------- 12.08.2004 ---------------------------------------------------------------------
// Version 3.06.06  MPI-version finished
//
// ---------------- 12.09.2004 ---------------------------------------------------------------------
// Version 3.07.00  implementation of 2 further solvers (PARMS):
//                  P_bcgstad() and P_fgmresd()
//
// ---------------- 26.09.2004 ---------------------------------------------------------------------
// Version 3.07.07  many code changes in preparation for Open Source Distribution
//
// ---------------- 24.10.2004 ---------------------------------------------------------------------
// Version 3.08.00  fixed bugs in MPI communication
//
// ---------------- 16.03.2005 ---------------------------------------------------------------------
// Version 3.09.00  input data files restyled
//
// ---------------- 21.10.2005 ---------------------------------------------------------------------
// Version 3.09.07  error in method TYPE::Friction() fixed (Colebrook/White)
//
// ---------------- 27.11.2005 ---------------------------------------------------------------------
// Version 3.09.08  Windows-Version with Icons und Setup
//
// ---------------- 17.01.2006 ---------------------------------------------------------------------
// Version 3.09.10  1. Changed sign for specification of boundary condition kInlet. Now
//                     the user has to specify positive values for q = H*U on inlet.
//                     Method: BCONLINE::GenBcon()
//                  2. Class TYPE now writes out a list of up to 100 material zones
//                     which are not specified in the material table.
//
// ---------------- 18.04.2006 ---------------------------------------------------------------------
// Version 3.10.02  some changes in methods GRID::DryRewet() and GRID::RewetDry()
//
// ---------------- 01.07.2005 ---------------------------------------------------------------------
// Version 4.00.00  computation of bed load (based on diploma-thesis of T. Krahm)
//
// ---------------- 21.10.2005 ---------------------------------------------------------------------
// Version 4.00.07  error in method TYPE::Friction() fixed (Colebrook/White)
//                  older version 3.09.06 updated to version 3.09.07
//
// ---------------- 09.10.2005 ---------------------------------------------------------------------
// Version 4.00.21  parallel version of bed load transport
//
// ---------------- 27.10.2005 ---------------------------------------------------------------------
// Version 4.00.23  implementation of parallel frontal solver (class FRONTM) corrected
//
// ---------------- 06.11.2005 ---------------------------------------------------------------------
// Version 4.00.24  class EQS_DZ: no bottom evolution of dry nodes (flag NODE::kMarsh)
//
// -------------------------------------------------------------------------------------------------
// Version 4.00.**  further developments and testing of bed load transport
//
// ---------------- 05.12.2007 ---------------------------------------------------------------------
// Version 4.00.44  change of iteration scheme for instationary flow computation
//                  computation will be continued if time step is not gained
//                  method: EQS_UVS2D::Execute()
//
// ---------------- 05.12.2007 ---------------------------------------------------------------------
// Version 4.00.45  new feature: read start time from time step file ($TM_SETTIME)
//
// ---------------- 07.01.2008 ---------------------------------------------------------------------
// Version 4.00.46  call of methods BCONSET::SetInflow() and BCONSET::SetOutflow() changed
//                  to use time at the end of timestep to interpolate boundary condtions
//
// ---------------- 10.01.2008 ---------------------------------------------------------------------
// Version 4.00.48  internal computation of the Defaults Manning's value "nDef" to gain a
//                  smooth transition of friction coefficients at nodes with h = ks
//
// ---------------- 14.01.2008 ---------------------------------------------------------------------
// Version 4.00.49  some minor changes to the convergence criterion in EQS_UVS2D:Execute()
//
// ---------------- 11.02.2008 ---------------------------------------------------------------------
// Version 4.00.52  statistics for MPI-Version
//
// ---------------- 06.05.2008 ---------------------------------------------------------------------
// Version 4.00.60  friction due to submerged vegetation / van Velzen's approach
//
// ---------------- 11.05.2008 ---------------------------------------------------------------------
// Version 4.00.70  friction due to form roughness according to van Rijn
//
// ---------------- 19.07.2008 ---------------------------------------------------------------------
// Version 4.00.75  revised format of roughness table; introducing TYPE::d50 and
//                  TYPE::d90 which replaces grain roughness ks = TYPE::rcoef[0].
//                  see additional switches in TYPE::Import_40000()
//
// ---------------- 07.10.2009 ---------------------------------------------------------------------
// Version 4.00.78  friction computation in case of small flow depth h < 2*ks corrected
//
// ---------------- 04.01.2010 ---------------------------------------------------------------------
// Version 4.00.81  first implementation of diffuse Sources or Sinks at nodes at
//                  BCONSET::InitBcon(), EQS_UVS2D::Region(), EQS_UVS2D::RegionAI()
//
// ---------------- 18.01.2010 ---------------------------------------------------------------------
// Version 4.00.83  a) change to the outlet boundary condition:
//                     nodes which belong to closed (kAutoSlip-condition) as well as to
//                     outlet boundary elements (kOutlet or kOpenBnd-condition) will be
//                     assigned a flow direction normal to the outlet boundary
//                     MODEL::SetRotation()
//                  b) implementation of periodic boundary conditions
//                     key words: kTM_PERIODIC_LINE, kTM_PERIODIC_NODE (see Timeint.h)
//                     methods:   TIMEINT::TIMEINT(), TIMEINT::Input_30900()
//                                PROJECT::Compute()
//                  c) error in method STATIST::Read() corrected
//
// ---------------- 20.02.2010 ---------------------------------------------------------------------
// Version 4.00.84  some modifications to the k-eps-modell
//                  
// ---------------- 25.03.2010 ---------------------------------------------------------------------
// Version 4.00.85  adaptions to the key words for input parameters of minmax-values
//
// ---------------- 29.03.2010 ---------------------------------------------------------------------
// Version 4.01.00  fundamental changes to the implementation of boundary conditions
//
// ---------------- 27.10.2010 ---------------------------------------------------------------------
// Version 4.01.07  introducing min-values for directional eddy viscosity in Elder model
//
// ---------------- 27.05.2011 ---------------------------------------------------------------------
// Version 4.01.15  input error in TYPE::Import_40000() corrected (line 607):
//                  Read turbulence parameters: sscanf( ...&type->estx, &type->esty... )
//                  instead of sscanf( ...&deflt.estx, &deflt.esty...)
//
// ---------------- 22.06.2011 ---------------------------------------------------------------------
// Version 4.02.00  binary output of time series implemented; keyword: $TIMESERIES
//
// ---------------- 27.01.2012 ---------------------------------------------------------------------
// Version 4.03.05  some changes to the input data format of the startup files
//
// ---------------- 31.12.2012 ---------------------------------------------------------------------
// Version 4.04.04  Different time stepping scheme implemented in class EQS_UVS2D_TM
//
// ---------------- 03.01.2013 ---------------------------------------------------------------------
// Version 4.05.00  Advanced usage of reporting levels.
//                  The amount of output to the report file and screen is depending on the
//                  chosen report level. A higher level will increase the amount. The output
//                  may be reduced to certain time steps defined with the key kREPORTTIME.
//
//                  level = 1: output for each time step and cycle at the end of iterations:
//                             time-step number, cycle, time and date, CPU, convergence data,
//                             discharge through control lines
//
//                  level = 2: output for each time step, cycle and each iteration:
//                             time-step number, cycle, time and date, CPU, convergence data,
//                             discharge through control lines
//
//                  level = 3: regular output for each time step, cycle and each iteration
//
//                  level = 4: more output for each time step, cycle and each iteration
//
//                  level = 5: full output of any available information
// 
// ---------------- 19.02.2013 ---------------------------------------------------------------------
// Version 4.05.09  Some changes ...
//                  1. Statistical parameters:
//                     - analysis restricted to wet nodes
//                     - new parameters: mean water elevation (meanS)
//                                       flooding rate        (fldRate)
//                  2. Dry-Rewet algorithm:
//                     The "steady-fill"-parameter within the dry-rewet-method 2 and 3 has been
//                     replaced by the "rewet-limit"-parameter. A rewetting node is initialized
//                     with the rewet-limit (rewetLimit >= dryLimit).
//                  3. Lower limit from eddy viscosity has been changed to a local parameter.
//                     minVt is no longer accepted; use the constant eddy viscosity parameter in
//                     the material table to set a lower limit for the eddy viscosity instead
//
// ---------------- 12.10.2013 ---------------------------------------------------------------------
// Version 4.05.17  New boundary conditions. Keywords: kTARGET_S, kTARGET_ST
// -------------------------------------------------------------------------------------------------

#include "Defs.h"
#include "Report.h"
#include "Scale.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Grid.h"
#include "Model.h"
#include "Section.h"
#include "Statist.h"

#include "Project.h"

#include "Time.h"


int main( int  argc, char *argv[] )
{
  // initialize MPI ----------------------------------------------------------------------

  int pid = 0;
  int npr = 1;

# ifdef _MPI_

  MPI_Init( &argc, &argv );

  MPI_Comm_size( MPI_COMM_WORLD, &npr );
  MPI_Comm_rank( MPI_COMM_WORLD, &pid );

# endif


  // set the standard output to unbuffered -----------------------------------------------

  setvbuf( stdout, NULL, _IONBF, 0 );


  // create the project ------------------------------------------------------------------

  PROJECT* project = new PROJECT;

  if( !project )  REPORT::rpt.Error( kMemoryFault, "can not allocate memory (main - 1)" );

  project->subdom.pid = pid;
  project->subdom.npr = npr;

  if( pid == 0 )  REPORT::rpt.Copyright( project->release, true );


  // get input file name -----------------------------------------------------------------

  if( argc > 1 )
  {
    strcpy( project->name.inputFile, argv[1] );
  }
  else
  {
    if( pid == 0 )
    {
      printf( "\n Give name of INPUT file: ");
      scanf( "%s", project->name.inputFile );
      printf( "\n" );
    }

#   ifdef _MPI_
    MPI_Bcast( &project->name.inputFile, project->name.Length(),
               MPI_CHAR, 0, MPI_COMM_WORLD );
#   endif
  }

  project->Compute();


  // finalize MPI ------------------------------------------------------------------------

# ifdef _MPI_
  MPI_Finalize();
# endif

  return 0;
}
