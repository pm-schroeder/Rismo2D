// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// M O D E L
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Model.h         : definition file of the class.
// Model.cpp       : implementation file of the class.
//
// Bound.cpp       : methods MODEL::Boundary()
//                           MODEL::Getbound()
// Contin.cpp      : method  MODEL::Continuity()
// ContinBSL.cpp   : method  MODEL::ContinuityBSL()
// Curv2D.cpp      : method  MODEL::Curv2D()
// DoDryRew.cpp    : methods MODEL::DoDryRewet()
//                           MODEL::MPI_Comm_Dry()
// DoFriction.cpp  : method  MODEL::DoFriction()
// LastNode.cpp    : method  MODEL::LastNode()
// Locate.cpp      : method  MODEL::SetLocation()
// Phi2D.cpp       : method  MODEL::Phi2D()
// ReorderElem.cpp : method  MODEL::ReorderElem()
// Rot2D.cpp       : method  MODEL::Rot2D()
// SetBdKD.cpp     : methods MODEL::SetBoundKD()
//                           MODEL::SetNodeKD()
//                           MODEL::SmoothKDBound()
// Transform.cpp   : methods MODEL::SetNormal()
//                           MODEL::SetRotation()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the Finite Element Model.
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
//  01.01.1998    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MODEL_INCL
#define MODEL_INCL

#include "Defs.h"
#include "Grid.h"

class PROJECT;
class SUBDOM;
class SED;
class TIME;
class TYPE;

struct KDCONST;
struct NAME;
struct TMSER;


class MODEL
{
  protected:
    ELEM**  boundList;
    int     init;             // counter increased each time when Initialize() was called

    double* phi;
    double* rot;
    double* curv;
    double* Lx;
    double* Ly;
    double* rc;
    double* d90;
    double* d50;
    double* Hr;
    double* Hd;
    double* kd;
    double* hp;
    double* dp;
    double* sp;

  public:
    GRID*   region;           // grid of region elements
    GRID*   control;          // grid of control elements
    GRID*   bound;            // grid of boundary elements

    int     np;
    NODE**  node;             // list of nodes
    int     ne;
    ELEM**  elem;             // list of elements including boundary elements

    ELEM*   list;

    SUBDOM* subdom;

  public:
    // Model.cpp ---------------------------------------------------------------------------
    MODEL();
    ~MODEL();

    void    Initialize();
    int     Getinit();
    void    Incinit();

    int     Input( PROJECT* project, NAME*, TIME* prevTime, SED* sed );
    void    Output( PROJECT* project, int timeStep );

    void    OutputSeriesHeader( PROJECT *project, TMSER* tmser );
    void    OutputSeries( PROJECT *project, int timeStep, TMSER* tmser, bool tmser_first_call );

    void    AttachOutput( PROJECT *project, int val );
    void    DetachOutput( PROJECT *project );

    void    UCDOutput( FILE* id, PROJECT* project, GRID* grid,
                       double* phi,  double* rot,  double* curv,
                       double* Lx,   double* Ly,
                       double* rc,   double* d90,   double* d50,
                       double* Hr,   double* Hd,    double* kd,
                       double* hp,   double* dp,    double* sp,
                       double* qbe,  double* Ls,    double* sx,   double* sy,
                       double* dzds, double* dzmx,  double* dhds );

//  The following methods have been moved to the class SUBDOM. SC, 16.04.2010
//  void    Mpi_assemble( double* vec, PROJECT* project );
//  void    Mpi_average( double* vec, PROJECT* project );

    // Bound.cpp ---------------------------------------------------------------------------
    void    Boundary();
    ELEM*   Getbound( int );

    // Contin.cpp --------------------------------------------------------------------------
    void    Continuity();

    // ContinBSL.cpp ------------------------------------------------------------------------
    void    ContinuityBSL();

    // DoDryRew.cpp ------------------------------------------------------------------------
    void    DoDryRewet( PROJECT* project, int* dried =NULL, int* wetted =NULL );
    void    MPI_Comm_Dry( PROJECT* project, int initialize );

    // DoFriction.cpp ----------------------------------------------------------------------
    void    DoFriction( PROJECT* );

    // LastNode.cpp ------------------------------------------------------------------------
    void    LastNode();

    // Locate.cpp --------------------------------------------------------------------------
    void    SetLocation();

    // ReorderElem.cpp ---------------------------------------------------------------------
    ELEM*   ReorderElem( int ns, SECTION* section );

    // Phi2D.cpp --------------------------------------------------------------------------
    double* Phi2D();

    // Curv2D.cpp --------------------------------------------------------------------------
    double* Curv2D();

    // Rot2D.cpp ---------------------------------------------------------------------------
    double* Rot2D();

    // SetBdKD.cpp -------------------------------------------------------------------------
    void    SetBoundKD( PROJECT* );
    void    SetNodeKD( PROJECT*, NODE*, int, int, double, double );
    void    SmoothKDBound( int );

    // Transform.cpp -----------------------------------------------------------------------
    void    SetNormal();
    void    SetRotation();
};

#endif
