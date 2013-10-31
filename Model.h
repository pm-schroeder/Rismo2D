// ======================================================================================
//                                      M O D E L
// ======================================================================================
// This class implements the Finite Element Model.
// ======================================================================================
//
// Copyright (C) 1992-2012  by  P.M. SCHROEDER
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
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.1998     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once ...

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
    void    OutputSeries( PROJECT *project, int timeStep, TMSER* tmser );

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
