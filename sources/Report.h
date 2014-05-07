// ======================================================================================
//                                    R E P O R T
// ======================================================================================
// This class implements the output of reports.
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
//    date               description
// ----------   ------   ----------------------------------------------------------------
// 01.01.1994     sc     first implementation / first concept
//
// ======================================================================================

// read this header file only once ...

#ifndef REPORT_INCL
#define REPORT_INCL

#include "Defs.h"

#define kUnexpectedFault  0
#define kMemoryFault      1
#define kVerifyFault      2
#define kParameterFault   3
#define kGridFault        4
#define kValueFault       5
#define kMethodFault      6
#define kRangeFault       7
#define kFileFault        8
#define kOpenFileFault    9
#define kReadFileFault   10
#define kWriteFileFault  11
#define kInputFault      12
#define kSolverFault     13
#define kUserFault       14
#define kMPI_Error       15

#define rptOPEN           1
#define rptCLOSE          2
#define rptSCREEN         4
#define rptFILE           8


class REPORT
{
  public:
    static REPORT  rpt;

    static clock_t cpuTime;
    static clock_t firstClock;
    static clock_t lastClock;
    static clock_t newClock;

    int    pid;
    int    level;
    FILE*  out;

    int    ntm;
    int*   tmlist;                        // array of time steps for output
    int    first;
    int    last;
    int    step;

  public:
    REPORT();
    virtual ~REPORT();

    void    Open( char* filename, int level, int pid );
    void    Close();

    void    Setlevel( int lv )  { level = lv; };

    void    Error( const char*, ... );
    void    Error( const int, const char*, ... );

    void    Warning( const int, const char*, ... );

    void    Message( int level, const char* text, ... );

    void    Output( const char* text, int level =0, int pid =0 );
    void    Screen( int level, const char* text, ... );

    void    Line1( int level );
    void    OutputLine1( int level=0 );
    void    ScreenLine1( int level );

    void    Line2( int level );
    void    OutputLine2( int level=0 );
    void    ScreenLine2( int level );

    char*   License();
    void    Copyright( int release, int screen );


  private:
    void    Text( const char* text );

  //----- Time.cpp -----------------------------------------------------------------------
  public:
    void    PrintTime( int );
    void    InitTheClock();
    char*   ConvertTime( clock_t );

  private:
    clock_t CpuClock();

};
#endif
