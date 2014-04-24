// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// R E P O R T
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Report.h   : definition file of the class.
// Report.cpp : implementation file of the class.
//
// Time.cpp   : methods REPORT::InitTheClock()
//                      REPORT::CpuClock()
//                      REPORT::ConvertTime()
//                      REPORT::PrintTime()
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements the output of reports.
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
//  01.01.1994    sc    first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

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
