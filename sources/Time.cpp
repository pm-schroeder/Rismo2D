// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// class REPORT
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


clock_t REPORT::cpuTime;
clock_t REPORT::firstClock;
clock_t REPORT::lastClock;
clock_t REPORT::newClock;


// ---------------------------------------------------------------------------------------

void REPORT::InitTheClock()
{
  cpuTime = (long) (clock()*100.0/CLOCKS_PER_SEC);  // in 10 milliseconds

  firstClock = lastClock = newClock = 0;
}


// ---------------------------------------------------------------------------------------

clock_t REPORT::CpuClock()
{
  clock_t cpu, cpuClock, diffCpu;
  double  resolution;

  cpuClock = clock();

  if ( cpuClock == ((clock_t) -1) )
    cpuClock = clock();                             // try again

  resolution = CLOCKS_PER_SEC / 100.0;

  cpu = (long) (cpuClock / resolution);             // in 10 milliseconds

  if ( cpu < cpuTime )
  {
    diffCpu = (long) (4294967296.0 / resolution     // SUN systems: wrap around
                      + cpu - cpuTime);
    cpuTime = cpu;
  }

  else
  {
    diffCpu = cpu - cpuTime;
    cpuTime = cpu;
  }

  return diffCpu;
}


// ---------------------------------------------------------------------------------------

char* REPORT::ConvertTime( clock_t msec )
{
  static char theTimer[20];
  long s[4];

  // determine number of milliseconds msec
                        s[0] = msec % 100;             // 10 milliseconds
  msec /= 100;          s[1] = msec % 60;              // seconds
  msec /=  60;          s[2] = msec % 60;              // minutes
  msec /=  60;          s[3] = msec;                   // hours

  sprintf( theTimer, "%4ld:%02ld:%02ld:%02ld", s[3], s[2], s[1], s[0] );

  return theTimer;
}


// ---------------------------------------------------------------------------------------

void REPORT::PrintTime( int level )
{
  clock_t theClock = CpuClock();

  newClock += theClock;


  if( level <= this->level )
  {
    clock_t totalCpu = newClock - firstClock;
    clock_t diffCpu  = newClock - lastClock;

    lastClock = newClock;

    char sCPU[20];    strcpy( sCPU,    ConvertTime(diffCpu) );
    char stotCPU[20]; strcpy( stotCPU, ConvertTime(totalCpu) );

    stotCPU[10] = '\0';

    char   date[26];
    time_t now  = time( NULL );
    strcpy ( date, ctime( &now ) );
    date[24] = '\0';

    REPORT::rpt.Message( level, "\n %23s                                    CPU: %13s\n %s %10s\n",
                                date, sCPU,
                                "                                                      total CPU:",
                                stotCPU );
  }
}
