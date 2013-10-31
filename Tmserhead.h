/* 
 * File:   tmserhead.h
 * Author: schroedr
 *
 * Created on 6. Juli 2011, 15:23
 */

#ifndef TMSERHEAD_H
#define	TMSERHEAD_H


enum { kTMSERHEADSize = 17*sizeof(int) };

union TMSERHEAD
{
  char buffer[kTMSERHEADSize];
  struct
  {
    int    np;
    int    ne;
    int    first;
    int    last;
    int    step;
    int    startTime[5];
    int    deltaTime[5];
    int    vcomp;                  // number of data components in time series file
    int    vdata;
  };
};

#endif	/* TMSERHEAD_H */

