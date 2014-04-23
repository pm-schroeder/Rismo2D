// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// S T A T I S T
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Statist.h   : definition file of the class.
// Statist.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// This class implements statistical functions for analysis of transient flow.
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

#ifndef STATIST_INCL
#define STATIST_INCL

class NODE;
class ELEM;
class SUBDOM;


class STATIST
{
  private:
    int     np;        // number of values (nodes)

    int     n;         // counter for summation
    int    *nwet;      // counter for summation of values at wet nodes

    double *U;         // sum of U
    double *V;         //        V
    double *S;         //        S
    double *H;         //        H
    double *Vt;        //        Vt

    double *UU;        //        U*U
    double *UV;        //        U*V
    double *VV;        //        V*V
    double *HH;        //        H*H
    double *VtVt;      //        Vt*Vt

  public:
    STATIST();
    ~STATIST();

    double GetMeanU( int no );
    double GetMeanV( int no );
    double GetMeanS( int no );
    double GetMeanH( int no );
    double GetMeanUs( int no );
    double GetMeanVt( int no );
    double GetVarU( int no );
    double GetVarUV( int no );
    double GetVarV( int no );
    double GetKinE( int no );
    double GetSdevH( int no );
    double GetVarVt( int no );
    double GetVtVt( int no );
    double GetFldRate( int no );

    void Init( int np );
    void Read( int np, char *fileName, SUBDOM *subdom );
    void Write( MODEL *model, int release, char *staFile,
                char *rgFile, int timeStep, PROJECT *project );
    void Sum( MODEL* );
    void Reset( MODEL* );
};

#endif
