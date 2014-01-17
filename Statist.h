// ======================================================================================
//                                    S T A T I S T
// ======================================================================================
// This class implements statistical functions for analysis of transient flow.
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

    double *maxTau;    // maximum of bottom friction
    double *maxUs;     // maximum of velocity

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
    double GetMaxTau( int no );
    double GetMaxUs( int no );

    void Init( int np );
    void Read( int np, char *fileName, SUBDOM *subdom );
    void Write( MODEL *model, int release, char *staFile,
                char *rgFile, int timeStep, PROJECT *project );
    void Sum( PROJECT*, MODEL* );
    void Reset( MODEL* );
};

#endif
