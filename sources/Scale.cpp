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
#include "Shape.h"
#include "Type.h"
#include "Node.h"
#include "Elem.h"
#include "Section.h"
#include "Model.h"
#include "Project.h"

#include "Scale.h"


SCALE::SCALE()
{
  this->type    = 0;
  this->lScale  =
  this->hScale  =
  this->kScale  =
  this->vScale  =
  this->tScale  =
  this->vtScale = 1.0;
}


SCALE::~SCALE()  { }


int    SCALE::gettype()    { return type; };
double SCALE::getlScale()  { return lScale; };
double SCALE::gethScale()  { return hScale; };
double SCALE::getkScale()  { return kScale; };
double SCALE::gettScale()  { return tScale; };
double SCALE::getvScale()  { return vScale; };
double SCALE::getvtScale() { return vtScale; };


void SCALE::init( int    type,
                  double lScale,
                  double hScale,
                  double kScale )
{
  this->type = type;

  switch( type )
  {
    case 2:
      this->lScale  = lScale;                    // scale of length
      this->hScale  = hScale;                    // scale of height
      this->kScale  = kScale;                    // scale of roughness height

      this->vScale  = sqrt( hScale );            // scale of velocities
      this->tScale  = lScale / sqrt( hScale );   // scale of time (horizontal)
      this->vtScale = hScale * hScale            // scale of viscosity
                      / sqrt( lScale );
//    this->vtScale = vScale * vScale * tScale;  // scale of viscosity (alternate)
      break;

    default:
      this->lScale  =
      this->hScale  =
      this->kScale  =
      this->vScale  =
      this->tScale  =
      this->vtScale = 1.0;
      break;
  }
}


void SCALE::init( int    type,
                  double lScale,
                  double vScale )
{
  this->type = type;

  switch( type )
  {
    case 1:
      this->lScale  = lScale;                    // scale of length
      this->hScale  = lScale;                    // scale of height
      this->kScale  = lScale;                    // scale of roughness height

      this->vScale  = vScale;                    // scale of velocities
      this->tScale  = lScale / vScale;           // scale of time
      this->vtScale = vScale * lScale;           // scale of viscosity
      break;

    default:
      this->lScale  =
      this->hScale  =
      this->kScale  =
      this->vScale  =
      this->tScale  =
      this->vtScale = 1.0;
      break;
  }
}


void SCALE::scale( PROJECT*   project )
{
  static int firstCall = true;

  switch( this->type )
  {
    case 1:
    case 2:   break;

    default:  return;
  }


  MODEL*   model    = project->M2D;
  GRID*    rg       = model->region;
  DRYREW*  dryRew   = &rg->dryRew;
  TIMEINT* timeint  = &project->timeint;
  BCONSET* bconSet  = timeint->bconSet;
  int      nSection = project->nSection;
  SECTION* section  = project->section;


  char text[200];

  if( firstCall )
  {
    REPORT::rpt.OutputLine1( 2 );

    sprintf( text, "  scaling of length:    %9.3lf\n", lScale );   REPORT::rpt.Output( text, 2 );
    sprintf( text, "             height:    %9.3lf\n", hScale );   REPORT::rpt.Output( text, 2 );
    sprintf( text, "             time:      %9.3lf\n", tScale );   REPORT::rpt.Output( text, 2 );
    sprintf( text, "             velocity:  %9.3lf\n", vScale );   REPORT::rpt.Output( text, 2 );
    sprintf( text, "             viscosity: %9.3lf\n", vtScale );  REPORT::rpt.Output( text, 2 );
  }

  for( int i=0; i<rg->Getnp(); i++ )
  {
    NODE* nd = rg->Getnode(i);

    nd->x       /= lScale;
    nd->y       /= lScale;
    nd->z       /= hScale;
    nd->zor     /= hScale;

    nd->dz      /= hScale;

    nd->vt      /= vtScale;
    nd->exx     /= vtScale;
    nd->exy     /= vtScale;
    nd->eyy     /= vtScale;
    nd->uu      /= vScale*vScale;
    nd->uv      /= vScale*vScale;
    nd->vv      /= vScale*vScale;

    nd->v.U     /= vScale;
    nd->v.V     /= vScale;
    nd->v.S     /= hScale;

    nd->v.K     /= vScale*vScale;
    nd->v.D     /= vScale*vScale/tScale;

    nd->v.C     /= 1.0;

    nd->v.dUdt  /= vScale/tScale;
    nd->v.dVdt  /= vScale/tScale;
    nd->v.dSdt  /= hScale/tScale;

    nd->v.dKdt  /= vScale*vScale/tScale;
    nd->v.dDdt  /= vScale*vScale/tScale/tScale;

    nd->vo.U    /= vScale;
    nd->vo.V    /= vScale;
    nd->vo.S    /= hScale;

    nd->vo.K    /= vScale*vScale;
    nd->vo.D    /= vScale*vScale/tScale;

    nd->vo.C    /= 1.0;

    nd->vo.dUdt /= vScale/tScale;
    nd->vo.dVdt /= vScale/tScale;
    nd->vo.dSdt /= hScale/tScale;

    nd->vo.dKdt /= vScale*vScale/tScale;
    nd->vo.dDdt /= vScale*vScale/tScale/tScale;
  }


  if( firstCall )
  {
    REPORT::rpt.OutputLine1( 2 );
    sprintf( text, "* element material...\n*\n" );
    REPORT::rpt.Output( text, 2 );

    sprintf( text, "%s\n*\n",
                   "*  no   typeB   rB   typeS   rS      dp1     sp1      hp2     dp2     sp2" );
    REPORT::rpt.Output( text, 2 );
  }

//ksScale = hScale * pow(10.0,-sqrt(lScale/hScale)/2.03);

  for( int i=0; i<TYPE::Get(); i++ )
  {
    TYPE* elt = TYPE::Getid( i+1 );

    switch( elt->rtype )
    {
      case TYPE::kCOWI:  elt->rcoef /= kScale;                        break;
      case TYPE::kMANN:  elt->rcoef /= pow(hScale,0.67)/sqrt(lScale); break;
      case TYPE::kCHEZ:  elt->rcoef /= sqrt(lScale/hScale);           break;
    }

    elt->dp[0] /= lScale;
    elt->sp[0] /= lScale;
    elt->hp    /= lScale;
    elt->dp[1] /= lScale;
    elt->sp[1] /= lScale;

    elt->vt   /= vtScale;

    if( firstCall )
    {
      sprintf( text, "  %3d   %1d %8.5lf   %7.4lf %7.4lf  %7.4lf %7.4lf %7.4lf\n",
               elt->no(0),
               elt->rtype, elt->rcoef,
               elt->dp[0], elt->sp[0], elt->hp, elt->dp[1], elt->sp[1] );
      REPORT::rpt.Output( text, 2 );
    }
  }


  if( firstCall )
  {
    REPORT::rpt.OutputLine1( 2 );
    sprintf( text, "%s\n*\n",
                   "*  no     vt        st    KDest    lm" );
    REPORT::rpt.Output( text, 2 );

    for( int i=0; i<TYPE::Get(); i++ )
    {
      TYPE* elt = TYPE::Getid( i+1 );

      sprintf( text, "  %3d  %9.6lf  %5.2lf  %5.2lf  %6.2lf\n",
               elt->no(0), elt->vt, elt->st, elt->KDest, elt->lm );
      REPORT::rpt.Output( text, 2 );
    }
  }


  project->timeint.deltaTime     /= tScale;
  project->timeint.relaxTimeFlow /= tScale;
  project->timeint.relaxTimeTurb /= tScale;

  project->hmin             /= hScale;

  project->maxDeltaUV       /= vScale;
  project->maxDeltaS        /= hScale;

  project->dep_minVt        /= vtScale;
  project->dep_minVtxx      /= vtScale;
  project->dep_minVtyy      /= vtScale;
  project->minVtKD          /= vtScale;
  project->minVtC           /= vtScale;
  project->minK             /= vScale*vScale;
  project->minD             /= vScale*vScale/tScale;
  project->maxUs            /= vScale;

  dryRew->dryLimit          /= hScale;
  dryRew->rewetLimit        /= hScale;


  if( firstCall )
  {
    REPORT::rpt.OutputLine1( 2 );
    sprintf( text, "* further values...\n*\n" );
    REPORT::rpt.Output( text, 2 );

    sprintf( text, "      hmin          = %9.3le\n", project->hmin );         REPORT::rpt.Output( text, 2 );

    sprintf( text, "      maxDeltaUV    = %9.3le\n", project->maxDeltaUV );   REPORT::rpt.Output( text, 2 );
    sprintf( text, "      maxDeltaS     = %9.3le\n", project->maxDeltaS );    REPORT::rpt.Output( text, 2 );

    sprintf( text, "      minVt (depr.) = %9.3le\n", project->dep_minVt );    REPORT::rpt.Output( text, 2 );
    sprintf( text, "      minVtKD       = %9.3le\n", project->minVtKD );      REPORT::rpt.Output( text, 2 );
    sprintf( text, "      minVtC        = %9.3le\n", project->minVtC );       REPORT::rpt.Output( text, 2 );
    sprintf( text, "      minK          = %9.3le\n", project->minK );         REPORT::rpt.Output( text, 2 );
    sprintf( text, "      minD          = %9.3le\n", project->minD );         REPORT::rpt.Output( text, 2 );
    sprintf( text, "      minC          = %9.3le\n", project->minC );         REPORT::rpt.Output( text, 2 );
    sprintf( text, "      maxUs         = %9.3le\n", project->maxUs );        REPORT::rpt.Output( text, 2 );

    sprintf( text, "      dryLimit      = %9.3le\n", dryRew->dryLimit );      REPORT::rpt.Output( text, 2 );
    sprintf( text, "      rewetLimit    = %9.3le\n", dryRew->rewetLimit );    REPORT::rpt.Output( text, 2 );
  }


  if( firstCall )
  {
    REPORT::rpt.OutputLine1( 2 );
    sprintf( text, "* %d section data...\n*\n", nSection );
    REPORT::rpt.Output( text, 2 );
  }

  for( int i=0; i<nSection; i++ )
  {
    section[i].scaleFR( lScale, hScale );

    if( firstCall )
    {
      sprintf( text, "  %10.4lf %10.4lf  %10.4lf %10.4lf  %10.4lf\n",
                     section[i].getxs(), section[i].getys(),
                     section[i].getxe(), section[i].getye(),
                     section[i].z );
      REPORT::rpt.Output( text, 2 );
    }
  }


  timeint->deltaTime     /= tScale;
  timeint->relaxTimeFlow /= tScale;
  timeint->relaxTimeTurb /= tScale;

  if( firstCall )
  {
    REPORT::rpt.Output( "\n", 2 );
    sprintf( text, "* scaled time...\n" );          REPORT::rpt.Output( text, 2 );

    sprintf( text, "    deltaTime:    %s\n",
                   timeint->deltaTime.Get() );      REPORT::rpt.Output( text, 2 );
    sprintf( text, "    relaxTimeFlow: %s\n",
                   timeint->relaxTimeFlow.Get() );  REPORT::rpt.Output( text, 2 );
    sprintf( text, "    relaxTimeTurb: %s\n",
                   timeint->relaxTimeTurb.Get() );  REPORT::rpt.Output( text, 2 );
  }


  for( int i=0; i<timeint->setsOfBcon; i++ )
  {
    if( firstCall )
    {
      REPORT::rpt.OutputLine1( 2 );
      sprintf( text, "* boundary condition set: %d\n", i + 1);
      REPORT::rpt.Output( text, 2 );
    }

    for( int j=0; j<bconSet[i].nofLines; j++ )
    {
      BCONLINE* bcLine = &(bconSet[i].bcLine[j]);

      if( firstCall )
      {
        sprintf( text, "*\n* boundary line: %d...\n*\n", bconSet[i].bcLine[j].no );
        REPORT::rpt.Output( text, 2 );
        sprintf( text, "%s%s\n",
                       "*      X         Y         qX        qY ",
                       "       S         K         D         C\n" );
        REPORT::rpt.Output( text, 2 );
      }

      for( int k=0; k<bcLine->ndat; k++ )
      {
        bcLine->x[k] /= lScale;
        bcLine->y[k] /= lScale;

        if( isFS(bcLine->kind, BCON::kQInlet) )
        {
          bcLine->U[k] /= lScale*hScale*vScale;
          bcLine->V[k] /= lScale*hScale*vScale;
        }

        else if( isFS(bcLine->kind, BCON::kInlet) )
        {
          bcLine->U[k] /= hScale*vScale;
          bcLine->V[k] /= hScale*vScale;
        }

        else if( isFS(bcLine->kind, BCON::kSetUV) )
        {
          bcLine->U[k] /= vScale;
          bcLine->V[k] /= vScale;
        }


        bcLine->S[k] /= hScale;
        bcLine->K[k] /= vScale*vScale;
        bcLine->D[k] /= vScale*vScale/tScale;


        if( isFS(bcLine->kind, BCON::kSetC) )
        {
          bcLine->C[k] /= 1.0;
        }

        else if( isFS(bcLine->kind, BCON::kRateC) )
        {
          bcLine->C[k] /= 1.0/lScale/hScale/vScale;
        }

        if( firstCall )
        {
          sprintf( text, "  %9.3le %9.3le %9.3le %9.3le %9.3le %9.3le %9.3le %9.3le\n",
                         bcLine->x[k], bcLine->y[k],
                         bcLine->U[k], bcLine->V[k], bcLine->S[k],
                         bcLine->K[k], bcLine->D[k], bcLine->C[k] );
          REPORT::rpt.Output( text, 2 );
        }
      }
    }


    if( firstCall  &&  bconSet[i].nofNodes > 0 )
    {
      REPORT::rpt.OutputLine1( 2 );
      sprintf( text, "* boundary nodes...\n*\n" );
      REPORT::rpt.Output( text, 2 );
      sprintf( text, "%s%s\n",
                     "*    no         qX        qY  ",
                     "      S         K         D         C\n" );
      REPORT::rpt.Output( text, 2 );
    }

    for( int j=0; j<bconSet[i].nofNodes; j++ )
    {
      BCON* bcNode = &(bconSet[i].bcNode[j]);

      if( isFS(bcNode->kind, BCON::kInlet) )
      {
        bcNode->val->U /= hScale*vScale;
        bcNode->val->V /= hScale*vScale;
      }

      else if( isFS(bcNode->kind, BCON::kSetUV) )
      {
        bcNode->val->U /= vScale;
        bcNode->val->V /= vScale;
      }


      bcNode->val->S /= hScale;
      bcNode->val->K /= vScale*vScale;
      bcNode->val->D /= vScale*vScale/tScale;


      if( isFS(bcNode->kind, BCON::kSetC) )
      {
        bcNode->val->C[0] /= 1.0;
      }

      else if( isFS(bcNode->kind, BCON::kRateC) )
      {
        bcNode->val->C[0] /= 1.0/lScale/hScale/vScale;
      }


      if( firstCall )
      {
        sprintf( text, " %5d %9.3le %9.3le %9.3le %9.3le %9.3le %9.3le\n",
                       bcNode->no,
                       bcNode->val->U, bcNode->val->V, bcNode->val->S,
                       bcNode->val->K, bcNode->val->D, bcNode->val->C[0] );
        REPORT::rpt.Output( text, 2 );
      }
    }
  }

  if( firstCall )
  {
    REPORT::rpt.OutputLine1( 2 );
    REPORT::rpt.Output( "\n\n", 2  );
  }

  firstCall = false;
}
