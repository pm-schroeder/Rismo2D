// ======================================================================================
//                                       T I M E
// ======================================================================================
// This class implements the time variable.
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
// 05.06.2005     sc     first implementation
// 07.12.2007     sc     change to iso time format: ddddddThh:mm:ss.ssssss
// 27.06.2011     sc     get and set time array for binary output / input
//                       void TIME::Get( int[5] ); void TIME::Set( int[5] );
//
// ======================================================================================

// read this header file only once ...

#ifndef TIME_INCL
#define TIME_INCL

class TIME
{
  public:
    const static TIME null;

  private:
    int micro;                // micro seconds (0 ... 1.000.000)
    int s;                    // seconds       (0 ... 60)
    int m;                    // minutes       (0 ... 60)
    int h;                    // hours         (0 ... 24)
    int d;                    // days
    char time[23];

  public:
    TIME()
    {
      Init();
    };

    //////////////////////////////////////////////////////////////////////////////////////

    TIME operator =( const TIME& t )
    {
      micro = t.micro;
      s     = t.s;
      m     = t.m;
      h     = t.h;
      d     = t.d;
      return *this;
    };

    TIME operator +( const TIME& t )
    {
      TIME r = *this;
      r.micro += t.micro;
      r.s     += t.s;
      r.m     += t.m;
      r.h     += t.h;
      r.d     += t.d;
      r.Reduce();
      return r;
    };

    TIME operator +=( const TIME& t )
    {
      this->micro += t.micro;
      this->s     += t.s;
      this->m     += t.m;
      this->h     += t.h;
      this->d     += t.d;
      this->Reduce();
      return *this;
    };

    TIME operator -( const TIME& t )
    {
      TIME r = *this;
      r.micro -= t.micro;
      r.s     -= t.s;
      r.m     -= t.m;
      r.h     -= t.h;
      r.d     -= t.d;
      r.Reduce();
      return r;
    };

    TIME operator -=( const TIME& t )
    {
      this->micro -= t.micro;
      this->s     -= t.s;
      this->m     -= t.m;
      this->h     -= t.h;
      this->d     -= t.d;
      this->Reduce();
      return *this;
    };

    TIME operator *( const double d )
    {
      TIME r = *this;
      double sec = r.Getsec() * d;
      r.Setsec( sec );
      return r;
    };

    TIME operator *=( const double d )
    {
      double sec = this->Getsec() * d;
      this->Setsec( sec );
      return *this;
    };

    TIME operator /( const double d )
    {
      TIME r = *this;
      double sec = r.Getsec() / d;
      r.Setsec( sec );
      return r;
    };

    TIME operator /=( const double d )
    {
      double sec = this->Getsec() / d;
      this->Setsec( sec );
      return *this;
    };
    //////////////////////////////////////////////////////////////////////////////////////

    int operator <( const TIME& t )
    {
      if(   (d <t.d)
         || (d==t.d && h <t.h)
         || (d==t.d && h==t.h && m <t.m)
         || (d==t.d && h==t.h && m==t.m && s <t.s)
         || (d==t.d && h==t.h && m==t.m && s==t.s && micro<t.micro) ) return true;
      else                                                            return false;
    };

    int operator >( const TIME& t )
    {
      if(   (d >t.d)
         || (d==t.d && h >t.h)
         || (d==t.d && h==t.h && m >t.m)
         || (d==t.d && h==t.h && m==t.m && s >t.s)
         || (d==t.d && h==t.h && m==t.m && s==t.s && micro>t.micro) ) return true;
      else                                                            return false;
    };

    int operator ==( const TIME& t )
    {
      if( d==t.d && h==t.h && m==t.m && s==t.s && micro==t.micro ) return true;
      else                                                         return false;
    };

    int operator !=( const TIME& t )
    {
      if( d!=t.d || h!=t.h || m!=t.m || s!=t.s || micro!=t.micro ) return true;
      else                                                         return false;
    };

    int operator <=( const TIME& t )
    {
      if( *this == t  ||  *this < t) return true;
      else                           return false;
    };

    int operator >=( const TIME& t )
    {
      if( *this == t  ||  *this > t) return true;
      else                           return false;
    };

    //////////////////////////////////////////////////////////////////////////////////////

    void Set( int t[5] )
    {
      Init();

      micro = t[0];
      s     = t[1];
      m     = t[2];
      h     = t[3];
      d     = t[4];

      Reduce();
    }

    //////////////////////////////////////////////////////////////////////////////////////

    void Get( int t[5] )
    {
      t[0] = micro;
      t[1] = s;
      t[2] = m;
      t[3] = h;
      t[4] = d;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    void Set( const char* t )
    {
      Init();

      strncpy( time, t, 22 );
      time[22] = '\0';

      int ncolon = 0;
      int nT     = 0;

      for( char* c=time; *c; c++ )
      {
        if( *c == ':' ) ncolon++;
        if( *c == 'T' ) nT++;
      }

      double sec;
      int    val[4] = { 0, 0, 0, 0 };

      if( nT == 0 )
      {
        switch( ncolon )
        {
          default:
          case 0:
            sscanf( time, "%lf", &sec );
            break;

          case 1:
            sscanf( time, "%d:%lf", &val[2], &sec );
            break;

          case 2:
            sscanf( time, "%d:%d:%lf", &val[1], &val[2], &sec );
            break;

          case 3:
            sscanf( time, "%d:%d:%d:%lf", &val[0], &val[1], &val[2], &sec );
            break;
        }
      }
      else
      {
        sscanf( time, "%dT%d:%d:%lf", &val[0], &val[1], &val[2], &sec );
      }

      val[3] = (int)(sec);
      micro  = (int)(1000000*sec - 1000000*val[3]);

      d = val[0];
      h = val[1];
      m = val[2];
      s = val[3];
/*
      int   n = 0;
      int   val[4];
      char* token;

      token = strtok( time, ":" );

      while( token  &&  n < 4 )
      {
        char* pos = strchr( token, '.' );

        if( pos )
        {
          double sec;
          sscanf( token, "%lf", &sec );

          val[n] = (int)(sec);
          micro  = (int)(1000000*sec - 1000000*val[n]);
          n++;
          break;
        }
        else
        {
          sscanf( token, "%d", &val[n] );
          n++;
        }

        token = strtok( NULL, ":" );
      }

      n--;  if( n >= 0 )  s = val[n];
      n--;  if( n >= 0 )  m = val[n];
      n--;  if( n >= 0 )  h = val[n];
      n--;  if( n >= 0 )  d = val[n];
*/
      Reduce();
    };

    char* Get()
    {
      if( d > 0 )  sprintf( time, "%dT%02d:%02d:%02d.%06d", d, h, m, s, micro );
      else         sprintf( time, "%02d:%02d:%02d.%06d", h, m, s, micro );
      return time;
    };

    void Setsec( double sec )
    {
      Init();

      s     = (int)(sec);
      micro = (int)(1000000.0 * (sec - s));

      Reduce();
    }

    double Getsec()
    {
      return   (double)d * 86400.0
             + (double)h * 3600.0
             + (double)m * 60.0
             + (double)s
             + (double)micro / 1000000.0;
    }

  protected:
    void Init()
    {
      micro = s = m = h = d = 0;
    };

    void Reduce()
    {
      int again = false;

      if( micro >= 1000000 )  { s += micro / 1000000;  micro %= 1000000; }
      if(     s >= 60      )  { m +=     s / 60;           s %= 60;      }
      if(     m >= 60      )  { h +=     m / 60;           m %= 60;      }
      if(     h >= 24      )  { d +=     h / 24;           h %= 24;      }

      if( micro < 0 )  { micro += 1000000;  s--;  again = true; }
      if(     s < 0 )  {     s += 60;       m--;  again = true; }
      if(     m < 0 )  {     m += 60;       h--;  again = true; }
      if(     h < 0 )  {     h += 24;       d--;  again = true; }

      if( again )  Reduce();
    };
};

#endif
