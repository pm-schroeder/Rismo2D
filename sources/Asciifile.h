// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// A S C I I F I L E
//
// /////////////////////////////////////////////////////////////////////////////////////////////////
//
// FILES
//
// Asciifile.h   : definition file of the class.
// Asciifile.cpp : implementation file of the class.
//
// -------------------------------------------------------------------------------------------------
//
// DESCRIPTION
//
// Supports a convenient way to read and write Ascii files.
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
//  01.01.1998     sc     first implementation / first concept
//
// /////////////////////////////////////////////////////////////////////////////////////////////////

#  ifndef ASCIIFILE_INCL
#  define ASCIIFILE_INCL

// ======================================================================================
// Include-Dateien...

#  include <stdio.h>

// ======================================================================================

class ASCIIFILE
{
  // ====================================================================================
  //                               V A R I A B L E N
  // ====================================================================================

  public:
  protected:

  private:
    FILE* id;            // Pointer auf den geoeffneten Stream
    char* name;          // Dateiname
    char* line;          // aktuell eingelesene Zeile
    int   length;        // maximale Laenge einer Zeie (default =1000)
    int   lineNo;        // aktuelle Zeilennummer

  // ====================================================================================
  //                                M E T H O D E N
  // ====================================================================================

  public:
    // ----------------------------------------------------------------------------------
    // Constructor der Klasse. Der Constructor oeffnet die Datei <fileName> entsprechend
    // der Vorgabe in <mode> zum Lesen, Schreiben, ...
    // <mode> ist identisch mit dem entsprechenden Parameter in der Funktion FILE* fopen()
    // (z.B.: "w" fuer Schreiben und "r" fuer Lesen).
    // <maxLineLength> ist ein optionaler Parameter, der die maximale Zeilenlaenge aendert;
    // die Zeilenlaenge hat einen Standardwert von 1000.
    ASCIIFILE( const char* fileName, const char* mode, int maxLineLength=1000 );

    // ----------------------------------------------------------------------------------
    // Destructor der Klasse. Schliesst den Stream und gibt Speicher frei.
    ~ASCIIFILE();

    // ----------------------------------------------------------------------------------
    // Liefert den Streampointer <id>. Diese Methode wird im allgemeinen zum Ueberpruefen
    // verwendet, ob die Datei erfolgreich geoeffnet wurde.
    FILE* getid();

    // ----------------------------------------------------------------------------------
    // Liefert die aktuelle Zeilennummer.
    int   getlineNo();

    // ----------------------------------------------------------------------------------
    // Liefert den Namen <name> der geoeffneten Datei.
    char* getname();

    // ----------------------------------------------------------------------------------
    // Spult an den Dateianfang zurueck (<lineNo> = 0).
    void  rewind();

    // ----------------------------------------------------------------------------------
    char* next();

    // ----------------------------------------------------------------------------------
    // Liest die naechste Zeile aus der Datei und gibt einen Pointer auf die eingelesene
    // Zeichenkette zurueck. Zeilen, die in der ersten Spalte ein "*" oder "#" haben, sind
    // Kommentarzeilen und werden ueberlesen.
    char* nextLine();

    // ----------------------------------------------------------------------------------
    // Schreibt auf die geoeffnete Datei. Die Methode arbeitet wie die Funktion printf()
    // aus der Standardbibliothek. Gemaess den Vorgaben in <format> werden weitere
    // Parameter erwartet.
    void  print( const char* format, ... );

    // ----------------------------------------------------------------------------------
    // Gleiche Methode wie <print(char* format,. ..)>. Zusaetzlich erfolgt jedoch eine
    // Ausgabe auf den Bildschirm.
    void  echoPrint( const char* format, ... );

  // ====================================================================================
  private:

  // ====================================================================================
  protected:
};

#endif
