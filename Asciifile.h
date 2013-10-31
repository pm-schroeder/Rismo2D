// ======================================================================================
//                                A S C I I F I L E
// ======================================================================================
// Die Klasse unterstuetzt den Zugriff auf ASCII-Dateien. Die Klasse ermoeglicht das
// Oeffnen und Schliessen sowie das zeilenweise Lesen und Schreiben von ASCII-Dateien.
// ======================================================================================
//
// Copyright (C) 1998-2005  by  P.M. SCHROEDER
//
// All rights reserved.
//
// This source code is part of the RISMO2D modelling software.
// As long as you have no contract (Source Code License Agreement
// for the Rismo2D Software / Version 1.0 or any later version")
// with the copyright holder, you are NOT ALLOWED to make any
// changes to this source code.
//
// Datum           :  27.02.2005
//
// ======================================================================================
//
//   Datum               Erlaeuterung
// ----------   ------   ----------------------------------------------------------------
// 01.01.1998     sc     erste Programmierung
//
// ======================================================================================

// Damit die Header-Datei nur ein einmal gelesen wird...

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
