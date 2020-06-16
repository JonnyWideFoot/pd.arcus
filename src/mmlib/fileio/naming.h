#ifndef __NAMING_H
#define __NAMING_H

// ----------------------------------------------------------------------------------------
// Naming.h
// DEPRECATED: Note this functionality should no longer be used!
// It has been replaced by functionality in sequence.h
// ----------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
// --- Naming constants ---
// ----------------------------------------------------------------------------------

//extern const char Nameconvention_AAOrder[22];
//extern char *PDB_Nameconvention[22][50];
//extern char *XPLOR_Nameconvention[22][50];
//extern char *PDAMBER_Nameconvention[22][50];

// --------------------------------------------------------------------------------
// Convertion from aminoacid numbers to letter and vice versa
// --------------------------------------------------------------------------------

char getAALetterFromFullName(char *resname);
char getAANumberFromFullName(char *resname);
char getAALetter(int aanum);
// The following are unused:
//const char *getAANameFull(int aanum);
//int getAANumber(char aachar);

#endif


