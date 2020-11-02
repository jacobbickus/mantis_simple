#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "g4root.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


class HistoManager
{
public:
HistoManager();
~HistoManager();

void finish();     // close root file
void Book();

private:
G4bool fFactoryOn;
};

#endif
