#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "g4root.hh"
#include <time.h>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


class HistoManager
{
public:
HistoManager();
~HistoManager();

void finish();     // close root file
void Book();
G4bool OnceAWhileSave(time_t);
G4bool checkPrint;

private:
G4bool fFactoryOn;
G4bool checkLast;
time_t startTime;

};

#endif
