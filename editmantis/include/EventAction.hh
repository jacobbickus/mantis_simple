#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include "G4Types.hh"
#include "G4ParticleGun.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"
#include "HistoManager.hh"
#include "EventMessenger.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class G4Event;
class EventMessenger;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(G4ParticleGun*,SteppingAction*);
   ~EventAction();

  public:
    void BeginOfEventAction(const G4Event* anEvent);
    void EndOfEventAction(const G4Event* anEvent);

    void SetNumPhotonsFlag(G4int val){drawNumPhotonsFlag = val;};

private:

    void ResetEverything();
    //void Weighting(const G4Event* anEvent);
    G4int drawNumPhotonsFlag;
    G4ParticleGun* particle_gun_local;
    SteppingAction* stepA_local;
    EventMessenger* eventM;


};

#endif
