#ifndef PCSensitiveDetector_h
#define PCSensitiveDetector_h 1
#include "G4VSensitiveDetector.hh"
#include "HistoManager.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "globals.hh"

class PCSensitiveDetector : public G4VSensitiveDetector
{
public:
  //constructor
  PCSensitiveDetector(G4String SDname);
  //destructor
  ~PCSensitiveDetector();

  G4bool ProcessHits(G4Step *step, G4TouchableHistory *hist);
  //void Initialize(G4HCofThisEvent* HCE);
  //void EndOfEvent(G4HCofThisEvent* HCE);
private:
  G4OpBoundaryProcessStatus fExpectedNextStatus;
  G4String procCount;
  G4double det_energy;
};

#endif
