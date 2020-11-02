#include "SteppingAction.hh"
#include "G4Cerenkov.hh"

extern G4bool output;

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{
}

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    // Run Logical Checks
    if(aStep->GetPostStepPoint()==NULL){
        return; // at the end of the world
    }
    else if(aStep->GetPostStepPoint()->GetPhysicalVolume()==NULL){
        return;
    }

    G4Track* theTrack = aStep->GetTrack();

     // testing brem target
     if(aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName().compare(0, 6 ,"Target") == 0
        && aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName().compare(0, 6, "Target") != 0 && output 
        && theTrack->GetParticleDefinition() == G4Gamma::Definition())
        {
          G4double energy_test = theTrack->GetKineticEnergy()/(MeV);
          manager->FillNtupleDColumn(0,0, energy_test);
          manager->AddNtupleRow(0);
        }

    
} // end of user steepping action function


bool SteppingAction::TrackMustDie(const G4Step *aStep){

    //Energy cuts
    if(aStep->GetTrack()->GetDefinition()->GetPDGEncoding()==EventGeneratorParticle &&
       aStep->GetTrack()->GetKineticEnergy()/MeV < LowEnergyCutoff){ //for the main track, kill it off below some limit
      return true;
    }

    return false;
}
