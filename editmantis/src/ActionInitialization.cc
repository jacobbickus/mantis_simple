#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"


ActionInitialization::ActionInitialization(const DetectorConstruction* det)
 : G4VUserActionInitialization(), fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{

    HistoManager* histo = new HistoManager();
    PrimaryGeneratorAction *primary_action=new PrimaryGeneratorAction();
    SetUserAction(primary_action);
    RunAction* run = new RunAction(histo);
    SetUserAction(run);

    SteppingAction *stepAction = new SteppingAction(fDetector,
      primary_action->GetParticleGun(), run);
    SetUserAction(stepAction);
    EventAction *eventAction = new EventAction(primary_action->GetParticleGun(),
    stepAction);
    SetUserAction(eventAction);
    StackingAction *stackAction = new StackingAction(fDetector, run);
    SetUserAction(stackAction);
}
