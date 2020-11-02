#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "HistoManager.hh"


ActionInitialization::ActionInitialization(const DetectorConstruction* det)
 : G4VUserActionInitialization()
{}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::Build() const
{

    HistoManager* histo = new HistoManager();
    SetUserAction(new PrimaryGeneratorAction());
    RunAction* run = new RunAction(histo);
    SetUserAction(run);
    SetUserAction(new SteppingAction());
}
