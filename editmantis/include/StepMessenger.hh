#ifndef StepMessenger_h
#define StepMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ApplicationState.hh"
#include "G4UIcmdWithAString.hh"
#include "SteppingAction.hh"
#include "G4UIdirectory.hh"

class SteppingAction;
class G4UIcmdWithAString;
class G4UIdirectory;


class StepMessenger: public G4UImessenger
{
public:
  StepMessenger(SteppingAction*);
  ~StepMessenger();

  void SetNewValue(G4UIcommand*, G4String);
private:
  SteppingAction* stepA;
  G4UIcmdWithAString* Cmd;
  G4UIdirectory *myDir;
};

#endif
