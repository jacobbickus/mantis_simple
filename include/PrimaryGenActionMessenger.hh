#ifndef PrimaryGenActionMessenger_h
#define PrimaryGenActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ApplicationState.hh"
#include "G4UIcmdWithADouble.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"

class PrimaryGeneratorAction;
class G4UIcmdWithADouble;
class G4UIdirectory;


class PrimaryGenActionMessenger: public G4UImessenger
{
public:
  PrimaryGenActionMessenger(PrimaryGeneratorAction*);
  ~PrimaryGenActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);
private:
  PrimaryGeneratorAction* genA;
  G4UIcmdWithADouble* Cmd;
  G4UIdirectory *myDir;
};

#endif
