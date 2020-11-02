#ifndef EventMessenger_h
#define EventMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ApplicationState.hh"
#include "G4UIcmdWithAString.hh"
#include "EventAction.hh"
#include "G4UIdirectory.hh"

class EventAction;
class G4UIcmdWithAString;
class G4UIdirectory;


class EventMessenger: public G4UImessenger
{
public:
  EventMessenger(EventAction*);
  ~EventMessenger();

  void SetNewValue(G4UIcommand*, G4String);
private:
  EventAction* eventA;
  G4UIcmdWithAString* Cmd;
  G4UIdirectory *myDir;
};

#endif
