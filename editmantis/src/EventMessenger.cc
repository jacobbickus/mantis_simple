#include "EventMessenger.hh"


EventMessenger::EventMessenger(EventAction* event)
:eventA(event)
{
  myDir = new G4UIdirectory("/event/");
  myDir->SetGuidance("Output Commands");
  Cmd = new G4UIcmdWithAString("/event/output",this);
  Cmd->SetGuidance("Choose if NumPhotons Ntuple desired");
  Cmd->SetGuidance("Choice: Yes, No(default)");
  Cmd->SetParameterName("choice",false);
  Cmd->SetDefaultValue("No");
  Cmd->SetCandidates("Yes No");
  Cmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle, G4State_EventProc);

}

EventMessenger::~EventMessenger()
{
  delete Cmd;
}


void EventMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ // this function only works if named SetNewValue
  if(command == Cmd)
  {
    G4String theCommand = newValue;
    if(theCommand == "Yes")
    {
      eventA->SetNumPhotonsFlag(1);
    }

  }
  else
  {
    G4cerr << "ERROR StepMessenger :: SetStepInputValue command != Cmd" << G4endl;
  }
}
