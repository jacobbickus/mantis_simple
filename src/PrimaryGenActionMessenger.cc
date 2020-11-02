#include "PrimaryGenActionMessenger.hh"


PrimaryGenActionMessenger::PrimaryGenActionMessenger(PrimaryGeneratorAction* genAction)
:genA(genAction)
{
  myDir = new G4UIdirectory("/input/");
  myDir->SetGuidance("Input Commands");
  Cmd = new G4UIcmdWithADouble("/input/energy",this);
  Cmd->SetGuidance("Choose Desired Energy Input");
  Cmd->SetParameterName("energy",false);
  Cmd->SetDefaultValue(-1);

}

PrimaryGenActionMessenger::~PrimaryGenActionMessenger()
{
  delete Cmd;
}


void PrimaryGenActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ // this function only works if named SetNewValue

  if(command == Cmd)
  {
    G4double theCommand = Cmd->GetNewDoubleValue(newValue);
    genA->SetEnergyValue(theCommand);

  }
  else
  {
    G4cerr << "ERROR PrimaryGenActionMessenger :: SetStepInputValue command != Cmd" << G4endl;
  }
}
