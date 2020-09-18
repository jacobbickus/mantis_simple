#include "PrimaryGenActionMessenger.hh"


PrimaryGenActionMessenger::PrimaryGenActionMessenger(PrimaryGeneratorAction* genAction)
:genA(genAction)
{
  myDir = new G4UIdirectory("/input/");
  myDir->SetGuidance("Input Commands");
  Cmd = new G4UIcmdWithADouble("/input/energy",this);
  Cmdres = new G4UIcmdWithAString("/input/ResData",this);
  Cmd->SetGuidance("Choose Desired Energy Input");
  Cmdres->SetGuidance("Choose if Resonance Sample Output desired");
  Cmd->SetParameterName("energy",false);
  Cmdres->SetParameterName("res",false);
  Cmd->SetDefaultValue(-1);
  Cmdres->SetCandidates("True true False false");
  Cmdres->SetDefaultValue("false");
  //Cmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);

}

PrimaryGenActionMessenger::~PrimaryGenActionMessenger()
{
  delete Cmd;
  delete Cmdres;
}


void PrimaryGenActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ // this function only works if named SetNewValue

  if(command == Cmd)
  {
    G4double theCommand = Cmd->GetNewDoubleValue(newValue);
    genA->SetEnergyValue(theCommand);

  }
  else if(command == Cmdres)
  {
    G4String theCommandRes = newValue;
    if(theCommandRes == "True" || theCommandRes == "true")
    {
      genA->SetResDataFlag(1);
    }
  }
  else
  {
    G4cerr << "ERROR PrimaryGenActionMessenger :: SetStepInputValue command != Cmd" << G4endl;
  }
}
