#include "StepMessenger.hh"


StepMessenger::StepMessenger(SteppingAction* stepAction)
:stepA(stepAction)
{
  myDir = new G4UIdirectory("/output/");
  myDir->SetGuidance("Output Commands");
  Cmd = new G4UIcmdWithAString("/output/myoutput",this);
  Cmd->SetGuidance("Choose Desired Outputs");
  Cmd->SetGuidance("Choice: WaterData, WaterIncidentData DetIncidentData, DetData, IntObjData, none (default)");
  Cmd->SetParameterName("choice",false);
  Cmd->SetDefaultValue("none");
  Cmd->SetCandidates("WaterData WaterIncidentData DetIncidentData DetData IntObjData none");
  //Cmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle, G4State_EventProc);

}

StepMessenger::~StepMessenger()
{
  delete Cmd;
}


void StepMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ // this function only works if named SetNewValue
  if(command == Cmd)
  {
    G4String theCommand = newValue;
    if(theCommand == "WaterData")
    {
      stepA->SetWaterDataFlag(1);
    }
    else if(theCommand == "WaterIncidentData")
    {
      stepA->SetIncWatDataFlag(1);
    }
    else if(theCommand == "DetIncidentData")
    {
      stepA->SetIncidentDataFlag(1);
    }
    else if(theCommand == "DetData")
    {
      stepA->SetDetDataFlag(1);
    }
    else if(theCommand == "IntObjData")
    {
      stepA->SetIntObjDataFlag(1);
    }
    else if(theCommand == "none")
    {
      G4cout << "No output requested." << G4endl;
    }
    else
    {
      G4cerr << "ERROR: StepMessenger::SetStepInputValue Command not found." << G4endl;
    }

  }
  else
  {
    G4cerr << "ERROR StepMessenger :: SetStepInputValue command != Cmd" << G4endl;
  }
}
