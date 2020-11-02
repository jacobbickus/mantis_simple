#include "HistoManager.hh"

extern G4String gOutName;
HistoManager::HistoManager(): fFactoryOn(false)
{}

HistoManager::~HistoManager()
{}

void HistoManager::Book()
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetVerboseLevel(0);
  manager->SetNtupleMerging(true);

  // open output file
  G4bool fileOpen = manager->OpenFile(gOutName);

  if(! fileOpen){
    G4cerr << "HistoManager::Book(): Cannot Open " <<manager->GetFileName()<<G4endl;
    return;
  }

    manager->CreateNtuple("ChopperData","Brem Distribution on Chopper");
    manager->CreateNtupleDColumn("E_incident");
    manager->CreateNtupleIColumn("isGamma");
    manager->FinishNtuple();

    fFactoryOn = true;
    G4cout << "Data Book Created." << G4endl;
    G4cout << "Output file is open in " << manager->GetFileName()<<"."
          << manager->GetFileType() << G4endl;

}

void HistoManager::finish()
{
    if(! fFactoryOn){
        G4cout << "ERROR HistoManager::finish: Failed to write to file" << G4endl;
        return;
    }
    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    manager->Write();
    manager->CloseFile();
    G4cout << "Ntuples are saved. " << G4endl;

    delete manager;
    fFactoryOn = false;
}
