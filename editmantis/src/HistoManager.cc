#include "HistoManager.hh"

extern G4String gOutName;
HistoManager::HistoManager():checkPrint(false), fFactoryOn(false), checkLast(false)
{}

HistoManager::~HistoManager()
{}

void HistoManager::Book()
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager->SetVerboseLevel(0);
  manager->SetNtupleMerging(true);

  // Create Directories
  //manager->SetNtupleDirectoryName("CherenkovData");

  // open output file
  G4bool fileOpen = manager->OpenFile(gOutName);

  if(! fileOpen){
    G4cerr << "HistoManager::Book(): Cannot Open " <<manager->GetFileName()<<G4endl;
    return;
  }

  // Create Ntuple for water data
    manager->CreateNtuple("WaterData","WaterData");

    manager->CreateNtupleDColumn("Energy");
    manager->CreateNtupleIColumn("isNRF");
    manager->CreateNtupleDColumn("x");
    manager->CreateNtupleDColumn("y");
    manager->CreateNtupleDColumn("z");
    manager->CreateNtupleDColumn("Theta");
    manager->CreateNtupleDColumn("Time");
    manager->FinishNtuple();

    // Create Ntuple for Incident Detector Data

    manager->CreateNtuple("IncidentDet","IncidentDet");

    manager->CreateNtupleDColumn("Energy");
    manager->CreateNtupleIColumn("isNRF");
    manager->CreateNtupleDColumn("x");
    manager->CreateNtupleDColumn("y");
    manager->CreateNtupleDColumn("z");
    manager->FinishNtuple();

    // Create Ntuple for Reactions within detector

    manager->CreateNtuple("DetProcess","DetProcess");

    manager->CreateNtupleSColumn("Process");
    manager->CreateNtupleDColumn("E_beam");
    manager->FinishNtuple();

    // Create Nutple for Energy if detected
    manager->CreateNtuple("Detected","Detected");
    manager->CreateNtupleDColumn("Energy");
    manager->CreateNtupleIColumn("isNRF");
    manager->CreateNtupleDColumn("x");
    manager->CreateNtupleDColumn("y");
    manager->CreateNtupleDColumn("z");
    manager->FinishNtuple();

    // Create Ntuple for Number of Photons Emitted in Water
    manager->CreateNtuple("NumPhotons","NumPhotons");
    manager->CreateNtupleDColumn("E_Beam");
    manager->CreateNtupleDColumn("NumPhotons");
    manager->FinishNtuple();

    // Create Ntuple for Interogation Object
    manager->CreateNtuple("IntObj","IntObj");
    manager->CreateNtupleDColumn("E_incident");
    manager->CreateNtupleIColumn("isNRF");
    manager->FinishNtuple();

    // Create Ntuple for incident water data 
    manager->CreateNtuple("IncidentWater","IncidentWater");
    manager->CreateNtupleDColumn("E_incident");
    manager->CreateNtupleIColumn("isNRF");
    manager->FinishNtuple();

#if defined (G4ANALYSIS_USE_ROOT)
    // Create Ntuple for Resonance Sampling
    manager->CreateNtuple("ResSample", "ResSample");
    manager->CreateNtupleIColumn("Res");
    manager->FinishNtuple();
#endif

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

G4bool HistoManager::OnceAWhileSave(time_t val)
{
  startTime = val;
  //std::cout << "here" <<std::endl;
  time_t Now = time(0);
  G4int testVal = (Now - startTime)%60;
  if(!checkLast && testVal > 1)
  {
    checkLast = true;
  }
  //std::cout << std::endl << testVal << std::endl;
  if((checkLast && (testVal == 1 || testVal == 0))) // every 60 seconds
  {
    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    manager->Write();
    checkLast = false;
    checkPrint = true;
    G4cout << "Data Saved on Event..." << G4endl;
  }
  else
  {
    checkPrint = false;
  }
  return checkPrint;
}
