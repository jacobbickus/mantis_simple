#include "PCSensitiveDetector.hh"

extern G4bool output;

PCSensitiveDetector::PCSensitiveDetector(G4String SDname)
:G4VSensitiveDetector(SDname)
{
  G4String DetName = SDname;
  fExpectedNextStatus = Undefined;
}
PCSensitiveDetector::~PCSensitiveDetector()
{}

G4bool PCSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
  G4StepPoint* thePoint = aStep-> GetPreStepPoint();
  G4TouchableHandle touchable = thePoint->GetTouchableHandle();
  G4int copyNo = touchable->GetVolume(0)->GetCopyNo();
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  G4Track* theTrack = aStep->GetTrack();
  const G4DynamicParticle* theParticle = theTrack->GetDynamicParticle();

  G4OpBoundaryProcessStatus theStatus = Undefined;
  G4ProcessManager* OpManager =
    G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  G4int MAXofPostStepLoops =
    OpManager->GetPostStepProcessVector()->entries();
  G4ProcessVector* postStepDoItVector =
    OpManager->GetPostStepProcessVector(typeDoIt);
if(output){
    manager->FillNtupleDColumn(1,0,theParticle->GetKineticEnergy()/(MeV));
    manager->FillNtupleIColumn(1,1, copyNo);
    manager->AddNtupleRow(1);
  }


  for (G4int i=0; i<MAXofPostStepLoops; ++i)
  {
      G4VProcess* currentProcess = (*postStepDoItVector)[i];
      //std::cout << G4Cerenkov::PostStepGetPhysicalInteractionLength(theTrack) << std::endl;

      G4OpBoundaryProcess* opProc = dynamic_cast<G4OpBoundaryProcess*>(currentProcess);

      if(opProc && output){
          theStatus = opProc->GetStatus();

          if(theStatus == Transmission)
          {
              //run->AddTransmission();
              procCount = "Trans";
              // G4cout << "transmission"<< G4endl;
          }
          else if(theStatus == FresnelRefraction)
          {
              //run->AddFresnelRefraction();
              procCount = "Refr";
              //G4cout << "fres refraction"<< G4endl;
          }
          else if (theStatus == TotalInternalReflection)
          {
            //run->AddTotalInternalReflection();
              procCount = "Int_Refl";
              //G4cout << "totalinternal" << G4endl;
          }
          else if (theStatus == LambertianReflection)
          {
            //run->AddLambertianReflection();
              procCount = "Lamb";
          }
          else if (theStatus == LobeReflection)
          {
            //run->AddLobeReflection();
              procCount = "Lobe";
          }
          else if (theStatus == SpikeReflection)
          {
            //run->AddSpikeReflection();
              procCount = "Spike";
          }
          else if (theStatus == BackScattering)
          {
            //run->AddBackScattering();
              procCount = "BackS";
          }
          else if (theStatus == Absorption)
          {
            //run->AddAbsorption();
              procCount = "Abs";
          }
          else if (theStatus == Detection)
          {
            //run->AddDetection();
              procCount = "Det";
              det_energy = aStep->GetTotalEnergyDeposit();
              manager->FillNtupleDColumn(3,0,det_energy);
              manager->FillNtupleIColumn(3,1,copyNo);
              manager->AddNtupleRow(3);
          }
          else if (theStatus == NotAtBoundary)
          {
              procCount = "NotAtBoundary";
          }
          else if (theStatus == SameMaterial)
          {
              procCount = "SameMaterial";
          }
          else if (theStatus == StepTooSmall)
          {
              procCount = "SteptooSmall";
          }
          else if (theStatus == NoRINDEX)
          {
              procCount = "NoRINDEX";
          }
          else
          {
              G4cout << "theStatus: " << theStatus
                     << " was none of the above." << G4endl;
              procCount = "noStatus";
            }
      } // for if opProc
      // fill process histogram
      if(output){
      manager->FillNtupleSColumn(2,0,procCount);
      manager->AddNtupleRow(2);
    }

  }// for for statement
  return true;
}
