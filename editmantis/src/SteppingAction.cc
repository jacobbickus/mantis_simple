#include "SteppingAction.hh"
#include "G4Cerenkov.hh"


SteppingAction::SteppingAction(const DetectorConstruction* det, G4ParticleGun* particle_gun, RunAction* localrun)
: G4UserSteppingAction(), drawIntObjDataFlag(0), drawWaterFlag(0), stepM(NULL)
{
    stepM = new StepMessenger(this);
    local_det = det;
    particle_gun_local = particle_gun;
    run = localrun;
}

SteppingAction::~SteppingAction()
{ delete stepM; }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    G4int isNRF = 0;

    G4String particleName = aStep->GetTrack()->GetDynamicParticle()
    ->GetParticleDefinition()->GetParticleName();
    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    // Run Logical Checks
    if(aStep->GetPostStepPoint()==NULL){
        return; // at the end of the world
    }
    else if(aStep->GetPostStepPoint()->GetPhysicalVolume()==NULL){
        return;
    }

    G4Track* theTrack = aStep->GetTrack();
    // kill photons past IntObj
    /*G4double EndIntObj = local_det->getEndIntObj();

    if(theTrack->GetPosition().z()/(cm) > EndIntObj/(cm))
    {
      // kill photons that go beyond the interrogation object
      theTrack->SetTrackStatus(fStopAndKill);
      run->AddStatusKilled();
    }
    else if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName().compare(0, 10, "Collimator") == 0)
    {
      // kill photons in collimator
      theTrack->SetTrackStatus(fStopAndKill);
      run->AddStatusKilled();
    }
    else if(theTrack->GetPosition().z()/(cm) < -1.*cm)
    {
      // Kill photons that go in behind beam origin
      theTrack->SetTrackStatus(fStopAndKill);
      run->AddStatusKilled();
    }
    */
    if (particleName == "opticalphoton") {
         const G4VProcess* pds = aStep->GetPostStepPoint()->
                                                      GetProcessDefinedStep();
         if (pds->GetProcessName() == "OpAbsorption") {
             run->AddOpAbsorption();

         }
         else if (pds->GetProcessName() == "OpRayleigh") {
             run->AddRayleigh();
         }
     }

     if(theTrack->GetCreatorProcess() !=0)
     {
       G4String CPName = theTrack->GetCreatorProcess()->GetProcessName();
       if(CPName == "NRF")
       {
         run->AddNRF();
         isNRF = 1;
       }
     }

     // Testing NRF Analysis
     // inside Interogation Object for first time
     if(drawIntObjDataFlag)
     {
       if(aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName().compare(0, 14 ,"IntObjPhysical") == 0
          && aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName().compare(0, 14, "IntObjPhysical") != 0)
       {
         G4double energy_IntObj = theTrack->GetKineticEnergy()/(MeV);
         manager->FillNtupleDColumn(5,0,energy_IntObj);
         manager->FillNtupleIColumn(5,1,isNRF);
         manager->AddNtupleRow(5);
       }
     }

    // Water Analysis
    // first time in detector determine incident water energies
    if(drawWaterIncDataFlag)
    {
      if(aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName().compare(0, 5 ,"Water") == 0
         && aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName().compare(0, 5, "Water") != 0)
         {
           G4double energy_inc_water = theTrack->GetKineticEnergy()/(MeV);
           manager->FillNtupleDColumn(6,0, energy_inc_water);
           manager->FillNtupleIColumn(6,1, isNRF);
           manager->AddNtupleRow(6);
         }
    }

    // Here I am inside the water
    if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName().compare(0,5,"Water")==0){

        // only care about secondaries that occur in water volume
        const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();

        if(secondaries->size()>0){
            for(unsigned int i=0; i<secondaries->size(); ++i){
                if(secondaries->at(i)->GetParentID()>0){
                    if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
                        if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Scintillation"){
                            G4double en = secondaries->at(i)->GetKineticEnergy();
                            run->AddScintillationEnergy(en);
                            run->AddScintillation();
                        }
                        if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
                            // for total run
                            G4double en = secondaries->at(i)->GetKineticEnergy();
                            run->AddCerenkovEnergy(en);
                            run->AddCerenkov();
                            //for step

                        }

                        // only care about cerernkov
                        // root file
                        X= secondaries->at(i)->GetPosition(); //take the position from the post step position
                        p = secondaries->at(i)->GetMomentum();
                        if(drawWaterFlag)
                        {
                          manager->FillNtupleDColumn(0,0, secondaries->at(i)->GetKineticEnergy()/(MeV));
                          manager->FillNtupleIColumn(0,1, isNRF);
                          manager->FillNtupleDColumn(0,2, X.x()/(cm));
                          manager->FillNtupleDColumn(0,3, X.y()/(cm));
                          manager->FillNtupleDColumn(0,4, X.z()/(cm));
                          manager->FillNtupleDColumn(0,5, asin(sqrt(pow(p.x(),2)+pow(p.y(),2))/p.mag()));
                          manager->FillNtupleDColumn(0,6, secondaries->at(i)->GetGlobalTime());
                          manager->AddNtupleRow(0);

                        }

                        Ev.push_back( secondaries->at(i)->GetKineticEnergy()/(MeV));// records for histogram


                    }
                }
            }
        } // end of optical photons if statement

    } // end of if loop while inside water

} // end of user steepping action function


bool SteppingAction::TrackMustDie(const G4Step *aStep){

    //Energy cuts
    if(aStep->GetTrack()->GetDefinition()->GetPDGEncoding()==EventGeneratorParticle &&
       aStep->GetTrack()->GetKineticEnergy()/MeV < LowEnergyCutoff){ //for the main track, kill it off below some limit
      return true;
    }

    return false;
}
