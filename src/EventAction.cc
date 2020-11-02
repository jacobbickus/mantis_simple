#include "EventAction.hh"


EventAction::EventAction(G4ParticleGun* particle_gun, SteppingAction* stepA)
: drawNumPhotonsFlag(0), eventM(NULL)
{
    eventM = new EventMessenger(this);
    particle_gun_local = particle_gun;
    stepA_local = stepA;
}

EventAction::~EventAction()
{delete eventM;}


void EventAction::BeginOfEventAction(const G4Event* anEvent)
{

    anEvent->GetEventID();
    ResetEverything(); //reset all the class variables

}


void EventAction::EndOfEventAction(const G4Event* anEvent)
{

    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    if(drawNumPhotonsFlag)
    {
      G4double E_beam=particle_gun_local->GetParticleEnergy();
      manager->FillNtupleDColumn(4,0,E_beam);
      manager->FillNtupleDColumn(4,1,stepA_local->Ev.size());
      manager->AddNtupleRow(4);
    }


    std::cout << "\r\tEvent and tracks:\t " << anEvent->GetEventID() << "\t" << stepA_local->Ev.size() <<std::flush;
}


void EventAction::ResetEverything(){

    stepA_local->Ev.clear();


}
