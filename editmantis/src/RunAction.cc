#include "G4Timer.hh"
#include "RunAction.hh"
#include "G4Run.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "EventAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"


RunAction::RunAction(HistoManager* histoAnalysis)
        : G4UserRunAction(),fTimer(nullptr), fHistoManager(histoAnalysis)
{

        fTimer = new G4Timer;
}

RunAction::~RunAction()
{
        delete fTimer;
}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{


        fHistoManager->Book();
        G4cout << "Writing Results to analysis file" << G4endl;
        aRun->GetNumberOfEvent(); // not necessary but gets rid of warning
        fTotalSurface = 0;
        fCerenkovCount = 0;
        fScintCount = 0;
        fRayleighCount = 0;
        fOpAbsorption = 0;
        fCerenkovEnergy = 0;
        fScintEnergy = 0;
        fNRF = 0;
        fStatusKilled = 0;

        fTimer->Start();
        setStartTime();
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
        fTimer->Stop();
        // save histograms
        G4cout << G4endl << "RunAction::EndOfRunAction: Histogram Saved." << G4endl;

        //G4cout << "Made to end of run" << G4endl;
        G4int TotNbofEvents = aRun->GetNumberOfEvent();
        std::ios::fmtflags mode = G4cout.flags();
        G4int prec = G4cout.precision(2);
        //G4double hits = G4double(fHitCount)/TotNbofEvents;
        //G4double hitsAbove = G4double(fPMTsAboveThreshold)/TotNbofEvents;
        G4cout << G4endl << "Run Summary" << G4endl;
        G4cout <<   "---------------------------------" << G4endl;
        G4cout << "Total Number of Events: " << TotNbofEvents << G4endl;
        G4cout << "Total number of Surface Events: " << fTotalSurface << G4endl;
        G4cout << "Total number of NRF Photons: " << fNRF << G4endl;
        G4cout << "Total number of Tracks Cut Based on Position: " << fStatusKilled << G4endl;
        G4cout <<
                "Average total energy of Cerenkov photons created in the Water per event: " << (fCerenkovEnergy/eV)/TotNbofEvents << " eV." << G4endl;
        G4cout << "Average number of Cerenkov photons created in the Water per event: " << fCerenkovCount/TotNbofEvents << G4endl;

        if (fCerenkovCount > 0) {
                G4cout << " Average Cherenkov Photon energy emitted in Water: " << (fCerenkovEnergy/eV)/fCerenkovCount << " eV." << G4endl;
        }

        G4cout <<
                "Average total energy of scintillation photons created per event in the Water: " << (fScintEnergy/eV)/TotNbofEvents << " eV." << G4endl;
        G4cout <<
                "Average number of scintillation photons created per event in the Water: " << fScintCount/TotNbofEvents << G4endl;

        if (fScintCount > 0) {
                G4cout << " Average Scintillation Photon energy emitted in the Water: " << (fScintEnergy/eV)/fScintCount << " eV." << G4endl;
        }

        G4cout << "Average number of OpRayleigh per event:   " << fRayleighCount/TotNbofEvents << G4endl;
        G4cout << "Average number of OpAbsorption per event: " << fOpAbsorption/TotNbofEvents << G4endl;

        //G4double bdry = G4double(fBoundaryAbsorptionCount)/TotNbofEvents;
        //G4cout << "Number of photons absorbed at boundary per event: " << bdry << G4endl;

        G4cout <<   "---------------------------------" << G4endl;

        G4cout.setf(mode, std::ios::floatfield);

        G4cout.precision(prec);
        // Also print to screen and not just to log file
        std::cout << std::endl << "Run Summary" << std::endl;
        std::cout <<   "---------------------------------" << std::endl;
        std::cout << "Total Number of Events: " << TotNbofEvents << std::endl;
        std::cout << "Total number of Surface Events: " << fTotalSurface << std::endl;
        std::cout << "Total number of NRF Photons: " << fNRF << std::endl;
        std::cout << "Total number of Tracks Cut Based on Position: " << fStatusKilled << std::endl;
        std::cout <<
                "Average total energy of Cerenkov photons created in the Water per event: " << (fCerenkovEnergy/eV)/TotNbofEvents << " eV." << std::endl;
        std::cout << "Average number of Cerenkov photons created in the Water per event: " << fCerenkovCount/TotNbofEvents << std::endl;

        if (fCerenkovCount > 0) {
                std::cout << " Average Cherenkov Photon energy emitted in Water: " << (fCerenkovEnergy/eV)/fCerenkovCount << " eV." << std::endl;
        }

        std::cout <<
                "Average total energy of scintillation photons created per event in the Water: " << (fScintEnergy/eV)/TotNbofEvents << " eV." << std::endl;
        std::cout <<
                "Average number of scintillation photons created per event in the Water: " << fScintCount/TotNbofEvents << std::endl;

        if (fScintCount > 0) {
                std::cout << " Average Scintillation Photon energy emitted in the Water: " << (fScintEnergy/eV)/fScintCount << " eV." << std::endl;
        }

        std::cout << "Average number of OpRayleigh per event:   " << fRayleighCount/TotNbofEvents << std::endl;
        std::cout << "Average number of OpAbsorption per event: " << fOpAbsorption/TotNbofEvents << std::endl;

        //G4double bdry = G4double(fBoundaryAbsorptionCount)/TotNbofEvents;
        //std::cout  << "Number of photons absorbed at boundary per event: " << bdry << std::endl;

        std::cout <<   "---------------------------------" << std::endl;

        fHistoManager->finish();

}
