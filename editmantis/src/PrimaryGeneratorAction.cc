#include "PrimaryGeneratorAction.hh"
#include "PrimaryGenActionMessenger.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(),
chosen_energy(-1), resDataFlag(0), genM(NULL),fParticleGun(0)
{

        genM = new PrimaryGenActionMessenger(this);
        G4int n_particle = 1;
        fParticleGun = new G4ParticleGun(n_particle);

        // Default Kinematics
        fParticleGun->SetParticleDefinition(G4Gamma::Definition());

        fParticleGun->SetParticleTime(0.0*ns);
        // Determine particle Beam Energy
        ReadInputSpectrumFile("input_spectrum.txt");


#if defined (G4ANALYSIS_USE_ROOT)

        // file contains the normalized brems distribution p(E), sampling distribution s(E),
        // and binary 0/1 for off/on resonance useful in weighting
        TFile *fin = TFile::Open("brems_distributions.root");
        hBrems  = (TH1D*) fin->Get("hBrems");
        hSample = (TH1D*) fin->Get("hSample");
        hBinary = (TH1D*) fin->Get("hBinary");
        if (hBrems && hSample && hBinary)
        {
          G4cout << "Imported brems and sampling distributions from " << fin->GetName() << G4endl << G4endl;
        }
        else
        {
          std::cerr << "Error reading from file " << fin->GetName() << std::endl;
          exit(1);
        }

#endif
        // Set Beam Position
        beam_offset_x = 0*cm;
        beam_offset_y = 0*cm;
        z0 = 0*cm;
        beam_size = 1;
        source_width=0; //by default the width along Z is zero

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
        delete fParticleGun;
        delete genM;

}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

        // Set Particle Energy (Must be in generate primaries)
#if defined (G4ANALYSIS_USE_ROOT)

        if(chosen_energy == -2)
        {
          energy = hSample->GetRandom()*MeV; // sample the resonances specified by hSample
        }
#else
        if(chosen_energy == -2)
        {
          std::cerr << "ERROR: G4ANALYSIS_USE_ROOT not defined in pre-compiler" << std::endl;
          std::cerr << "SYSTEM EXITING" << std::endl;
          exit(100);
        }

#endif
        if(chosen_energy == -1)
        {
          random=G4UniformRand()*N[N.size()-1];

          // first order interpolation
          for(unsigned int i=0; i<N.size(); i++)
                  if(N[i]>random && i>0)
                  {
                      G4float f = (random - N[i-1]) / (N[i] - N[i-1]);
                      energy= f*e[i] + (1-f)*e[i-1];
                      break;
                  }
        }
        else if(chosen_energy != -2 && chosen_energy != -1)
        {
          energy = chosen_energy;
        }

        fParticleGun->SetParticleEnergy(energy*MeV);


        const float pi=acos(-1);

        // Set beam position
        r = beam_size*acos(G4UniformRand())/pi*2.;
        ph = 360.*G4UniformRand()*CLHEP::deg;
        x_r = r*cos(ph);
        y_r = r*sin(ph);
        z_r = source_width*(G4UniformRand()-0.5);
        fParticleGun->SetParticlePosition(G4ThreeVector(x_r+beam_offset_x,y_r+beam_offset_y,z_r+z0));


        // Set beam momentum
        theta = 0.*CLHEP::deg;
        phi = 0.*CLHEP::deg;
        G4ThreeVector vDir;
        vDir = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        fParticleGun->SetParticleMomentumDirection(vDir);

        fParticleGun->GeneratePrimaryVertex(anEvent);

        // if histogram exists find whether it was a resonance sample or not
#if defined (G4ANALYSIS_USE_ROOT)
        G4bool onRes;
        if(hBinary)
        {
          onRes = hBinary->GetBinContent(hBinary->GetXaxis()->FindBin(energy));
          G4int theRes;
          if(onRes)
          {
            theRes = 1;
          }
          else
          {
            theRes = 0;
          }
          if(resDataFlag)
          {
            G4AnalysisManager* manager = G4AnalysisManager::Instance();
            manager->FillNtupleIColumn(7,0,theRes);
            manager->AddNtupleRow(7);
          }
        }

#endif

}


void PrimaryGeneratorAction::ReadInputSpectrumFile(std::string filename){

        std::ifstream f(filename);
        if(f.is_open()) { //check that the file is open
                G4cout<<"oooooooooooooo Reading input file for beam energies oooooooooooo"<<G4endl;
                N.push_back(0);
                while(!f.eof()) {
                        float a,b;
                        f>>a>>b;
                        e.push_back(a);
                        dNde.push_back(b);
                        if(dNde.size()>1) N.push_back(N.at(N.size()-1)+dNde.at(dNde.size()-1)*(e.at(e.size()-1)-e.at(e.size()-2)));
                }
        }
        else {
                std::cerr<<"Input file expected, however open failed, exiting."<<std::endl;
                exit(8);
        }

        f.close();
}
