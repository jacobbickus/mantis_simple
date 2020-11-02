#include "PrimaryGeneratorAction.hh"
#include "PrimaryGenActionMessenger.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(),
chosen_energy(-1), genM(NULL),fParticleGun(0)
{

        genM = new PrimaryGenActionMessenger(this);
        G4int n_particle = 1;
        fParticleGun = new G4ParticleGun(n_particle);

        // Default Kinematics
        fParticleGun->SetParticleDefinition(G4Electron::Definition());

        fParticleGun->SetParticleTime(0.0*ns);
   
        // Set Beam Position
        beam_offset_x = 0*cm;
        beam_offset_y = 0*cm;
        z0 = 0*cm;
        beam_size = 1.3*mm;
        source_width=0; //by default the width along Z is zero

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
        delete fParticleGun;
        delete genM;

}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   energy = chosen_energy;

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

}
