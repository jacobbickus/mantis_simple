#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;
class G4UserLimits;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
//private:
//    G4VPhysicalVolume *World;

public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    //G4VPhysicalVolume* GetDet() const {return physWater;}

    virtual G4VPhysicalVolume* Construct();
        //return physWorld;

    std::vector<G4ThreeVector> PmtPositions;
    void SetWaterX(G4double val)
    {
      water_size_x = val;
      water_size_x = water_size_x*cm;
    }
    void SetWaterY(G4double val)
    {
      water_size_y = val;
      water_size_y = water_size_y*cm;
    }
    void SetWaterZ(G4double val)
    {
      water_size_z = val;
      water_size_z = water_size_z*cm;
    }
    void SetPC_radius(G4double val)
    {PMT_rmax = val;
      PMT_rmax = PMT_rmax*cm;
    };


private:
  G4LogicalVolume *logicalLinac;
  G4LogicalVolume *logicBremTarget;
  G4LogicalVolume *logicaltarget;
    G4LogicalVolume* logicPC;
    G4UserLimits* fStepLimit;
    G4LogicalVolume* logicPMT;
    G4VPhysicalVolume* physPC;
    G4VPhysicalVolume* physPMT;
    G4VPhysicalVolume* physWater;
    G4double water_size_x, water_size_y, water_size_z,PMT_rmax;
    DetectorMessenger* detectorM;
};


#endif
