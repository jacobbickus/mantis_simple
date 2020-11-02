#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "PCSensitiveDetector.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SubtractionSolid.hh"

#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4OpRayleigh.hh"

#include "G4Box.hh" //for box
//#include "G4Cons.hh" //for cone
//#include "G4Orb.hh" //for orb
//#include "G4Sphere.hh" //for sphere
#include "G4Tubs.hh" // for tube/cylinder

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4PhysicalConstants.hh"
// for color attributes
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "DetectorMessenger.hh"
#include "G4UserLimits.hh"


DetectorConstruction::DetectorConstruction()
        : G4VUserDetectorConstruction(),water_size_x(100*cm),water_size_y(100*cm), water_size_z(50*cm),PMT_rmax(50*cm), detectorM(NULL)
{
  detectorM = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
  delete detectorM;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  // Set World Dimensions
  G4double world_size_x = 300*cm;
  G4double world_size_y = 300*cm;
  G4double world_size_z = 300*cm;

  // Make Solid World

  G4Box* solidWorld = new G4Box("World", world_size_x, world_size_y, world_size_z);
  G4Element *elN = new G4Element("Nitrogen", "N2", 7, 14.01*g/mole);
  G4Element *elO = new G4Element("Oxygen", "O2", 8, 16.0*g/mole);
  G4Element* elC = new G4Element("Carbon", "C", 6, 12.0*g/mole);
  G4Element* elH = new G4Element("Hydrogen","H",1, 1.0*g/mole);
  //G4Element* elF = new G4Element("Fluorine", "F", 9, 0.001696*g/mole);
  G4Material *myAir = new G4Material("Air", 1.290*mg/cm3, 2);
  myAir->AddElement(elN, 0.7);
  myAir->AddElement(elO, 0.3);
  G4Material *myVacuum = new G4Material("Vacuum", 1.e-5*g/cm3, 1, kStateGas, 273.15, 2.e-2*bar);
  myVacuum->AddMaterial(myAir,1);
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material *tungsten = nist->FindOrBuildMaterial("G4_W");

  // Make Logical Volume
  G4LogicalVolume* logicWorld =
          new G4LogicalVolume(solidWorld, //its solid
                              air, //its material
                              "World"); //its name

  // Make Physical volume ** NEVER CHANGE THIS **
  G4VPhysicalVolume* physWorld =
          new G4PVPlacement(0,         //no rotation
                            G4ThreeVector(), //at (0,0,0)
                            logicWorld, //its logical volume
                            "World",   //its name
                            0,         //its mother  volume
                            false,     //no boolean operation
                            0,         //copy number
                            false);    //overlaps checking

  // Make linac
  G4bool testBrem = true;
  if(testBrem)
  {
    G4Tubs *solidLinac = new G4Tubs("Linac",0, 2*cm, 3*cm, 0*deg, 360*deg);
    logicalLinac = new G4LogicalVolume(solidLinac, tungsten, "Linac");
    new G4PVPlacement(0, G4ThreeVector(0,0, 3*cm), logicalLinac, "Linac", logicWorld, false, 0, checkOverlaps);
    G4Tubs *solidVacuum = new G4Tubs("Vacuum", 0, 10*mm, 3*cm, 0*deg, 360*deg);
    G4LogicalVolume *logicalVacuum = new G4LogicalVolume(solidVacuum, myVacuum, "Vacuum");
    new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalVacuum, "Vacuum", logicalLinac, false,0,checkOverlaps);
    // Make Brem target
    G4Box *solidBremTarget = new G4Box("Brem", 2*mm, 2*mm, 0.1*mm);
    logicBremTarget = new G4LogicalVolume(solidBremTarget, tungsten, "Brem");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -2*cm),logicBremTarget,"Brem", logicalVacuum, false, 0, checkOverlaps);

  // Tub of water
  std::cout << "The Water Tank X was set to: " << water_size_x/(cm)<< " cm" << std::endl;
  std::cout << "The Water Tank Y was set to: " << water_size_y/(cm)<< " cm" << std::endl;
  std::cout << "The Water Tank Z was set to: " << water_size_z/(cm) << " cm" << std::endl;

  // Make test target
  G4Box *solidTarget = new G4Box("target",5*cm, 5*cm, 1*cm);
  logicaltarget = new G4LogicalVolume(solidTarget, tungsten,"Target");
  new G4PVPlacement(0,G4ThreeVector(0, 0, 25*cm), logicaltarget, "Target", logicWorld, false, 0, checkOverlaps);
}

// Make Physical volume
G4double water_x_pos = 0*cm;
G4double water_y_pos =0*cm;
G4double water_z_pos =100*cm;

// Make Solid Water casing
G4Box* solidCasing = new G4Box("Encasing", water_size_x, water_size_y, water_size_z);
G4Material *plexiglass = new G4Material("Plexiglass", 1.19*g/cm3,3);
plexiglass->AddElement(elC, 5);
plexiglass->AddElement(elH, 8);
plexiglass->AddElement(elO, 2);
G4LogicalVolume* logicCasing = new G4LogicalVolume(solidCasing, plexiglass, "Encasing");
new G4PVPlacement(0,
  G4ThreeVector(water_x_pos,water_y_pos,water_z_pos), logicCasing,
  "Encasing", logicWorld, false, 0, checkOverlaps);

G4double plexiThickness = 0.18*mm; //0.18*mm;
// Make Teflon tape wrap
G4double tapeThick = 0.01*cm;
G4VSolid* solidTape = new G4Box("Tape", water_size_x-plexiThickness, water_size_y-plexiThickness, water_size_z-plexiThickness);
//G4VSolid* subtractFrontTape = new G4Box("SubtractTape", water_size_x+1*cm, water_size_y+1*cm, tapeThick);
//G4SubtractionSolid* solidTape = new G4SubtractionSolid("Tape", solidWholeTape, subtractFrontTape, 0, G4ThreeVector(0,0, -water_size_z)); // name, first solid, subtraction solid, rotation, translation
G4Material *teflonTape = nist->FindOrBuildMaterial("G4_TEFLON");
G4LogicalVolume* logicTape = new G4LogicalVolume(solidTape, teflonTape, "Tape");
G4PVPlacement* physTape = new G4PVPlacement(0,G4ThreeVector(0, 0, 0), logicTape, "Tape", logicCasing, false, 0, checkOverlaps);

// Make Solid Water
G4Box* solidWater = new G4Box("Water", water_size_x - plexiThickness - tapeThick,
water_size_y- plexiThickness - tapeThick, water_size_z - plexiThickness - tapeThick);
G4Material* Water = nist->FindOrBuildMaterial("G4_WATER");

// Make Logical Volume
G4LogicalVolume* logicWater =
        new G4LogicalVolume(solidWater, //its solid
                            Water, //its material
                            "Water"); //its name

physWater = new G4PVPlacement(0,                 //no rotation
                              G4ThreeVector(0,0,0),
                              logicWater, //its logical volume
                              "Water", //its name
                              logicTape, //its mother logical volume
                              false, //no boolean operation
                              0,     //copy number
                              checkOverlaps); //overlaps checking

  // PMT Cylinder(Tube) face

  G4double PMT_rmin = 0*cm;
  //PMT_rmax = 100*cm; // made it 10 times larger for now
  std::cout << "The PC Radius was set to " << PMT_rmax/(cm) << " cm" << std::endl;
  G4double PMT_z = 6*2.54*cm;
  G4double PMT_start_phi = 0.*deg;
  G4double PMT_phi = 360.*deg;// may not work with pi

  // Make Solid PMT

  G4Tubs* solidPMT = new G4Tubs("PMT", PMT_rmin, PMT_rmax, PMT_z, PMT_start_phi, PMT_phi);
  G4Material* PMT_mat = nist->FindOrBuildMaterial("G4_Pyrex_Glass"); // technically PMT glass is a special borosilicate glass calle k-free glass but pyrex is close enough as a borosilicate glass

  // Make Logical PMT

  logicPMT = new G4LogicalVolume(solidPMT, PMT_mat, "PMT");
  // Set as Sensitive detector



  // Make Physical PMT

  G4double PMT_x_pos = 0*cm;
  G4double PMT_y_pos = 0*cm;
  G4double PMT_z_pos = (water_z_pos + water_size_z + PMT_z);
  PmtPositions.push_back(G4ThreeVector(PMT_x_pos,PMT_y_pos,PMT_z_pos));

  physPMT = new G4PVPlacement(0,
                              G4ThreeVector(0,0,water_size_x/2 - PMT_z - plexiThickness - tapeThick),
                              logicPMT,
                              "PMT",
                              logicWater,
                              false,
                              0,
                              checkOverlaps);

  // Make Solid PhotoCathode and set tracking limits
  G4double PC_z = 20.*nm; // here if set below 1cm results vary greatly
  G4Tubs* solidPhotoCathode = new G4Tubs("PC", PMT_rmin, PMT_rmax, PC_z, PMT_start_phi, PMT_phi);
  G4Material* PC_mat = nist->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");

  //Make Logical Photocathode volume
  logicPC = new G4LogicalVolume(solidPhotoCathode, PC_mat, "PC");
  G4double PMT_window_thickness = 3*mm;
  // Make physical volume
  physPC = new G4PVPlacement(0,
                             G4ThreeVector(0,0,-PMT_z + PMT_window_thickness), // now position is relation to mother volume(PMT)
                             logicPC,
                             "PC",
                             logicPMT, // daughter of PMT logical
                             false,
                             0,
                             checkOverlaps);

  // These didn't actually help just slowed code down
  //G4double maxPCStep = 1*nm;
  //fStepLimit = new G4UserLimits(maxPCStep);
  //logicPC->SetUserLimits(fStepLimit);
  //logicPMT->SetUserLimits(fStepLimit);
  //logicWater->SetUserLimits(fStepLimit);

  // Visualization Attributes

  //logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible()); // sets world invisible in visualization
  G4VisAttributes * yellow= new G4VisAttributes( G4Colour(255/255.,255/255.,51/255. ));
  G4VisAttributes * red= new G4VisAttributes( G4Colour(255/255., 0/255., 0/255. ));
  G4VisAttributes * blue= new G4VisAttributes( G4Colour(0/255., 0/255.,  255/255. ));
  G4VisAttributes * grayc= new G4VisAttributes( G4Colour(128/255., 128/255.,  128/255. ));
  G4VisAttributes * lightGray= new G4VisAttributes( G4Colour(178/255., 178/255.,  178/255. ));
  G4VisAttributes * green= new G4VisAttributes( G4Colour(0/255., 255/255.,  0/255. ));

  yellow->SetVisibility(true);
  yellow->SetForceSolid(true);
  red->SetVisibility(true);
  red->SetForceWireframe(true);
  blue->SetVisibility(true);
  green->SetVisibility(true);
  green->SetForceWireframe(true);
  grayc->SetVisibility(true);
  grayc->SetForceSolid(true);
  lightGray->SetVisibility(true);
  lightGray->SetForceSolid(true);

  // Set Visual colors
  if(testBrem)
  {
    logicaltarget->SetVisAttributes(blue);
    logicBremTarget->SetVisAttributes(yellow);
  }
  logicCasing->SetVisAttributes(red);
  logicTape->SetVisAttributes(yellow);
  logicWater->SetVisAttributes(blue);
  logicPMT->SetVisAttributes(green);
  logicPC->SetVisAttributes(red);


  //
  // ------------ Generate & Add Material Properties Table ------------
  //
  G4double photonEnergy[] =
  { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
    2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
    2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
    2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
    2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
    3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
    3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
    3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  G4double photonEnergy2[] =
  {1.72*eV, 1.746197183*eV, 1.771142857*eV, 1.796811594*eV,1.823235294*eV,
    1.850447761*eV, 1.878484848*eV, 1.907384615*eV, 1.9371875*eV,
    1.967936508*eV, 1.999677419*eV, 2.032459016*eV, 2.066333333*eV,
    2.101355932*eV,2.137586207*eV, 2.175087719*eV, 2.213928571*eV, 2.254181818*eV,
   2.295925926*eV, 2.339245283*eV, 2.384230769*eV, 2.430980392*eV, 2.4796*eV,
   2.530204082*eV, 2.582916667*eV, 2.63787234*eV, 2.695217391*eV,
    2.755111111*eV, 2.817727273*eV, 2.883255814*eV, 2.951904762*eV,
    3.023902439*eV, 3.0995*eV, 3.178974359*eV, 3.262631579*eV,
    3.542285714*eV, 4.132666667*eV, 4.508363636*eV, 4.9592*eV
  };


  //
  // Water need to add rayleigh scattering! -- added automatically if material name = Water
  //
  G4double refractiveIndex1[] =
  { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
    1.346,  1.3465, 1.347,  1.3475, 1.348,
    1.3485, 1.3492, 1.35,   1.3505, 1.351,
    1.3518, 1.3522, 1.3530, 1.3535, 1.354,
    1.3545, 1.355,  1.3555, 1.356,  1.3568,
    1.3572, 1.358,  1.3585, 1.359,  1.3595,
    1.36,   1.3608};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] =
  {0.812347685*m, 1.209189843*m, 1.602564103*m, 1.937984496*m, 2.150537634*m,
   2.277904328*m, 2.43902439*m, 2.941176471*m, 3.217503218*m, 3.429355281*m,
   3.629764065*m, 3.78214826*m, 4.496402878*m, 7.4019245*m, 11.16071429*m, 14.38848921*m,
   16.15508885*m, 17.69911504*m, 21.09704641*m, 23.04147465*m, 24.44987775*m,
   30.76923077*m, 49.01960784*m,66.66666667*m,78.74015748*m,94.33962264*m,
   102.145046*m, 108.4598698*m, 157.480315*m, 202.020202*m, 220.2643172*m,
   211.4164905*m, 150.8295626*m, 117.5088132*m, 87.95074758*m,
   60.60606061*m,18.11594203*m,11.79245283*m,9.652509653*m};

  assert(sizeof(absorption) == sizeof(photonEnergy2));

  G4double scintilFast[] =
  { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4double scintilSlow[] =
  { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
    7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
    3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
    4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
    7.00, 6.00, 5.00, 4.00 };

  assert(sizeof(scintilSlow) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1,nEntries)->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH",photonEnergy2, absorption, nEntries)->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast, nEntries)->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow, nEntries)->SetSpline(true);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

  G4double energy_water[] = {
          1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
          1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
          1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
          1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
          1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
          2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
          2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
          2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
          2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
          2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
          3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
          3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
          3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
          4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
          5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
  };

  const G4int numentries_water = sizeof(energy_water)/sizeof(G4double);

  G4double energy_ray_water[] = {
      1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
      1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
      1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
      1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
      1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
      2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
      2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
      2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV
     }; // 32 entries



  const G4int numentries_water_ray = sizeof(energy_ray_water)/sizeof(G4double);

  G4double mie_water[] = {
          167024.4*m, 158726.7*m, 150742  *m,
          143062.5*m, 135680.2*m, 128587.4*m,
          121776.3*m, 115239.5*m, 108969.5*m,
          102958.8*m, 97200.35*m, 91686.86*m,
          86411.33*m, 81366.79*m, 76546.42*m,
          71943.46*m, 67551.29*m, 63363.36*m,
          59373.25*m, 55574.61*m, 51961.24*m,
          48527.00*m, 45265.87*m, 42171.94*m,
          39239.39*m, 36462.50*m, 33835.68*m,
          31353.41*m, 29010.30*m, 26801.03*m,
          24720.42*m, 22763.36*m, 20924.88*m,
          19200.07*m, 17584.16*m, 16072.45*m,
          14660.38*m, 13343.46*m, 12117.33*m,
          10977.70*m, 9920.416*m, 8941.407*m,
          8036.711*m, 7202.470*m, 6434.927*m,
          5730.429*m, 5085.425*m, 4496.467*m,
          3960.210*m, 3473.413*m, 3032.937*m,
          2635.746*m, 2278.907*m, 1959.588*m,
          1675.064*m, 1422.710*m, 1200.004*m,
          1004.528*m, 833.9666*m, 686.1063*m
  };

  assert(sizeof(mie_water) == sizeof(energy_water));

  // gforward, gbackward, forward backward ratio
  G4double mie_water_const[3]={0.99,0.99,0.8};

  myMPT1->AddProperty("MIEHG",energy_water,mie_water,numentries_water)->SetSpline(true);
  myMPT1->AddConstProperty("MIEHG_FORWARD",mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD",mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",mie_water_const[2]);

  G4double rayleigh_water[] = {
      2335.5*m, 2210.5*m, 2090.8*m, 1976.4*m, 1866.8*m, 1762.1*m,
      1662.1*m, 1566.6*m, 1475.4*m, 1388.5*m, 1305.6*m, 1224.6*m,
      1146.8*m, 1075.6*m, 1007.8*m, 941.23*m, 880.95*m, 821.12*m,
      766.36*m, 714.50*m, 665.40*m, 618.97*m, 575.11*m, 533.69*m,
      493.46*m, 457.10*m, 421.51*m, 389.06*m, 358.57*m, 329.98*m,
      303.17*m, 277.43*m

     }; // 32 entries


  // Note to Self: Never have spline set to true if rayleigh water goes to zero. Creates seg fault

  assert(sizeof(rayleigh_water) == sizeof(energy_ray_water));

  myMPT1->AddProperty("RAYLEIGH",energy_ray_water,rayleigh_water,numentries_water_ray)->SetSpline(true);

  //G4cout << "Water G4MaterialPropertiesTable" << G4endl;
  //myMPT1->DumpTable();

  Water->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator

  Water->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // APPLY OPTICAL QUALITIES

  /* ****************************************************************************************************************************************************** */

G4MaterialPropertiesTable* casingMPT = new G4MaterialPropertiesTable();
G4MaterialPropertiesTable* casingOPMPT = new G4MaterialPropertiesTable();
// add all of the properties here for plexiglass
G4double ephotonPlexi[] = {1.1808067*eV,1.2336786*eV,1.298269*eV,1.3699967*eV,1.4501135*eV,1.5401826*eV,1.6421815*eV,1.7586482*eV,1.8928962*eV,2.049334*eV,
                           2.2339586*eV,2.4551426*eV, 2.7249385*eV,2.850223*eV,3.0*eV,4.0*eV};
const G4int numPlexi = sizeof(ephotonPlexi)/sizeof(G4double);
G4double plexiReflectivity[] = {0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02}; // transparent
G4double plexiTrans[] = {.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98,.98};
G4double plexiIOF[] = {1.4813,1.482,1.4827,1.4832,1.4837,1.4842,1.485,1.4862,1.488,1.4904,1.4933,1.4962,1.4996,1.5027,1.5027,1.5027};
G4double plexiEff[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
G4double plexiAbs[] = {1000.0*m,1000.0*m,1000.0*m,1000.0*m,1000.0*m,1000.0*m,1000.0*m,1000.0*m,30.0*m,30.0*m,11.4*m,11.4*m,7.3*m,7.3*m, 5.0*m, 1.0*mm}; // pretty much transparent until UV energy
assert(sizeof(plexiIOF) == sizeof(ephotonPlexi));
assert(sizeof(plexiReflectivity) == sizeof(ephotonPlexi));
assert(sizeof(plexiTrans) == sizeof(ephotonPlexi));
assert(sizeof(plexiEff) == sizeof(ephotonPlexi));
assert(sizeof(plexiAbs) == sizeof(ephotonPlexi));
G4OpticalSurface *casing_opsurf = new G4OpticalSurface("casingSurface", glisur, polished, dielectric_dielectric);
casingOPMPT->AddProperty("REFLECTIVITY", ephotonPlexi, plexiReflectivity, numPlexi);
casingOPMPT->AddProperty("TRANSMITTANCE", ephotonPlexi, plexiTrans, numPlexi);
casingMPT->AddProperty("RINDEX", ephotonPlexi, plexiIOF, numPlexi);
casingOPMPT->AddProperty("RINDEX", ephotonPlexi, plexiIOF, numPlexi);
casingOPMPT->AddProperty("EFFICIENCY", ephotonPlexi, plexiEff, numPlexi);
casingMPT->AddProperty("ABSLENGTH", ephotonPlexi, plexiAbs, numPlexi);
casingOPMPT->AddProperty("ABSLENGTH", ephotonPlexi, plexiAbs, numPlexi);

plexiglass->SetMaterialPropertiesTable(casingMPT);
casing_opsurf->SetMaterialPropertiesTable(casingOPMPT);
new G4LogicalSkinSurface("casing_surf", logicCasing, casing_opsurf);

// Add properties for Teflon Tape
G4MaterialPropertiesTable* tapeMPT = new G4MaterialPropertiesTable();
G4MaterialPropertiesTable* tapeOPMPT = new G4MaterialPropertiesTable();
G4double ephotonTape[] = {1.1427161*eV,1.2915073*eV,1.4586435*eV,1.6984206*eV,2.0160114*eV,2.5303*eV,3.262755*eV, 4.0*eV};
const G4int numTape = sizeof(ephotonTape)/sizeof(G4double);
G4double tapeReflectivity[] = {.98,.98,.98,.98,.98,.98,.98,.98};
G4double tapeTrans[] = {0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02};
//G4double backScatter[] = {1.,1.,1.,1.,1.,1.,1.,1.};
G4double tapeEff[] = {0.,0.,0.,0.,0.,0.,0.,0.};
G4double tapeIOF[] = {1.32,1.32,1.32,1.32,1.32,1.32,1.33,1.33};
G4double tapeAbs[] = {1000*m, 1000*m,1000*m,1000*m,1000*m,1000*m,1000*m,1000*m}; // setting absorption length arbitrarily large, don't want teflon to kill opt photons
assert(sizeof(tapeReflectivity) == sizeof(ephotonTape));
//assert(sizeof(backScatter) == sizeof(ephotonTape));
assert(sizeof(tapeEff) == sizeof(ephotonTape));
assert(sizeof(tapeIOF) == sizeof(ephotonTape));
assert(sizeof(tapeAbs) == sizeof(ephotonTape));
assert(sizeof(tapeTrans) == sizeof(ephotonTape));
G4OpticalSurface *tape_opsurf = new G4OpticalSurface("tapeSurface",glisur, polishedfrontpainted, dielectric_dielectric);
tapeOPMPT->AddProperty("REFLECTIVITY", ephotonTape, tapeReflectivity, numTape);
tapeOPMPT->AddProperty("TRANSMITTANCE", ephotonTape, tapeTrans, numTape);
//tapeOPMPT->AddProperty("BACKSCATTERCONSTANT", ephotonTape, backScatter, numTape);
//tapeMPT->AddProperty("REFLECTIVITY", ephotonTape, tapeReflectivity, numTape)->SetSpline(true);
//tapeOPMPT->AddProperty("RINDEX", ephotonTape, tapeIOF, numTape)->SetSpline(true);
tapeMPT->AddProperty("RINDEX", ephotonTape, tapeIOF, numTape);
tapeOPMPT->AddProperty("EFFICIENCY", ephotonTape, tapeEff, numTape);
tapeMPT->AddProperty("ABSLENGTH", ephotonTape, tapeAbs, numTape);

tape_opsurf->SetMaterialPropertiesTable(tapeOPMPT);
teflonTape->SetMaterialPropertiesTable(tapeMPT);
//new G4LogicalSkinSurface("tapeSurf", logicTape, tape_opsurf);
new G4LogicalBorderSurface("tape_surf", physWater, physTape, tape_opsurf);

  // PMT
  G4double ephotonPMT[] = {1.7711*eV, 1.9074*eV, 2.0663*eV, 2.2542*eV, 2.4796*eV, 2.7551*eV, 3.0995*eV, 3.5423*eV, 4.133*eV};
  const G4int num = sizeof(ephotonPMT)/sizeof(G4double);

  // PMT Surface Properties
  G4double reflectivity[] = {0.038, 0.0392,  0.040,  0.041,
    0.0415,  0.0428,  0.0433,  0.0440, 0.0451 };
  assert(sizeof(reflectivity) == sizeof(ephotonPMT));
  G4double transmittance[] = {0.962, 0.961, 0.96, 0.959, 0.959, 0.958, 0.957,0.956,0.955};
  G4double abslength[] = {1000*m,1000*m,1000*m,1000*m,1000*m,1000*m,1000*m,1000*m,1000*m};
  assert(sizeof(abslength) == sizeof(ephotonPMT));
  G4double efficiency[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  assert(sizeof(efficiency) == sizeof(ephotonPMT));
  G4double pmt_iof[] = {1.486, 1.4945, 1.5013, 1.5075, 1.5118, 1.5214, 1.5253, 1.5308, 1.539};
  assert(sizeof(pmt_iof) == sizeof(ephotonPMT));

  G4OpticalSurface* PMT_opsurf = new G4OpticalSurface("PMTSurface",unified,polished,dielectric_dielectric);
  G4MaterialPropertiesTable* PMTopt = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* PMTabs = new G4MaterialPropertiesTable();
  PMTopt->AddProperty("REFLECTIVITY", ephotonPMT, reflectivity, num);
  PMTopt->AddProperty("TRANSMITTANCE",ephotonPMT, transmittance, num);
  PMTopt->AddProperty("EFFICIENCY", ephotonPMT, efficiency, num);
  PMTopt->AddProperty("RINDEX", ephotonPMT, pmt_iof, num);
  PMTabs->AddProperty("RINDEX", ephotonPMT, pmt_iof, num);
  PMTabs->AddProperty("ABSLENGTH", ephotonPMT, abslength, num);
  PMT_opsurf->SetMaterialPropertiesTable(PMTopt);
  PMT_mat->SetMaterialPropertiesTable(PMTabs);

  new G4LogicalBorderSurface("pmt-water_surf", physWater, physPMT, PMT_opsurf);

  // Photocathode surface properties
  G4double photocath_EFF[]={0.25*perCent,36.246*perCent,39.8*perCent,40.0*perCent,36.0*perCent,30.0*perCent,
    24.0*perCent, 15.0*perCent,4.8*perCent}; //Enables 'detection' of photons
  assert(sizeof(photocath_EFF) == sizeof(ephotonPMT));
  G4double photocath_ReR[]={3.3817,3.3970,3.4392,3.5246,3.6880,4.0462,4.3354,3.6361,3.5688};
  assert(sizeof(photocath_ReR) == sizeof(ephotonPMT));
  G4double photocath_ImR[]={1.7862,1.7538,1.7171,0.8109, 0.565079, 0.41516, 0.3074, 0.252807,0.24902};
  assert(sizeof(photocath_ImR) == sizeof(ephotonPMT));

  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY",ephotonPMT,photocath_EFF,num);
  photocath_mt->AddProperty("REALRINDEX",ephotonPMT,photocath_ReR,num);
  photocath_mt->AddProperty("IMAGINARYRINDEX",ephotonPMT,photocath_ImR,num);

  G4OpticalSurface* photocath_opsurf= new G4OpticalSurface("photocath_opsurf"); // name, G4OpticalSurfaceModel, G4OpticalSurfaceFinish, G4SurfaceType, G4double value = 1.0
  photocath_opsurf->SetType(dielectric_metal);
  photocath_opsurf->SetFinish(polished);
  //photocath_opsurf->SetPolish(1.0);
  photocath_opsurf->SetModel(glisur);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  // Create Logical Skin surfaces (might need to change these to G4LogicalBorderSurfaces)
  new G4LogicalSkinSurface("photocath_surf", logicPC, photocath_opsurf); // name, physical volume of surface, phsical volume of world?, G4optical surface
  //
  // Air
  //
  G4double refractiveIndex2[] =
  { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries);

  //G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  //myMPT2->DumpTable();

  air->SetMaterialPropertiesTable(myMPT2);

  // Setup the sensitive detector
  PCSensitiveDetector* sensitive = new PCSensitiveDetector("PC");
  G4SDManager* sdman = G4SDManager::GetSDMpointer();
  sdman->AddNewDetector(sensitive);
  logicPC->SetSensitiveDetector(sensitive);


  //
  //always return the physical World!!!
  //
  return physWorld;
}

/* ************************************************************************************ */
