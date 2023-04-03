//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: NuballDetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file NuballDetectorConstruction.cc
/// \brief Implementation of the NuballDetectorConstruction class

#include "NuballDetectorConstruction.hh"
#include "globals.hh"

#include "TMath.h"

//Materials
#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VUserDetectorConstruction.hh"

#include "NuballPhaseISD.hh"
#include "NuballCloverSD.hh"
#include "NuballLaBr3SD.hh"

#include "NuballClover.hh"
#include "NuballCloverBGO.hh"
#include "NuballPhaseI.hh"
#include "NuballPhaseIBGO.hh"
#include "NuballLaBr3.hh"
#include "NuballLaBr3Valencia.hh"
#include "LICORNE_Chamber.hh"
#include "OUPS_Chamber.hh"
//#include "NuballUTarget.hh"
#include "NuballThTarget.hh"

#include<fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuballDetectorConstruction::NuballDetectorConstruction()
: G4VUserDetectorConstruction(),
  fMakeLICORNE(false),
  fMakeOUPS(false),
  fMakeNEUTRONS(false),
  fMakePhaseI(true),
  fMakeClovers(true),
  fMakeLaBr3(true),
  fMakeFATIMA(false),
  fMakeTarget(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuballDetectorConstruction::~NuballDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* NuballDetectorConstruction::Construct()
{
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool Overlap = false;

  G4double coneangle = 20;

  // Constructing all the materials needed
  G4NistManager* nistManager = G4NistManager::Instance();

  G4cout << " Will Construct Materials" << G4endl;
  //Constructing DurAl
  std::vector<G4String> elDurAl = {"Cu", "Mg", "Mn", "Al"};
  std::vector<G4double> weightDurAl = {3.5*perCent, 0.5*perCent, 0.6*perCent, 95.4*perCent};
  nistManager->ConstructNewMaterial("DurAl", elDurAl, weightDurAl, 2.8*g/cm3);

  //Constructing Laboratory Vacuum
  G4double pressure = 1.e-5*pascal;  // 1e-6 mbar
  G4double temperature = (273.15+25.) * kelvin;
  nistManager->ConstructNewGasMaterial("LabVacuum", "G4_AIR", temperature, pressure, true);

  //Constructing BGO
  std::vector<G4String> elBGO = {"Bi", "Ge", "O"};
  std::vector<G4int> weightBGO = {4, 3, 12};
  nistManager->ConstructNewMaterial("BGO", elBGO, weightBGO, 7.13*g/cm3);

  //Constructing Hevimet
  std::vector<G4String> elHevimet = {"W", "Ni", "Fe"};
  std::vector<G4double> weightHevimet = {96.0*perCent, 3.5*perCent, 0.5*perCent};
  nistManager->ConstructNewMaterial("Hevimet", elHevimet, weightHevimet, 17.0*g/cm3);

  //Constructing LaBr3
  std::vector<G4String> elLaBr3 = {"La", "Br"};
  std::vector<G4int> weightLaBr3 = {1, 3};
  nistManager->ConstructNewMaterial("LaBr3", elLaBr3, weightLaBr3, 5.06*g/cm3);

  //Constructing Glass
  std::vector<G4String> elGlass = {"H", "C", "O"};
  std::vector<G4int> weightGlass = {8, 5, 2};
  nistManager->ConstructNewMaterial("Glass", elGlass, weightGlass, 1.18*g/cm3);

  //Constructing MgO
  std::vector<G4String> elMgO = {"Mg", "O"};
  std::vector<G4int> weightMgO = {1, 1};
  nistManager->ConstructNewMaterial("MgO", elMgO, weightMgO, 3.58*g/cm3);

  //Construction muMet
  std::vector<G4String> elMuMet = {"C", "Ni", "Mn", "Mo", "Si", "Fe"};
  std::vector<G4double> weightMuMet = {0.0002, 0.80, 0.005, 0.042, 0.0035, 0.1493};
  nistManager->ConstructNewMaterial("MuMet", elMuMet, weightMuMet, 8.747*g/cm3);

  G4cout << "Materials Created" << G4endl;


  G4RotationMatrix rotMat;
  G4double theta(0.), phi(0.);

  //Creating the list of detectors
  //fListOfDetectors = new DetectorRegister();

  //------------------------------------------------------------------//
  //                                                                  //
  //                 Experimental Setup environment                   //
  //                                                                  //
  //------------------------------------------------------------------//

  //------------------------------------------------------------------
  // World Volume
  //------------------------------------------------------------------
  world_size = 120.*cm;
  experimentalHall_solid = new G4Sphere("experimentalHall_solid", 0, world_size, 0, 360, 0, 360);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_solid, nistManager->FindOrBuildMaterial("G4_AIR"), "experimentalHall_log", 0, 0, 0);
  experimentalHall_phys = new G4PVPlacement(0, G4ThreeVector(), experimentalHall_log, "experimentalHall_phys", 0, false, 0);

  G4ThreeVector chamber_translation;
  if(fMakeLICORNE){
    //------------------------------------------------------------------
    // LICORNE Chamber
    //------------------------------------------------------------------
    // LICORNE cell's distance to the center
    chamberdistance = 21.9*cm;

    auto LICORNE = new LICORNE_Chamber();

    // Calculation of the distance between the center of LICORNE control volume and the sphere center
    cell_length = LICORNE -> GetCellThickness();
    chamber_length = LICORNE -> GetChamberThickness();
    translationdistance = chamberdistance + (chamber_length + cell_length)/2. ;

    //reset angles
    rotMat.set(0,0,0);
    rotMat.rotateY(theta);
    rotMat.rotateZ(phi);
    rotMat.invert();

    chamber_translation.set(0,0,-translationdistance);

    LICORNE->SetRotation(rotMat);
    LICORNE->SetPosition(chamber_translation);
    LICORNE->Placement(1, experimentalHall_phys, Overlap);
  }

  if(fMakeNEUTRONS)
  {
    //------------------------------------------------------------------
    // Neutron Cone
    //------------------------------------------------------------------
    G4double cone_outerradius = TMath::Tan(coneangle*TMath::Pi()/180.)*(chamberdistance+cell_length)*2;
    G4cout << cone_outerradius << " & " << chamberdistance+cell_length<< G4endl;
    neutron_cone_solid = new G4Cons("neutron_cone_solid", 0, 0, 0, cone_outerradius, (chamberdistance+cell_length)*2, 0, 360.);
    neutron_cone_log = new G4LogicalVolume(neutron_cone_solid, nistManager->FindOrBuildMaterial("G4_AIR"), "neutron_cone_log", 0, 0, 0);
    neutron_cone_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,(chamberdistance+cell_length)), "neutron_cone_phys", neutron_cone_log, experimentalHall_phys, false,0,0);
  }

  if(fMakeOUPS)
  {
    //------------------------------------------------------------------
    // OUPS Chamber
    //------------------------------------------------------------------
    // OUPS cell's distance to the center
    chamberdistance = 0*cm;

    OUPS_Chamber* OUPS = new OUPS_Chamber();

    //reset angles
    rotMat.set(0,0,0);
    rotMat.rotateY(theta);
    rotMat.rotateZ(phi);
    rotMat.invert();

    chamber_translation.set(0,0,-chamberdistance);

    OUPS->SetRotation(rotMat);
    OUPS->SetPosition(chamber_translation);
    OUPS->Placement(1, experimentalHall_phys, Overlap);

    G4cout << "Chamber created and placed " << G4endl;
  }

  //------------------------------------------------------------------//
  //                                                                  //
  //                        Target Geometry                           //
  //                                                                  //
  //------------------------------------------------------------------//
  if(fMakeTarget)
  {
    TargetTh* target = new TargetTh();
    rotMat.set(0,0,0);
    target->SetRotation(rotMat);
    target->SetPosition(G4ThreeVector(0, 0, 0));
    target->Placement(0, experimentalHall_phys, Overlap);
  }

  //------------------------------------------------------------------//
  //                                                                  //
  //                        Detection Setup                           //
  //                                                                  //
  //------------------------------------------------------------------//

  //------------------------------------------------------------------
  // Germanium detectors
  //------------------------------------------------------------------

  //-------------- Phase I detectors + BGOshield -------------------//

  if(fMakePhaseI)
  {
    PhaseI* phase1s = new PhaseI();
    PhaseIBGO* phase1bgos = new PhaseIBGO();

    std::fstream filePhaseI("PhaseIPosition.dat", std::ios_base::in);
    std::fstream filePhaseIBGO("PhaseIBGOPosition.dat", std::ios_base::in);
    if(!filePhaseI.is_open() || !filePhaseIBGO.is_open())
    {
      G4cerr << "ERROR: Input files for PhaseI geometry not found !" << G4endl;
      return 0;
    }

    for(int i = 0; i < 10; i++)
    {
      G4double rI, thetaI, phiI; //PhaseI coord + resolution parameters
      G4double rBGO, thetaBGO, phiBGO;//BGO coord + resolution parameters

      G4String PhaseIName, PhaseIAlveolus;
      G4String BGOName, BGOAlveolus;
      G4int PhaseILabel, BGOLabel;

      G4RotationMatrix rotMatp1;
      G4RotationMatrix rotMatBGO;
      filePhaseI >> PhaseIName >> PhaseIAlveolus >>  rI >> thetaI >> phiI >> PhaseILabel;
      filePhaseIBGO >> BGOName >> BGOAlveolus >> rBGO >> thetaBGO >> phiBGO >> BGOLabel;

      thetaI*=deg;
      phiI*=deg;
      thetaBGO*=deg;
      phiBGO*=deg;


      //--Phase I--//
      rotMatp1.set(0,0,0);
      rotMatp1.rotateY(thetaI);
      rotMatp1.rotateZ(phiI);

      G4double radius = (rI + 23.8)*mm;
      //G4double radius = rI*mm;
      G4ThreeVector translation(radius*sin(thetaI)*cos(phiI),
      radius*sin(thetaI)*sin(phiI),
      radius*cos(thetaI));

      //phase1s->SetRotation(rotMatp1);
      //phase1s->SetPosition(translation);
      phase1s->Placement(PhaseILabel, experimentalHall_phys, Overlap, translation, rotMatp1, i);

      //--BGO--//
      rotMatBGO.set(0,0,0);
      //rotMatBGO.set(10*i*deg,0,0);
      rotMatBGO.rotateZ(((i%2)*18.5*deg+(i%2-1)*198.5*deg)); //to try to reduce the overlap between attenant BGOs
      rotMatBGO.rotateY(thetaBGO);
      rotMatBGO.rotateZ(phiBGO);

      //rotMatBGO.rotateX(phiBGO);
      //rotMatBGO.rotateY(thetaBGO);
      ///rotMatBGO.invert();
      //G4double sradius = (rBGO - 35.)*mm;
      G4double sradius = (rBGO+23.8)*mm;

      G4ThreeVector stranslation(sradius*sin(thetaBGO)*cos(phiBGO),
      sradius*sin(thetaBGO)*sin(phiBGO),
      sradius*cos(thetaBGO));
      //G4ThreeVector stranslation(0, 0, sradius);
      //stranslation.rotateX(-phiBGO);
      //stranslation.rotateY(-thetaBGO);

      phase1bgos->SetHeavyMet(false);
      //phase1bgos->SetRotation(rotMatBGO);
      //phase1bgos->SetPosition(stranslation);
      phase1bgos->Placement(BGOLabel, experimentalHall_phys, Overlap, stranslation, rotMatBGO, i);
    }
    filePhaseI.close();
    filePhaseIBGO.close();
  }

  //-------------- Clover detectors + BGOshield--------------------//
  if(fMakeClovers)
  {
    Clover* clovers = new Clover();
    CloverBGO* cloverbgos = new CloverBGO();
    std::fstream fileClover("CloverPosition.dat", std::ios_base::in);
    std::fstream fileCloverBGO("CloverBGOPosition.dat", std::ios_base::in);

    if(!fileClover.is_open() || !fileCloverBGO.is_open())
    {
      G4cerr << "ERROR: Input files for clover geometry not found !" << G4endl;
      return 0;
    }
    for(int i = 0; i < 24; i++)
    {
      G4double rC, thetaC, phiC;
      G4double rBGO, thetaBGO, phiBGO;
      G4String cloverName, cloverAlveolus;
      G4String cloverBGOName, cloverBGOAlveolus;
      G4int cloverLabel, cloverBGOLabel;

      rotMat.set(0,0,0);
      fileClover >> cloverName >> cloverAlveolus >> rC >> thetaC >> phiC >> cloverLabel;
      fileCloverBGO >> cloverBGOName >> cloverBGOAlveolus >> rBGO >> thetaBGO >> phiBGO >> cloverBGOLabel;

      thetaC*=deg;
      phiC*=deg;
      thetaBGO*=deg;
      phiBGO*=deg;
      //--Clover--//
      rotMat.rotateY(thetaC);
      rotMat.rotateZ(phiC);
      G4double radius = (rC + 34.5)*mm;

      G4ThreeVector translation(radius*sin(thetaC)*cos(phiC),
      radius*sin(thetaC)*sin(phiC),
      radius*cos(thetaC));

      clovers->Placement(cloverLabel, experimentalHall_phys, Overlap, translation, rotMat, i);

      //--BGO--//
      rotMat.set(0,0,0);
      rotMat.rotateZ(-45*deg+180*deg*(i%2-1));
      rotMat.rotateY(thetaBGO);
      rotMat.rotateZ(phiBGO);

      //G4double sradius = (rBGO - 35.)*mm;
      G4double sradius = (rBGO)*mm;
      G4ThreeVector stranslation(sradius*sin(thetaBGO)*cos(phiBGO),
      sradius*sin(thetaBGO)*sin(phiBGO),
      sradius*cos(thetaBGO));

      cloverbgos->SetHeavyMet(false);
      cloverbgos->Placement(cloverBGOLabel, experimentalHall_phys, Overlap, stranslation, rotMat, i);

    }
    fileClover.close();
    fileCloverBGO.close();
  }

  if(fMakeLaBr3)
  {
    LaBr3_valencia* LaBrValencia = new LaBr3_valencia();
    LaBr3* LaBr = new LaBr3();
    std::fstream file("LaBr3Position.dat", std::ios_base::in);

    if(!file.is_open())
    {
      G4cerr << "ERROR: Input files for LaBr3 geometry not found !" << G4endl;
      return 0;
    }


    G4double r1(0.),theta1(0.),phi1(0.);
    G4double r2(0.),theta2(0.),phi2(0.);
    for(int i = 0; i < 10; i++)
    {
      G4String valenciaName, valenciaAlveolus;
      G4String UKName, UKAlveolus;
      G4int valenciaLabel, UKLabel;

      rotMat.set(0,0,0);

      file >> valenciaName >> valenciaAlveolus >> r1 >> theta1 >> phi1 >> valenciaLabel;
      file >> UKName >> UKAlveolus >> r2 >> theta2 >> phi2 >> UKLabel;

      theta1*=deg;
      phi1*=deg;
      theta2*=deg;
      phi2*=deg;

      G4double radius1 = (r1+128.75)*mm;
      G4double radius2 = r2*mm;

      //G4cout << "r = " << r << " theta = " << theta << " phi = " << phi << G4endl;
      rotMat.rotateY(theta1);
      rotMat.rotateZ(phi1);

      G4ThreeVector translation1(radius1*sin(theta1)*cos(phi1),
      radius1*sin(theta1)*sin(phi1),
      radius1*cos(theta1));

      LaBrValencia->Placement(valenciaLabel, experimentalHall_phys, Overlap, translation1, rotMat, i);

      rotMat.set(0,0,0);
      rotMat.rotateX(180*deg);
      rotMat.rotateY(theta1); //same orientation than 1
      rotMat.rotateZ(phi1);
      //rotMat.invert();

      G4ThreeVector translation2(radius2*sin(theta2)*cos(phi2),
      radius2*sin(theta2)*sin(phi2),
      radius2*cos(theta2));

      LaBr->Placement(UKLabel, experimentalHall_phys, Overlap, translation2, rotMat,i);

    }
    file.close();
  }


  //-----------------------------------------------------------------------//
  //                    Visual Attributes in the visualisation             //
  //-----------------------------------------------------------------------//

  //------------------------------------------Experimental Hall
  G4VisAttributes* WorldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));  //
  //WorldVisAtt -> SetForceSolid(true);
  WorldVisAtt -> SetVisibility(false);
  experimentalHall_log -> SetVisAttributes(WorldVisAtt);

  // Print all materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4cout << "Exiting NuballDetectorConstruction::Construct()" << G4endl;

  return experimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix NuballDetectorConstruction::StringToRotationMatrix(G4String rotation)
{
  // We apply successive rotations OF THE OBJECT around the FIXED
  // axes of the parent's local coordinates; rotations are applied
  // left-to-right (rotation="r1,r2,r3" => r1 then r2 then r3).

  G4RotationMatrix rot;

  unsigned int place = 0;

  while (place < rotation.size())
  {

    G4double angle;
    char* p;

    const G4String tmpstring=rotation.substr(place+1);

    angle = strtod(tmpstring.c_str(),&p) * deg;

    if (!p || (*p != (char)',' && *p != (char)'\0'))
    {
      G4cerr << "Invalid rotation specification: " <<
      rotation.c_str() << G4endl;
      return rot;
    }

    G4RotationMatrix thisRotation;

    switch(rotation.substr(place,1).c_str()[0])
    {
      case 'X': case 'x':
      thisRotation = G4RotationMatrix(CLHEP::HepRotationX(angle));
      break;
      case 'Y': case 'y':
      thisRotation = G4RotationMatrix(CLHEP::HepRotationY(angle));
      break;
      case 'Z': case 'z':
      thisRotation = G4RotationMatrix(CLHEP::HepRotationZ(angle));
      break;
      default:
      G4cerr << " Invalid rotation specification: "
      << rotation << G4endl;
      return rot;
    }

    rot = thisRotation * rot;
    place = rotation.find(',',place);
    if (place > rotation.size()) break;
    ++place;
  }

  return rot;
}

void NuballDetectorConstruction::ConstructSDandField()
{
  if(fMakePhaseI)
  {
  NuballPhaseISD* fPhaseISD = new NuballPhaseISD("PhaseISD", "PhaseIHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(fPhaseISD);
  for(int i = 0; i < 10; i++)
  {
    SetSensitiveDetector(("PhaseIGe"+std::to_string(i)).c_str(), fPhaseISD);
    SetSensitiveDetector(("PhaseIBGOCrystal"+std::to_string(i)).c_str(), fPhaseISD);
  }
}
if(fMakeClovers)
{
  NuballCloverSD* fCloverSD = new NuballCloverSD("CloverSD", "CloverHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(fCloverSD);
  for(int i = 0; i < 24; i++)
  {
    SetSensitiveDetector(("CloverLeaf"+std::to_string(i)).c_str(), fCloverSD);
    SetSensitiveDetector(("CloverBGOCrystal"+std::to_string(i)).c_str(), fCloverSD);
  }
}
  if(fMakeLaBr3)
  {
  NuballLaBr3SD* fLaBr3SD = new NuballLaBr3SD("LaBr3SD", "LaBr3HitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(fLaBr3SD);


  for(int i = 0; i < 10; i++)
  {
    SetSensitiveDetector(("LaBr3Crystal_log"+std::to_string(i)).c_str(), fLaBr3SD);
    SetSensitiveDetector(("LaBr3ValenciaCrystal_log"+std::to_string(i)).c_str(), fLaBr3SD);
  }

}
}
