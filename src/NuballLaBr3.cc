//---------------------------------------------------------------------
// Create the solids defining an FATIMA LaBr3 detector
//---------------------------------------------------------------------
#include "NuballLaBr3.hh"

//Materials
#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "globals.hh"

#include "G4Transform3D.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// Definition of sensitive detectors

static const G4double inch = 2.54*cm;

LaBr3::LaBr3()
{

  // Definition of crystals/PM sizes
  fCrystalD = 1.5 * inch;
  fCrystalThickness = 2. * inch;

  fPMinnerD = (60.-2.)*mm;
  fPMouterD = (60.)*mm;
  fPMThickness = 200.*mm;

  flightguideThickness = 0.4 * cm;
  fCapD = fCrystalD + 2.* mm;
  fCapThickness = fCrystalThickness + flightguideThickness + 2*mm;

  // Calculation of radii
  fCrystalR = fCrystalD / 2.;
  fCapR = fCapD /2.;
  fPMinnerR = fPMinnerD/2.;
  fPMouterR = fPMouterD/2.;
  flightguideR = fCrystalR;

  //create the solids.....
  CreateSolids();

}

//Destructor
LaBr3::~LaBr3() { }

//---------------------------------------------------------------------
// Create the solids defining Phase-II LaBr3s
//---------------------------------------------------------------------
void  LaBr3::CreateSolids()
{
  // Calculation of crystal shift
  lappingSize = (fCrystalThickness-(fCapThickness-2*mm))/2.;
  LGlappingSize = (fCrystalThickness+ flightguideThickness)/2. + lappingSize;
  PMlappingSize = (fPMThickness)/2.;
  alulappingsize = fCapThickness/2.;

  //An approximate LaBr3

  //---------------------------------------------------------
  // Control Box
  //---------------------------------------------------------
  LaBr_solid = new G4Box("LaBr3_solid",
			 fPMouterR+0.1*mm,
			 fPMouterR+0.1*mm,
			 (fPMThickness+fCapThickness)/2+0.1*mm);           // Final Angle



  //---------------------------------------------------------
  // Aluminum CAP
  //---------------------------------------------------------
  Alcover_solid = new G4Tubs("LaBr3Alcover_solid",
			     0.,       // Inner radius
			     fCapR,           // Outer Radius
			     fCapThickness/2.,   // Height
			     0.,              // Initial Angle
			     360.);           // Final Angle


  //---------------------------------------------------------
  // LaBr3 crystal
  //---------------------------------------------------------
  LaBr3_solid = new G4Tubs("LaBr3Crystal_solid",
			   0.,                  // Inner radius
			   fCrystalR,           // Outer Radius
			   fCrystalThickness/2.,   // Height
			   0.,                  // Initial Angle
			   360.);               // Final Angle

  //---------------------------------------------------------
  // Light Guide
  //---------------------------------------------------------
  lightguide_solid = new G4Tubs("LaBr3LG_solid",
				0.,                  // Inner radius
				flightguideR,           // Outer Radius
				flightguideThickness/2.,   // Height
				0.,                  // Initial Angle
				360.);               // Final Angle


  //---------------------------------------------------------
  // Photomultiplier
  //---------------------------------------------------------
  PM_solid = new G4Tubs("LaBr3PM_solid",
			fPMinnerR,                  // Inner radius
			fPMouterR,           // Outer Radius
			fPMThickness/2.,   // Height
			0.,                  // Initial Angle
			360.);               // Final Angle



}

//------------------------------------------------------------------
void LaBr3::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps, G4ThreeVector position, G4RotationMatrix rotation, int i)
{

  //---------------------------------
 //make the required materials and assign them
 G4NistManager* nistManager = G4NistManager::Instance();

 //assign default materials.....
 G4Material* CapMaterial = nistManager->FindOrBuildMaterial("DurAl");

 G4Material* crystalMaterial = nistManager->FindOrBuildMaterial("LaBr3");

 G4Material* lightguideMaterial  = nistManager->FindOrBuildMaterial("Glass");

 G4Material* PMMaterial = lightguideMaterial;

 G4Material* ControleMaterial = nistManager->FindOrBuildMaterial("G4_AIR");

G4LogicalVolume* LaBr_log = new G4LogicalVolume(LaBr_solid,
        ControleMaterial,
        ("LaBr3_log"+std::to_string(i)).c_str(),
        0, 0, 0);

G4LogicalVolume* Alcover_log = new G4LogicalVolume(Alcover_solid,
           CapMaterial,
           ("LaBr3Alcover_log"+std::to_string(i)).c_str(),
           0, 0, 0);

G4LogicalVolume* LaBr3_log = new G4LogicalVolume(LaBr3_solid,
         crystalMaterial,
         ("LaBr3Crystal_log"+std::to_string(i)).c_str(),
         0, 0, 0);

G4LogicalVolume* lightguide_log= new G4LogicalVolume(lightguide_solid,
   lightguideMaterial,
   ("LaBr3LG_log"+std::to_string(i)).c_str(),
   0, 0, 0);

G4LogicalVolume* PM_log= new G4LogicalVolume(PM_solid,
           PMMaterial,
           ("LaBr3PM_log"+std::to_string(i)).c_str(),
           0, 0, 0);

//-----------------------------------------------------------------------//
 //                    Visual Attributes in the visualisation             //
 //-----------------------------------------------------------------------//
 //---------------------------------LaBr3 Al cover
 G4VisAttributes* CoverVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));   // gray
 //CoverVisAtt->SetForceSolid(true);
 CoverVisAtt ->SetForceWireframe(true);
 Alcover_log->SetVisAttributes(CoverVisAtt);

 G4VisAttributes* ControlVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));   // gray
 ControlVisAtt->SetVisibility(false);
 LaBr_log->SetVisAttributes(ControlVisAtt);


 //------------------------------------------LaBr3 crystal
 G4VisAttributes* LaBr3VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0));  // Magenta
 LaBr3VisAtt -> SetForceSolid(true);
 //LaBr3VisAtt -> SetForceWireframe(true);
 LaBr3_log -> SetVisAttributes(LaBr3VisAtt);

 G4VisAttributes* IonizationctrlAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  // cyan
 IonizationctrlAtt -> SetForceSolid(true);
 PM_log -> SetVisAttributes(IonizationctrlAtt);


 //-----------------------------------------------------------------------//
 //             Preparing all the "detectors" for tracking                //
 //-----------------------------------------------------------------------//
 fLogical.push_back(LaBr3_log);

 // Definition of the rotation matrix that might be used for placement
  G4RotationMatrix rmC(0,0,0);

auto LaBr_phys = new G4PVPlacement(G4Transform3D(rotation, position),
        "LaBr3_phys",                 //its name
        LaBr_log,     //its logical volume
        physiMother,        //its mother
        true,                       //no boolean operat
        copyNo,                     //copy number
        0);

  //=================================================================================
  // Definition of the Alu cap Log and phys
  //=================================================================================

  auto Alcover_phys = new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,alulappingsize-(fPMThickness+fCapThickness)/4)),
				   "LaBr3AlCaps_phys",                 //its name
				   Alcover_log,     //its logical volume
				   LaBr_phys,        //its mother
				   true,                       //no boolean operat
				   copyNo,                     //copy number
				   checkOverlaps);


  //=================================================================================
  // Definition of the LaBr3 crystal Log and phys
  //=================================================================================

  new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,-lappingSize)),
				 "LaBr3Crystal",
				 LaBr3_log,
				 Alcover_phys,
				 true,
				 copyNo,
				 checkOverlaps);


  //=================================================================================
  // Definition of the Light Guide Log and phys
  //=================================================================================

  new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,-LGlappingSize)),
				    "LaBr3LG_phys",
				    lightguide_log,
				    Alcover_phys,
				    true,
				    copyNo,
				    checkOverlaps);

  //=================================================================================
  // Definition of the Photomultiplier Log and phys
  //=================================================================================

  new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,-PMlappingSize-(fPMThickness+fCapThickness)/4)),
			    "LaBr3PM_phys",
			    PM_log,
			    LaBr_phys,
			    true,
			    copyNo,
			    checkOverlaps);
}
