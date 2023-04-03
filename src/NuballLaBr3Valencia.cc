//---------------------------------------------------------------------
// Create the solids defining a FATIMA LaBr3 detector
//---------------------------------------------------------------------
#include "NuballLaBr3Valencia.hh"

//Materials
#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "TMath.h"

static const G4double inch = 2.54*cm;

LaBr3_valencia::LaBr3_valencia():fLogical(0)
{

  // Definition of crystal parameters
  fsmallerdiameter   = 25.4* mm; //19
  fgreaterdiameter   = 38.1*mm; //38
  fconethickness     = 38.1*mm;//30.5?

  //Definition of Aluminium Cap
  fcapsmallerdiameter    = 32.25*mm;// was 32.25 mm
  fcapgreaterdiameter    = 48.0*mm;// was 48.0 mm
  fcapconethickness      = 45.5*mm; // was 45.5 mm
  ffrontthickness        = 0.5*mm;

  // Definition of Reflector parameters
  fRsmallerdiameter   = (fcapsmallerdiameter-0.5*2)* mm;
  fRgreaterdiameter   = (fcapgreaterdiameter-1.3*2)*mm;
  fRconethickness     = (fconethickness+0.625)*mm;

  //Definition of Lead Shielding
  fshieldinnersmallerdiameter    = (fcapsmallerdiameter+0.1)*mm;
  fshieldsmallerdiameter    = (fcapsmallerdiameter+2.2)*mm;
  fshieldinnergreaterdiameter    = (fcapgreaterdiameter+0.1)*mm;
  fshieldgreaterdiameter    = (fcapgreaterdiameter+2.2)*mm;
  fshieldconethickness      = (fcapconethickness+1.0)*mm;

  // Definition of Light Guide parameters
  flightguideThickness = 5.*mm; // was 6.5 mm
  flightguideDiameter  = fgreaterdiameter;

  // Definition of PM coating parameters
  fPMcoatD = (60.+2*0.8)*mm;
  fPMcoatThickness = 200.*mm; // was 180 mm

  // Definition of PM parameters
  fPMinnerD = (60.-2.)*mm;
  fPMouterD = (60.)*mm;
  fPMThickness = 200.*mm; // was 180 mm

  // Definition of crystals/PM sizes
  fCrystalD            = fsmallerdiameter;
  fCrystalR            = fCrystalD/2.;
  fCrystalThickness    = fconethickness;
  fCapD                = fcapsmallerdiameter;
  fCapThickness        = fcapconethickness;


  // Definition of Control Volume
  fcontroldiameter  = fPMcoatD+1*mm;
  fcontrolthickness = fcapconethickness + fPMThickness+1.5*mm;

  // Calculation of radii
  fCrystalR = fCrystalD / 2.;
  fCapR = fCapD /2.;
  fPMR = fPMcoatD/2.;

  //create the solids.....
  CreateSolids();

}

//Destructor
LaBr3_valencia::~LaBr3_valencia() { }

//---------------------------------------------------------------------
// Create the solids defining Phase-II LaBr3_valencias
//---------------------------------------------------------------------
void  LaBr3_valencia::CreateSolids()
{


  // Calculation of crystal shift
  lappingSize = 6*mm;
  Refllappingsize = flightguideThickness/2.-0.15;
  LGlappingSize = (fCrystalThickness/2 + flightguideThickness/2.-ffrontthickness*5.5)*mm;
  PMlappingSize = (fcapconethickness+0.5*mm)/2.;
  alulappingsize = (fPMThickness+0.5*mm)/2.;

  G4RotationMatrix rotMat;
  //set up the angles for the detectors...
  rotMat.set(0,0,0);

  //An approximate LaBr3
  //G4cout << " ----> Constructing archetypal LaBr3" << G4endl;

  //---------------------------------------------------------
  // Control Box
  //---------------------------------------------------------
  LaBr_solid = new G4Tubs("LaBr3Valencia_solid",
			  0.,                  // Inner radius
			  fcontroldiameter/2.,           // Outer Radius
			  fcontrolthickness/2.,   // Height
			  0.,                  // Initial Angle
			  360.*deg);               // Final Angle


  //---------------------------------------------------------
  // CAP Lead shield
  //---------------------------------------------------------
  Pbcover_cone = new G4Cons("LaBr3ValenciaPbcover_cone",
			    0.,
			    fshieldsmallerdiameter/2.,
			    0.,
			    fshieldgreaterdiameter/2.,
			    fshieldconethickness/2.,
			    0,
			    2*TMath::Pi());

  //HepPolyhedronCone(s)/Tube(s): error in input parameters (radiuses)
  ///Rmn1=16.125 Rmx1=16.13 Rmn2=24 Rmx2=17.13 Dz=22.75 Phi1=0 Dphi=6.28319

  //---------------------------------------------------------
  // Aluminum CAP
  //---------------------------------------------------------
  Alcover_cone = new G4Cons("LaBr3ValenciaAlcover_cone",
			    0.,
			    fcapsmallerdiameter/2.,
			    0.,
			    fcapgreaterdiameter/2.,
			    fcapconethickness/2.,
			    0,
			    2*TMath::Pi());

  Alcover_front = new G4Tubs("LaBr3ValenciaAlcover_front",
			     0.,
			     fcapsmallerdiameter/2.,
			     ffrontthickness/2.,   // Height
			     0.,              // Initial Angle
			     360.*deg);           // Final Angle

  Alcover_solid  = new G4UnionSolid("LaBr3ValenciaAlcover_solid1",
				    Alcover_cone,
				    Alcover_front,
				    &rotMat,
				    G4ThreeVector(0,0,-(ffrontthickness/2.+fcapconethickness/2.)));


  //---------------------------------------------------------
  // Reflector Cone
  //---------------------------------------------------------
  Refl_solid = new G4Cons("LaBr3ValenciaAir_cone",
			  0.,
			  fRsmallerdiameter/2.,
			  0.,
			  fRgreaterdiameter/2.,
			  fRconethickness/2.,
			  0,
			  2*TMath::Pi());

  //---------------------------------------------------------
  // LaBr3 crystal
  //---------------------------------------------------------
  LaBr3_solid = new G4Cons("LaBr3ValenciaCrystal_cone",
			   0,
			   fsmallerdiameter/2.,
			   0,
			   fgreaterdiameter/2.,
			   fconethickness/2.,
			   0,
			   2*TMath::Pi());

  //---------------------------------------------------------
  // Light Guide
  //---------------------------------------------------------
  lightguide_solid = new G4Tubs("LaBr3ValenciaLG_solid",
				0.,                  // Inner radius
				flightguideDiameter/2.,           // Outer Radius
				flightguideThickness/2.,   // Height
				0.,                  // Initial Angle
				360.*deg);               // Final Angle

  //---------------------------------------------------------
  // Photomultiplier mumetal coat
  //---------------------------------------------------------
  PMCap_solid = new G4Tubs("LaBr3ValenciaPMcoat_solid",
			   0.,              // Inner radius
			   fPMcoatD/2.,     // Outer Radius
			   fPMcoatThickness/2., // Height
			   0.,              // Initial Angle
			   360.*deg);           // Final Angle

  //---------------------------------------------------------
  // Photomultiplier
  //---------------------------------------------------------
  PM_solid = new G4Tubs("LaBr3ValenciaPM_solid",
			fPMinnerD/2.,        // Inner radius
			fPMouterD/2.,           // Outer Radius
			fPMThickness/2.,   // Height
			0.,                // Initial Angle
			360.*deg);             // Final Angle



}
//------------------------------------------------------------------
void LaBr3_valencia::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps, G4ThreeVector position, G4RotationMatrix rotation, int i)
{

  //make the required materials and assign them
    G4NistManager* nistManager = G4NistManager::Instance();

    //assign default materials.....
    G4Material* Reflectormaterial = nistManager->FindOrBuildMaterial("MgO");

    G4Material* Shieldmaterial = nistManager->FindOrBuildMaterial("G4_Pb");

    G4Material* PMcoatmaterial = nistManager->FindOrBuildMaterial("MuMet");

    G4Material* CapMaterial = nistManager->FindOrBuildMaterial("DurAl");

    G4Material* crystalMaterial = nistManager->FindOrBuildMaterial("LaBr3");

    G4Material* lightguideMaterial  = nistManager->FindOrBuildMaterial("Glass");

    G4Material* PMMaterial = lightguideMaterial;

    G4Material* ControleMaterial  = nistManager->FindOrBuildMaterial("G4_AIR");

    //=================================================================================
    // Definition of the control box Log and phys
    //=================================================================================
    G4LogicalVolume* LaBr_log = new G4LogicalVolume(LaBr_solid,
  				 ControleMaterial,
  				 ("LaBr3Valencia_log"+std::to_string(i)).c_str(),
  				 0, 0, 0);

    //=================================================================================
    // Definition of the Lead Shield Log and phys
    //=================================================================================
    G4LogicalVolume* Shieldcover_log = new G4LogicalVolume(Pbcover_cone,
  					Shieldmaterial,
  					("LaBr3ValenciaShieldcover_log"+std::to_string(i)).c_str(),
  					0, 0, 0);

    //=================================================================================
    // Definition of the Alu cap Log and phys
    //=================================================================================
    G4LogicalVolume* Alcover_log = new G4LogicalVolume(Alcover_solid,
  				    CapMaterial,
  				    ("LaBr3ValenciaAlcover_log"+std::to_string(i)).c_str(),
  				    0, 0, 0);

    //=================================================================================
    // Definition of the Alu cap Log and phys
    //=================================================================================
    G4LogicalVolume* Refl_log = new G4LogicalVolume(Refl_solid,
  				 Reflectormaterial,
  				 ("LaBr3ValenciaRefl_log"+std::to_string(i)).c_str(),
  				 0, 0, 0);

    //=================================================================================
    // Definition of the LaBr3 crystal Log and phys
    //=================================================================================
    G4LogicalVolume* LaBr3Val_log = new G4LogicalVolume(LaBr3_solid,
  				  crystalMaterial,
  				  ("LaBr3ValenciaCrystal_log"+std::to_string(i)).c_str(),
  				  0, 0, 0);

    //=================================================================================
    // Definition of the Light Guide Log and phys
    //=================================================================================
    G4LogicalVolume* lightguide_log= new G4LogicalVolume(lightguide_solid,
  				      lightguideMaterial,
  				      ("LaBr3ValenciaLG_log"+std::to_string(i)).c_str(),
  				      0, 0, 0);
    //=================================================================================
    // Definition of the PM mu metal coating Log and phys
    //=================================================================================
    G4LogicalVolume* PMcoat_log= new G4LogicalVolume(PMCap_solid,
  				  PMcoatmaterial,
  				  ("LaBr3ValenciaPMcoat_log"+std::to_string(i)).c_str(),
  				  0, 0, 0);


    //=================================================================================
    // Definition of the Photomultiplier Log and phys
    //=================================================================================
    G4LogicalVolume* PM_log= new G4LogicalVolume(PM_solid,
  			      PMMaterial,
  			      ("LaBr3ValenciaPM_log"+std::to_string(i)).c_str(),
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
    //Refl_log->SetVisAttributes(ControlVisAtt);
    PMcoat_log->SetVisAttributes(ControlVisAtt);

    //------------------------------------------LaBr3 crystal
    G4VisAttributes* LaBr3VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0));  // Magenta
    LaBr3VisAtt -> SetForceSolid(true);
    //LaBr3VisAtt -> SetForceWireframe(true);
    LaBr3Val_log -> SetVisAttributes(LaBr3VisAtt);

    G4VisAttributes* IonizationctrlAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  // cyan
    IonizationctrlAtt -> SetForceSolid(true);
    PM_log -> SetVisAttributes(IonizationctrlAtt);


    G4VisAttributes* ShieldAtt= new G4VisAttributes(G4Colour(0.74,0.74,0.74));
    ShieldAtt -> SetForceWireframe(true);//SetForceSolid(true);
    //ShieldAtt->SetVisibility(false);
    Shieldcover_log -> SetVisAttributes(ShieldAtt);

    G4VisAttributes* CtrlAtt= new G4VisAttributes(G4Colour(1.0,0.0,0));
    CtrlAtt -> SetForceSolid(true);
    CtrlAtt ->SetForceWireframe(true);
    Refl_log -> SetVisAttributes(CtrlAtt);

    //-----------------------------------------------------------------------//
    //             Preparing all the "detectors" for tracking                //
    //-----------------------------------------------------------------------//
    fLogical.push_back(LaBr3Val_log);
    //---------------------------------

  // Definition of the rotation matrix that might be used for placement
  G4RotationMatrix rmC(0,0,0);

  auto LaBr_phys = new G4PVPlacement(G4Transform3D(rotation, position),
				"LaBr3Valencia_phys",                 //its name
				LaBr_log,     //its logical volume
				physiMother,        //its mother
				true,                       //no boolean operat
				copyNo,                     //copy number
				0);

  auto Shieldcover_phys = new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,-alulappingsize)),
				       "LaBr3ValenciaShieldcover_phys",                 //its name
				       Shieldcover_log,     //its logical volume
				       LaBr_phys,        //its mother
				       true,                       //no boolean operat
				       copyNo,                     //copy number
				       checkOverlaps);

  auto Alcover_phys = new G4PVPlacement(G4Transform3D(rmC,G4ThreeVector(0,0,0)),
				   "LaBr3ValenciaAlcover_phys",                 //its name
				   Alcover_log,     //its logical volume
				   Shieldcover_phys,        //its mother
				   true,                       //no boolean operat
				   copyNo,                     //copy number
				   0);


  auto Refl_phys = new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,-Refllappingsize)),
				"LaBr3ValenciaRefl_phys",                 //its name
				Refl_log,     //its logical volume
				Alcover_phys,        //its mother
				true,                       //no boolean operat
				copyNo,                     //copy number
				0);

  new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,0)),
				 "LaBr3ValenciaCrystal_phys",
				 LaBr3Val_log,
				 Refl_phys,
				 true,
				 copyNo,
				 checkOverlaps);

  new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,LGlappingSize)),
				    "LaBr3ValenciaLG_phys",
				    lightguide_log,
				    Alcover_phys,
				    true,
				    copyNo,
				    0);


  auto PMcoat_phys=new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,+PMlappingSize)),
				"LaBr3ValenciaPMcoat_phys",
				PMcoat_log,
				LaBr_phys,
				true,
				copyNo,
				checkOverlaps);


  new G4PVPlacement(G4Transform3D(rmC, G4ThreeVector(0,0,0)),
			    "LaBr3ValenciaPM_phys",
			    PM_log,
			    PMcoat_phys,
			    true,
			    copyNo,
			    0);

}
