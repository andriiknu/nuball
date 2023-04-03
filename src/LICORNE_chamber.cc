//---------------------------------------------------------------------
// Create the solids defining an FATIMA LICORNE_Chamber detector
//---------------------------------------------------------------------
#include "LICORNE_Chamber.hh"

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

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

static const G4double inch = 2.54*cm;

LICORNE_Chamber::LICORNE_Chamber()
{

    // LICORNE chamber parameters
    chamber_diameter = 17.*cm;
    chamber_length = 17.*cm;
    gascell_diameter = 2.*cm;
    gascell_length   = 3.5*cm;
    shielding_length=gascell_length+5.*mm;
    shielding_diameter= gascell_diameter + 3.*cm;
    shielding_innerdiameter = gascell_diameter+0.5*cm;
    control_length = chamber_length + shielding_length + 0.5*mm;
    control_diameter = chamber_diameter + 2.*mm;

    //---------------------------------
    //Definition of the materials
    G4NistManager* nistManager = G4NistManager::Instance();

    //assign default materials.....
    ControlMaterial        = nistManager->FindOrBuildMaterial("G4_AIR");

    G4Material* DurAl = new G4Material("DurAl", 2.8*g/cm3, 4);
    DurAl->AddElement(nistManager->FindOrBuildElement("G4_Cu"), 3.5*perCent);
    DurAl->AddElement(nistManager->FindOrBuildElement("G4_Mg"), 0.5*perCent);
    DurAl->AddElement(nistManager->FindOrBuildElement("G4_Mn"), 0.6*perCent);
    DurAl->AddElement(nistManager->FindOrBuildElement("G4_Al"), 95.4*perCent);

    G4double pressure = 1.e-5*pascal;  // 1e-6 mbar
    G4double density  = 1.2928*mg/cm3 * pressure / (1.01325 * 100.e3*pascal);
    G4double temperature = (273.15+25.) * kelvin;
    G4Material* vacuumMaterial = new G4Material("Vacuum", density, 1, kStateGas, temperature, pressure);
    vacuumMaterial->AddMaterial(nistManager->FindOrBuildMaterial("G4_AIR"), 1.);

    ChamberMaterial        = DurAl;
    innerChamberMaterial   = vacuumMaterial;
    cellMaterial           = DurAl;

    // DiHydrogen gas at 1.15bar (Target chamber volume)
    pressure = 115000*pascal; // 100kPa = 1 bar
    temperature = (273.15+20.)*kelvin;
    density = (pressure*2.0159*g/mole)/(8.3145*temperature)*mg/cm3;
    G4Material* HGas = new G4Material("HGas",density,1,kStateGas,temperature,pressure);
    HGas->AddElement(nistManager->FindOrBuildElement("G4_H"),2);

    innercellMaterial      = HGas;
    innershieldingMaterial = nistManager->FindOrBuildMaterial("G4_Pb");
    shieldingMaterial      = nistManager->FindOrBuildMaterial("G4_Cu");


    //create the solids.....
    CreateSolids();

}

//Destructor
LICORNE_Chamber::~LICORNE_Chamber() { }

void LICORNE_Chamber::SetPosition(G4ThreeVector thisPos) {
    position = thisPos*mm;
    G4cout << " ----> A LICORNE_Chamber will be placed at " << position/mm << " mm" << G4endl;
}

void LICORNE_Chamber::SetRotation(G4RotationMatrix thisRot) { rotation = thisRot; }

//---------------------------------------------------------------------
// Create the solids defining Phase-II LICORNE_Chambers
//---------------------------------------------------------------------
void  LICORNE_Chamber::CreateSolids()
{

    //An approximate LICORNE_Chamber
    G4cout << " ----> Constructing archetypal LICORNE_Chamber" << G4endl;

    //------------ CONTROL VOLUME ------------//
    // Definition of the geometric volume
    controlchamber_solid = new G4Tubs("controlchamber_solid",
                                      0.,                  // Inner radius
                                      control_diameter/2., // Outer Radius
                                      control_length/2.,   // Height
                                      0.,                  // Initial Angle
                                      360.);               // Final Angle

    //------------ LICORNE CHAMBER ------------//
    // Definition of the geometric volume
    chamber_solid = new G4Tubs("chamber_solid",
                               0.,                  // Inner radius
                               chamber_diameter/2., // Outer Radius
                               chamber_length/2.,   // Height
                               0.,                  // Initial Angle
                               360.);               // Final Angle

    innerchamber_solid = new G4Tubs("innerchamber_solid",
                                    0.,                       // Inner radius
                                    chamber_diameter/2.-1*cm, // Outer Radius
                                    chamber_length/2.-1*cm,   // Height
                                    0.,                       // Initial Angle
                                    360.);                    // Final Angle

    //------------ LICORNE GAS Cell ------------//
    // Definition of the geometric volume
    cell_solid = new G4Tubs("cell_solid",
                            0.,                  // Inner radius
                            gascell_diameter/2., // Outer Radius
                            gascell_length/2.,   // Height
                            0.,                  // Initial Angle
                            360.);               // Final Angle

    innercell_solid = new G4Tubs("innercell_solid",
                                 0.,                       // Inner radius
                                 gascell_diameter/2.-1*mm, // Outer Radius
                                 gascell_length/2.-1*mm,   // Height
                                 0.,                       // Initial Angle
                                 360.);                    // Final Angle


    //------------ LICORNE Shielding ------------//
    // Definition of the geometric volume
    shielding_solid = new G4Tubs("shielding_solid",
                                 shielding_innerdiameter/2.,                  // Inner radius
                                 shielding_diameter/2., // Outer Radius
                                 shielding_length/2.,   // Height
                                 0.,                  // Initial Angle
                                 360.);               // Final Angle

    innershielding_solid = new G4Tubs("innershielding_solid",
                                      shielding_innerdiameter/2.+0.5*mm,                       // Inner radius
                                      shielding_diameter/2.-.5*mm, // Outer Radius
                                      shielding_length/2.-.5*mm,   // Height
                                      0.,                       // Initial Angle
                                      360.);                    // Final Angle


}


//------------------------------------------------------------------
void LICORNE_Chamber::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

    //------------ CONTROL VOLUME ------------//
    // Definition of the control logical and physical volumes
    controlchamber_log = new G4LogicalVolume(controlchamber_solid,
                                             ControlMaterial,
                                             "controlchamber_log",
                                             0, 0, 0);
    controlchamber_phys = new G4PVPlacement(&rotation,
                                            G4ThreeVector(position.x(),position.y(),position.z()),
                                            "controlchamber_phys",  //its name
                                            controlchamber_log,     //its logical volume
                                            physiMother,  //its mother
                                            true,                   //no boolean operat
                                            copyNo,                      //copy number
                                            checkOverlaps);


    //------------ LICORNE CHAMBER ------------//
    // Definition of the control logical and physical volumes
    chamber_log = new G4LogicalVolume(chamber_solid,
                                      ChamberMaterial,
                                      "chamber_log",
                                      0, 0, 0);

    chamber_phys = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,-shielding_length/2.),
                                     "chamber_phys",  //its name
                                     chamber_log,     //its logical volume
                                     controlchamber_phys,  //its mother
                                     true,                   //no boolean operat
                                     copyNo,                      //copy number
                                     checkOverlaps);

    // Definition of the control logical and physical volumes
    innerchamber_log = new G4LogicalVolume(innerchamber_solid,
                                           innerChamberMaterial,
                                           "innerchamber_log",
                                           0, 0, 0);

    innerchamber_phys = new G4PVPlacement(0,
                                          G4ThreeVector(0,0,0),
                                          "innerchamber_phys",  //its name
                                          innerchamber_log,     //its logical volume
                                          chamber_phys,         //its mother
                                          true,                 //no boolean operat
                                          copyNo,                    //copy number
                                          checkOverlaps);

    //------------ LICORNE GAS Cell ------------//
    // Definition of the control logical and physical volumes
    cell_log = new G4LogicalVolume(cell_solid,
                                   cellMaterial,
                                   "cell_log",
                                   0, 0, 0);

    cell_phys = new G4PVPlacement(0,
                                  G4ThreeVector(0,0,chamber_length/2.),
                                  "cell_phys",  //its name
                                  cell_log,     //its logical volume
                                  controlchamber_phys,  //its mother
                                  true,                   //no boolean operat
                                  copyNo,                      //copy number
                                  checkOverlaps);

    innercell_log = new G4LogicalVolume(innercell_solid,
                                        innercellMaterial,
                                        "innercell_log",
                                        0, 0, 0);

    innercell_phys = new G4PVPlacement(0,
                                       G4ThreeVector(0,0,0),
                                       "innercell_phys",  //its name
                                       innercell_log,     //its logical volume
                                       cell_phys,  //its mother
                                       true,                   //no boolean operat
                                       copyNo,                      //copy number
                                       checkOverlaps);

    //------------ LICORNE GAS Cell Shielding ------------//
    // Definition of the control logical and physical volumes
    shielding_log = new G4LogicalVolume(shielding_solid,
                                        shieldingMaterial,
                                        "shielding_log",
                                        0, 0, 0);

    shielding_phys = new G4PVPlacement(0,
                                       G4ThreeVector(0,0,chamber_length/2.),
                                       "shielding_phys",  //its name
                                       shielding_log,     //its logical volume
                                       controlchamber_phys,  //its mother
                                       true,                   //no boolean operat
                                       copyNo,                      //copy number
                                       checkOverlaps);

    innershielding_log = new G4LogicalVolume(innershielding_solid,
                                             innershieldingMaterial,
                                             "innershielding_log",
                                             0, 0, 0);

    innershielding_phys = new G4PVPlacement(0,
                                            G4ThreeVector(0,0,0),
                                            "innershielding_phys",  //its name
                                            innershielding_log,     //its logical volume
                                            shielding_phys,  //its mother
                                            true,                   //no boolean operat
                                            copyNo,                      //copy number
                                            checkOverlaps);



    //-----------------------------------------------------------------------//
    //                    Visual Attributes in the visualisation             //
    //-----------------------------------------------------------------------//
    //---------------------------------LICORNE_Chamber Control volume
    G4VisAttributes* ControlVisAtt= new G4VisAttributes(G4Colour(0.66,0.66,0.66));   // gray
    ControlVisAtt->SetVisibility(false);
    controlchamber_log -> SetVisAttributes(ControlVisAtt);

    //---------------------------------LICORNE_Chamber Al (chamber and cell)
    G4VisAttributes* CoverVisAtt= new G4VisAttributes(G4Colour(0.66,0.66,0.66));   // gray
    CoverVisAtt->SetForceSolid(true);
    CoverVisAtt->SetForceAuxEdgeVisible(false);
    //CoverVisAtt ->SetForceWireframe(true);
    cell_log->SetVisAttributes(CoverVisAtt);
    chamber_log->SetVisAttributes(CoverVisAtt);

    //------------------------------------------LICORNE_Chamber copper shielding
    G4VisAttributes* LICORNE_ChamberShieldingVisAtt= new G4VisAttributes(G4Colour(0.72,0.40,0.23));  // Copper
    LICORNE_ChamberShieldingVisAtt -> SetForceSolid(true);
    LICORNE_ChamberShieldingVisAtt->SetForceAuxEdgeVisible(false);
    shielding_log -> SetVisAttributes(LICORNE_ChamberShieldingVisAtt);


}
