//---------------------------------------------------------------------
// Create the solids defining an Ortec Phase-I detector
//---------------------------------------------------------------------
#include "NuballPhaseI.hh"

//Materials
#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "globals.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <iostream>

PhaseI::PhaseI()
{

    fAlCap2Ge_PhaseI = 15.1 * mm; //fAlCap2Ge_Ortec;  //Distance from Al outer face to Ge
    fGeTaperL_PhaseI = 30.0 * mm; //fGeTaperL_Ortec;
    fTotalGeL_PhaseI = 70.0 * mm; //fTotalGeL_Ortec;
    fGeOuterD_PhaseI = 66.0 * mm; //fGeOuterD_Ortec;
    fGeInnerD_PhaseI = 64.0 * mm; //fGeInnerD_Ortec;   //dR = 10./100.*taperL

    fEndCapFrontR_PhaseI    =  34.0*mm;
    fEndCapBackR_PhaseI     =  44.0*mm;
    fEndCapTaperL_PhaseI    = 100.0*mm; //=> tan(theta) = 10/100
    fEndCapTubeL_PhaseI     = 260.0*mm;
    fEndCapThickness_PhaseI =   1.2*mm;
    CreatePhaseISolids();
}


//------------------------------------------------------------------
//Destructor
PhaseI::~PhaseI() {
    //
}
/*
void PhaseI::SetPosition(G4ThreeVector thisPos) {

    position = thisPos*mm;
    //G4cout << " ----> A Phase-I will be placed at " << position/mm << " mm" << G4endl;

}

void PhaseI::SetRotation(G4RotationMatrix thisRot) {

    rotation = thisRot;

    //G4cout << " ----> Single Ge will be placed at " << position/mm << " mm" << G4endl;

}
*/

//---------------------------------------------------------------------
//---------------------------------------------------------------------
// Create the solids defining Phase-I detectors.......
//---------------------------------------------------------------------
//---------------------------------------------------------------------
void PhaseI::CreatePhaseISolids()
{
    //------------------------------------------
    // If it has been set up correctly then the
    // different dimesions for the Ortec, Eurisys
    // and Shorty have already been accounted for
    //-------------------------------------------

    G4double phiStartGe =   0.0*deg; // full circle...
    G4double phiTotalGe = 360.0*deg;

    G4double endCapFrontR    = fEndCapFrontR_PhaseI;
    G4double endCapBackR     = fEndCapBackR_PhaseI;
    G4double endCapTaperL    = fEndCapTaperL_PhaseI;
    G4double endCapCylinderL = fEndCapTubeL_PhaseI;

    //------------------------------
    // Al container tube
    const G4int numZPlanesAl=3;      // no. polycone planes
    G4double zPlaneAl[numZPlanesAl]; // positions of planes
    zPlaneAl[0] = 0.0*mm;
    zPlaneAl[1] = endCapTaperL;
    zPlaneAl[2] = endCapTaperL + endCapCylinderL;
    G4double rInnerAl[numZPlanesAl]; // interior radii
    rInnerAl[0] = 0.0*mm;
    rInnerAl[1] = 0.0*mm;
    rInnerAl[2] = 0.0*mm; // I think this is supposed to be 0, to avoid geometry warnings (JK) //was endCapBackR
    G4double rOuterAl[numZPlanesAl]; // exterior radii
    rOuterAl[0] = endCapFrontR;
    rOuterAl[1] = endCapBackR;  //=> tan(theta) = 10/100
    rOuterAl[2] = endCapBackR;

    solidAlCap_PhaseI = new G4Polycone("PhaseIcapAl", phiStartGe, phiTotalGe, numZPlanesAl, zPlaneAl, rInnerAl, rOuterAl);


    //------------------------------
    // Vacuum tube inside container : make Al walls 1.2 mm
    G4double AlThick = fEndCapThickness_PhaseI;

    const G4int numZPlanesVac=3;       // no. polycone planes
    G4double zPlaneVac[numZPlanesVac]; // positions of planes
    zPlaneVac[0] = zPlaneAl[0] + AlThick;                             //AlThick
    zPlaneVac[1] = zPlaneAl[1] + 10.*AlThick*(1./cos(atan(0.1)) - 1.);// AlThick;         //
    zPlaneVac[2] = zPlaneAl[2] - AlThick;
    G4double rInnerVac[numZPlanesVac]; // interior radii
    rInnerVac[0] = 0.0*cm;
    rInnerVac[1] = 0.0*cm;
    rInnerVac[2] = 0.0*cm;
    G4double rOuterVac[numZPlanesVac]; // exterior radii
    rOuterVac[0] = rOuterAl[0] + 0.1*AlThick - AlThick/cos(atan(0.1));
    rOuterVac[1] = rOuterAl[1] - AlThick;
    rOuterVac[2] = rOuterAl[2] - AlThick;

    solidVacuum_PhaseI = new G4Polycone("PhaseIvacuum", phiStartGe, phiTotalGe, numZPlanesVac, zPlaneVac, rInnerVac, rOuterVac);

    //could consider subtracting the above vacuum from the AlCap => one less solid to place

    //-------------------------------------------------------------------------------------
    // Germanium.....
    //G4double Gap     =  3.5 * mm;
    G4double delta   = 2.0 * mm;  //rounded edges

    const G4int numZPlanesGe=4;      // no. polycone planes
    G4double zPlaneGe[numZPlanesGe]; // positions of planes
    zPlaneGe[0] = 0.0 * mm;
    zPlaneGe[1] = delta;
    zPlaneGe[2] = fGeTaperL_PhaseI;
    zPlaneGe[3] = fTotalGeL_PhaseI;
    G4double rInnerGe[numZPlanesGe]; // interior radii
    rInnerGe[0] = 0.0*cm;
    rInnerGe[1] = 0.0*cm;
    rInnerGe[2] = 0.0*cm;
    rInnerGe[3] = 0.0*cm;
    G4double rOuterGe[numZPlanesGe]; // exterior radii
    rOuterGe[0] = fGeInnerD_PhaseI/2. - delta;     //a little rounding...
    rOuterGe[1] = fGeInnerD_PhaseI/2. + 0.1*delta; //tan(theta) = 0.1
    rOuterGe[2] = fGeOuterD_PhaseI/2.;
    rOuterGe[3] = fGeOuterD_PhaseI/2.;


    solidGe_PhaseI = new G4Polycone("PhaseIgermanium", phiStartGe, phiTotalGe, numZPlanesGe, zPlaneGe, rInnerGe, rOuterGe);

    //------------------------------------------------------------------
    // Inner bore hole + lithium contact + passivated Ge
    G4double GeDepth      =  5.2 * mm;           //treated => could make hole 0.2 mm longer than nominal 5.2 mm
    G4double contactThick =  0.50*mm;
    G4double passiveThick =  0.30*mm;

    //G4double zPasGeLength = fTotalGeL_PhaseI - GeDepth; // Passivated Ge : inner contact 0.3 mm thick

    //G4double innerRHole =  0.00*mm;
    //G4double holeR      =  5.00*mm;            //treated => could make it 0.5 mm bigger than the nominal 5.0*mm
    G4double holeR      =  5.50*mm;            //treated => could make it 0.5 mm bigger than the nominal 5.0*mm fudging the code DC
    G4double contactR   = holeR + contactThick;
    G4double passiveR   = contactR + passiveThick;
    G4double holeL      = fTotalGeL_PhaseI - GeDepth;
    G4double tubeL      = holeL - holeR;

    //save this for placement
    fContact_dZ_PhaseI = fTotalGeL_PhaseI - tubeL/2.;

    //the same translation works for all the following.
    G4ThreeVector transHole( 0.001*mm, 0.001*mm, -tubeL/2.-0.001*mm); //if offsets are 0. it does not display !!

    //now add a passivated layer
    G4Sphere* passivatedSphere = new G4Sphere("PhaseIpassSphere", 0.0*mm, passiveR,           0.*degree, 360.*degree, 0.*degree, 180.*degree);
    G4Tubs*   passivatedTube   = new G4Tubs(  "PhaseIpassTube",   0.0*mm, passiveR, tubeL/2., 0.*degree, 360.*degree);
    solidPassivated_PhaseI = new G4UnionSolid("PhaseIpassivatedGe",passivatedTube,passivatedSphere,0,transHole);

    //and the Li contact
    G4Sphere* contactSphere = new G4Sphere("PhaseIsphere1", 0.0*mm, contactR,           0.*deg, 360.*deg, 0.*deg, 180.*deg);
    G4Tubs*   contactTube   = new G4Tubs(  "PhaseItube1",   0.0*mm, contactR, tubeL/2., 0.*deg, 360.*deg);
    solidContact_PhaseI = new G4UnionSolid("PhaseIliContact",contactTube,contactSphere,0,transHole);

    //bore out a hole
    G4Sphere* boreSphere = new G4Sphere("PhaseIboreSphere",0.0*mm, holeR,           0.*degree, 360.*degree, 0.*degree, 180.*degree);
    G4Tubs*   boreTube   = new G4Tubs(  "PhaseIboreTube",  0.0*mm, holeR, tubeL/2., 0.*degree, 360.*degree);
    solidBoreHole_PhaseI = new G4UnionSolid("PhaseIhole+sphere",boreTube,boreSphere,0,transHole);


}


//------------------------------------------------------------------
//------------------------------------------------------------------
void PhaseI::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps, G4ThreeVector position, G4RotationMatrix rotation, int i)
{

  //--------------------------------------------------------------------------------------------------//
  //Definition of the materials

  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* geMaterial = nistManager->FindOrBuildMaterial("G4_Ge");

  G4Material* endCapMaterial = nistManager->FindOrBuildMaterial("DurAl");

  G4Material* vacuumMaterial = nistManager->FindOrBuildMaterial("LabVacuum");

  G4Material* contactMaterial = nistManager->FindOrBuildMaterial("G4_Li");

  //--------------------------------------------------------------------------------------------------//

  //Creation of the logical volumes

  G4LogicalVolume* logicAlCap_PhaseI = new G4LogicalVolume(solidAlCap_PhaseI, endCapMaterial, ("PhaseIcapAl"+std::to_string(i)).c_str(), 0, 0, 0);

  G4LogicalVolume* logicVacuum_PhaseI = new G4LogicalVolume(solidVacuum_PhaseI, vacuumMaterial, ("PhaseIvacuum"+std::to_string(i)).c_str(), 0, 0, 0);

  G4LogicalVolume* logicGe_PhaseI = new G4LogicalVolume(solidGe_PhaseI, geMaterial, ("PhaseIGe"+std::to_string(i)).c_str(), 0, 0, 0);

  G4LogicalVolume* logicPassivated_PhaseI = new G4LogicalVolume(solidPassivated_PhaseI, geMaterial, ("PhaseIpassiveGe"+std::to_string(i)).c_str(), 0, 0, 0);

  G4LogicalVolume* logicContact_PhaseI = new G4LogicalVolume(solidContact_PhaseI, contactMaterial, ("PhaseIlithium"+std::to_string(i)).c_str(), 0, 0, 0);

  G4LogicalVolume* logicBoreHole_PhaseI = new G4LogicalVolume(solidBoreHole_PhaseI, vacuumMaterial, ("PhaseIhole"+std::to_string(i)).c_str(), 0, 0, 0);

  //Creation of the visualisation attribute

  G4VisAttributes* visAttActive = new G4VisAttributes( G4Colour(1.0,1.0,0.0) );
  //visAttActive->SetForceWireframe(true);
  visAttActive->SetVisibility(true);

  G4VisAttributes* visAttHole = new G4VisAttributes( G4Colour(1.0,0.0,1.0) );
  visAttHole->SetVisibility(false);
  //visAttHole->SetVisibility(false);

  G4VisAttributes* visAttGeVac = new G4VisAttributes( G4Colour(0.9,1.0,0.9) );
  visAttGeVac->SetForceWireframe(true);
  //visAttGeVac->SetVisibility(true);
  visAttGeVac->SetVisibility(false);

  G4VisAttributes* visAttPassive = new G4VisAttributes(G4Colour(0.0,1.0,1.0) );
  //visAttPassive->SetVisibility(true);
  visAttPassive->SetVisibility(false);

  G4VisAttributes* visAttAlCap = new G4VisAttributes( G4Colour(0.9,1.0,1.0) );
  visAttAlCap->SetVisibility(true);
  //visAttAlCap->SetForceWireframe(true);

  G4VisAttributes* visAttBGO = new G4VisAttributes(G4Colour(0.5,1.0,1.0) );
  visAttBGO->SetVisibility(false);
  //visAttBGO->SetForceWireframe(true);

  logicGe_PhaseI->SetVisAttributes(visAttActive);
  logicBoreHole_PhaseI->SetVisAttributes(visAttPassive);
  logicPassivated_PhaseI->SetVisAttributes(visAttPassive);
  logicContact_PhaseI->SetVisAttributes(visAttHole);
  logicVacuum_PhaseI->SetVisAttributes(visAttGeVac);
  logicAlCap_PhaseI->SetVisAttributes(visAttAlCap);


  //Definition of the scoring volume
  fScoringVolume.push_back(logicGe_PhaseI);

    //Placement of one copy of a Phase I
    std::cout << "PhaseI " << copyNo << " being placed" << std::endl;
    position*=mm;

    auto physiAlCap_PhaseI = new G4PVPlacement(G4Transform3D(rotation, position),
                                      ("PhaseICapAl"+std::to_string(copyNo)).c_str(),           //its name
                                      logicAlCap_PhaseI, //its logical volume
                                      physiMother,       //its mother
                                      false,             //no boolean operat
                                      copyNo,            //copy number
                                      checkOverlaps);             //overlap check
    //vacuum between AlCap and Ge
    auto physiVacuum_PhaseI = new G4PVPlacement(0,
                                          G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm),
                                           ("PhaseIVacuum"+std::to_string(copyNo)).c_str(),           //its name
                                           logicVacuum_PhaseI, //its logical volume
                                           physiAlCap_PhaseI,  //its mother
                                           false,              //no boolean operat
                                           copyNo,            //copy number
                                           checkOverlaps);             //overlap check

    //the germanium
    G4double geZ = fAlCap2Ge_PhaseI + fEndCapThickness_PhaseI;//-delta
    //G4cout << "Distance of Ge cristal to the center " << (position.getR() - geZ/2.)/cm << " cm" << G4endl;
    //G4cout << "Opening Angle is : " << TMath::ATan((fGeInnerD_PhaseI/2)/(position.getR() - geZ/2.))*TMath::RadToDeg()*2 << G4endl;
    auto physiGe_PhaseI = new G4PVPlacement(0,
                                        G4ThreeVector(0.0*mm, 0.0*mm, geZ),
                                       ("PhaseIGermanium"+std::to_string(copyNo)).c_str(),         //its name
                                       logicGe_PhaseI,        //its logical volume
                                       physiVacuum_PhaseI,    //its mother
                                       false,             //no boolean operat
                                       copyNo,            //copy number
                                       checkOverlaps);             //overlap check

    //passivated region between Ge-dector and Li-contact
    auto physiPassivated_PhaseI = new G4PVPlacement(0,
                                               G4ThreeVector(0.0*mm,0.0*mm,fContact_dZ_PhaseI),
                                               ("PhaseIPassiveGe"+std::to_string(copyNo)).c_str(),    //its name
                                               logicPassivated_PhaseI, //its logical volume
                                               physiGe_PhaseI,        //its mother
                                               false,          //no boolean operat
                                               copyNo,              //copy number
                                               checkOverlaps);          //overlap check

    //the lithium contact
    auto physiContact_PhaseI = new G4PVPlacement(0, G4ThreeVector(0.0*mm,0.0*mm,0.0*mm), //dz = contact-thick ?
                                            ("PhaseILiContact"+std::to_string(copyNo)).c_str(),         //its name
                                            logicContact_PhaseI, //its logical volume
                                            physiPassivated_PhaseI,        //its mother (I think it should be physiGe_PhaseI)
                                            false,          //no boolean operat
                                            copyNo,              //copy number
                                            checkOverlaps);          //overlap check


    //the lithium contact
    new G4PVPlacement(0, G4ThreeVector(0.0*mm,0.0*mm,0.0*mm), //dz = contact-thick ?
                                             ("PhaseIBoreHole"+std::to_string(copyNo)).c_str(),         //its name
                                             logicBoreHole_PhaseI, //its logical volume
                                             physiContact_PhaseI,
                                             false,          //no boolean operat
                                             copyNo,              //copy number
                                             checkOverlaps);          //overlap check
}
