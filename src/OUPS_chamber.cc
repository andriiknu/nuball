//---------------------------------------------------------------------
// Create the solids defining an FATIMA OUPS_Chamber detector
//---------------------------------------------------------------------
#include "OUPS_Chamber.hh"

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

OUPS_Chamber::OUPS_Chamber()
{

    // OUPS chamber parameters
    // First Part
    cp1 = 10;
    z= new G4double[cp1];
    G4double buffz[10] = {-490*mm,
		      -478*mm,
		      -478*mm,
		      -140.82*mm,//   -135.4493118526888*mm,
		      -130*mm,
		      -121.54*mm,
		      -100*mm,
		      -80.51*mm,
		      -58.1464*mm,
		      52*mm};
    for(int i=0; i<cp1;i++) z[i]=buffz[i];
    router= new G4double[cp1];
    G4double buffrouter[10]={50*mm,
        50*mm,
        18.5*mm,
        18.5*mm,
        21.68*mm,
        33.23*mm,
        76.67*mm,
        96.92*mm,
        105*mm,
        105*mm};
    for(int i=0; i<cp1;i++) router[i]=buffrouter[i];
    rinner= new G4double[cp1];
    G4double buffrinner[10] ={15.5*mm,
        15.5*mm,
        15.5*mm,
        15.5*mm,
        18.21*mm,//31.85486009350359*mm ,//<-
        25.97*mm,
        71.61*mm,
        92.98*mm,
        102*mm,
        102*mm};
    for(int i=0; i<cp1;i++) rinner[i]=buffrinner[i];

    // Second Part
    cp2 = 4;
    z2= new G4double[cp2]; G4double buffz2[4] = {0,12*mm,12*mm,200*mm};for(int i=0; i<cp2;i++) z2[i]=buffz2[i];
    router2= new G4double[cp2];G4double buffrouter2[4]={65*mm,65*mm,40.5*mm,40.5*mm};for(int i=0; i<cp2;i++) router2[i]=buffrouter2[i];
    rinner2= new G4double[cp2];G4double buffrinner2[4]={37.5*mm,37.5*mm,37.5*mm,37.5*mm};for(int i=0; i<cp2;i++) rinner2[i]=buffrinner2[i];

    // Plastic Part
    cp3 = 8;
    z3= new G4double[cp3];
    G4double buffz3[8]={52*mm,
        60.6218*mm,
        77.9423*mm,
        105.258*mm,
        126.143*mm,
        185*mm,
        185*mm,
        200*mm};for(int i=0; i<cp3;i++) z3[i]=buffz3[i];
    router3= new G4double[cp3];
    G4double buffrouter3[9] = {105*mm,
        105*mm,
        99*mm,
        69.26*mm,
        58*mm,
        58*mm,
        82.5*mm,
        82.5*mm};for(int i=0; i<cp3;i++) router3[i]=buffrouter3[i];
    rinner3= new G4double[cp3];
    G4double buffrinner3[9] = {97*mm,
        97*mm,
        91*mm,
        61*mm,
        50*mm,
        50*mm,
        50*mm,
        50*mm};for(int i=0; i<cp3;i++) rinner3[i]=buffrinner3[i];



    //---------------------------------
    //make the required materials and assign them
    G4NistManager* nistManager = G4NistManager::Instance();

    //assign default materials.....
    plastic        = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");

    matCapsule = new G4Material("DurAl", 2.8*g/cm3, 4);
    matCapsule->AddElement(nistManager->FindOrBuildElement("G4_Cu"), 3.5*perCent);
    matCapsule->AddElement(nistManager->FindOrBuildElement("G4_Mg"), 0.5*perCent);
    matCapsule->AddElement(nistManager->FindOrBuildElement("G4_Mn"), 0.6*perCent);
    matCapsule->AddElement(nistManager->FindOrBuildElement("G4_Al"), 95.4*perCent);

    //create the solids.....
    CreateSolids();

}

//Destructor
OUPS_Chamber::~OUPS_Chamber() { }

void OUPS_Chamber::SetPosition(G4ThreeVector thisPos) {
    position = thisPos*mm;
    G4cout << " ----> A OUPS_Chamber will be placed at " << position/mm << " mm" << G4endl;
}

void OUPS_Chamber::SetRotation(G4RotationMatrix thisRot) { rotation = thisRot; }

//---------------------------------------------------------------------
// Create the solids defining Phase-II OUPS_Chambers
//---------------------------------------------------------------------
void  OUPS_Chamber::CreateSolids()
{

    //An approximate OUPS_Chamber
    G4cout << " ----> Constructing archetypal OUPS_Chamber volumes" << G4endl;

    //------------ CONTROL VOLUME ------------//
    // Definition of the geometric volume
    // First Part

    //G4cout << "First Part" << G4endl;
    //for(int i=0; i< cp1; i++) G4cout << "z["<< i << "] = " << z[i] << "; rinner[" << i << "] = " << rinner[i] << "; router[" << i << "] = " << router[i] << G4endl;
    solchamber1 = new G4Polycone("chamber1sol",0*deg, 360.0*deg,
                                 cp1,z,rinner,router);

    subtub1 = new G4Tubs("subtub1",0*mm,42.5*mm,50*mm,0,360*deg);
    G4RotationMatrix rm1;
    rm1.rotateX(-90*deg);

    sub1chamber =
    new G4SubtractionSolid("sub1chamber",solchamber1,subtub1,
                           G4Transform3D(rm1,G4ThreeVector(0,-100*mm,0)));

    // Second Part
    //G4cout << "First Part" << G4endl;
    solchamber2 = new G4Polycone("chamber2sol",0*deg, 360.0*deg,
                                 cp2,z2,rinner2,router2);
    G4Tubs *subtub2 = new G4Tubs("subtub2",0*mm,102*mm,52*mm,0,360*deg);
    G4RotationMatrix rm2;
    rm2.rotateX(90*deg);
    G4SubtractionSolid *sub2chamber =
    new G4SubtractionSolid("sub2chamber",solchamber2,subtub2,
                           G4Transform3D(rm2,G4ThreeVector(0,0,137*mm)));
    G4RotationMatrix rm3;
    rm3.rotateX(-90*deg);
    rm3.rotateY(180*deg);

    // Final solid for the chamber
    union1chamber = new G4UnionSolid("union1chamber",sub2chamber,sub1chamber,
                                                   G4Transform3D(rm3,G4ThreeVector(0,0,136.99*mm)));

    // Plastic Part
    //G4cout << "Plastic Part" << G4endl;
    //G4cout << cp3 << G4endl;
    //for(int i=0; i< cp3; i++) G4cout << "z3["<< i << "] = " << z3[i] << "; rinner3[" << i << "] = " << rinner3[i] << "; router3[" << i << "] = " << router3[i] << G4endl;
    solplastic = new G4Polycone("plasticsol",0*deg, 360.0*deg,cp3,
                                            z3,rinner3,router3);

    //G4cout << "DONE" << G4endl;

}


//------------------------------------------------------------------
void OUPS_Chamber::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{
    G4cout << " ----> Constructing archetypal OUPS_Chamber logical and physical volumes" << G4endl;

    // Definition of the rotation matrix that might be used for placement
    G4RotationMatrix* rmC = new G4RotationMatrix(0,0,0);

    //------------ Chamber VOLUME ------------//
    // Definition of the control logical and physical volumes
    logchamber1 = new G4LogicalVolume(union1chamber,matCapsule,"chambre1log",0,0,0);
    G4RotationMatrix rm4;
    rm4.rotateX(90*deg);
    rm4.rotateZ(180*deg);
    physchamber1 = new G4PVPlacement(G4Transform3D(rm4,
                                     G4ThreeVector(0,-137*mm,0)),
                                     "chamber1phys",logchamber1,
                                     physiMother, false, checkOverlaps);

    //------------ Plastic VOLUME ------------//
    // Definition of the control logical and physical volumes
    logplastic  = new G4LogicalVolume(solplastic,plastic,"plasticlog");
    physplastic = new G4PVPlacement(&rotation,
                                    G4ThreeVector(position.x(),position.y(),position.z()),
                                    "chamber2phys",logplastic,
                                    physiMother, false, checkOverlaps);



    //-----------------------------------------------------------------------//
    //                    Visual Attributes in the visualisation             //
    //-----------------------------------------------------------------------//
    //---------------------------------OUPS_Chamber First Part
    G4VisAttributes *vischal = new G4VisAttributes(G4Colour(.6, .6, .6));
    logchamber1->SetVisAttributes(vischal);

    //---------------------------------OUPS_Chamber Plastic
    G4VisAttributes *visplasitc = new G4VisAttributes(G4Colour(1., 1., 1.));
    logplastic->SetVisAttributes(visplasitc);

}
