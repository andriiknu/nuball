//---------------------------------------------------------------------
// Create the solids defining an Eurogam Clover ("Phase-II") detector
//---------------------------------------------------------------------
#include "NuballClover.hh"

//Materials
#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include<string>

using namespace std;

Clover::Clover()
{

  lappingSize = 3.0;

  fTotalGeL     = 70.00 * mm;  //playing around to fudge DC
  fCrystalR     = (25-lappingSize)* mm;  //playing around to fudge DC
  fEndCap2Ge    = 20.00 * mm;  //Distance from Al outer face to Ge

  //added to fudge PTG's efficiency measurements for STUK (Karl)
  fFudge = 5.0*mm; // Fudge factor set to 0, 4.9.2013 JK
  fEndCap2Ge += fFudge;
  fHoleR            =  5.5 * mm; //was 5.0
  fPassiveThick     =  0.5 * mm;
  fContactThick     =  0.5 * mm;
  fGapBetweenLeaves =  0.6 * mm;

  //create the solids.....
  CreateSolids();

}

//Destructor
Clover::~Clover() { }

//---------------------------------------------------------------------
// Create the solids defining Phase-II Clovers
//---------------------------------------------------------------------
void  Clover::CreateSolids()
{


  //An approximate CloverII
  //G4cout << " ----> Constructing archetypal Clover" << G4endl;

  //---------------------------------------------------------
  //end-cap
  G4double endCapFrontThickness = 1.5*mm;//1.5 in the Duchene paper //1.2*mm; in another simulation
  G4double endCapTaperThickness = 1.5*mm;
  G4double endCapSideThickness  = 1.5*mm;

  G4double GeGap      =  fEndCap2Ge; //outsideface of endcap to Ge
  G4double taperAngle =  7.0*degree;

  G4double endCapTotalL = fTotalGeL + GeGap + endCapFrontThickness + 5.*mm; //+ Gap at rear end
  G4double endCapFrontD = 43.5*mm;
  G4double endCapBackD  = 50.5*mm;
  G4double endCapTaperL = 55.0*mm;

  G4double endCapBoxL   = endCapTotalL - endCapTaperL;

  //the tapered part
  G4Trap* solidTaperedCloverEC
    = new G4Trap("taperedCloverEC",
		 endCapTaperL/2., //Half z-length [pDz]
		 0.00*degree,     //Polar angle of line joining centres of the faces @ -/+pDz
		 2.*taperAngle,   //equivalent azimuthal angle
		 endCapFrontD,    //pDy1 half y length at -pDz
		 endCapFrontD,    //pDx1 half x length at -pDz, -pDy1
		 endCapFrontD,    //pDx2 half x length at -pDz, +pDy1
		 0.00*degree,//pAlpha1 wrt y-axis from the centre of the side (lower endcap)
		 endCapBackD,    //pDy2 half y length at +pDz
		 endCapBackD,    //pDx3 half x length at +pDz, -pDy2
		 endCapBackD,    //pDx4 half x length at +pDz, +pDy2
		 0.00*degree); //pAlpha2 wrt y-axis from the centre of the side (upper endcap)

  //the rectangular part.....
  G4Box*        endCapBox  = new G4Box("CloverendCapBox",endCapBackD,endCapBackD,endCapBoxL/2.);
  G4ThreeVector transECBox(   0.*mm, 0.*mm, endCapTaperL/2.+endCapBoxL/2.);

  //add the two together
  solidEndCap = new G4UnionSolid("CloverBox+Taper",solidTaperedCloverEC,endCapBox,0,transECBox);
  //need the taperL for placement
  fEndCapTaperL = endCapTaperL;


  //---------------------------------------------------------
  //end-cap inner vacuum
  G4double endCapDelta_1 = endCapTaperThickness/cos(taperAngle) - endCapFrontThickness*tan(taperAngle);

  G4double endCapVacTaperL = endCapTaperL - endCapFrontThickness;// - endCapDelta_2;
  G4double endCapVacBoxL   = endCapBoxL   - endCapFrontThickness;
  G4double endCapVacTotalL = endCapVacBoxL + endCapVacTaperL;
  G4double endCapVacFrontD = endCapFrontD - endCapDelta_1;
  G4double endCapVacBackD  = endCapBackD  - endCapSideThickness;

  //position of vacuum wrt end-cap
  fVacuumPosZ = (-endCapTotalL + endCapVacTotalL )/2. + 1.5*endCapFrontThickness;

  //tapered part...
  G4Trap* solidTaperVac
    = new G4Trap("CloverTaperVac",
		 endCapVacTaperL/2.,    //Half z-length [pDz]
		 0.00*degree, //Polar angle of line joining centres of the faces @ -/+pDz
		 14.0*degree,   //aequivalent zimuthal angle
		 endCapVacFrontD,    //pDy1 half y length at -pDz
		 endCapVacFrontD,    //pDx1 half x length at -pDz, -pDy1
		 endCapVacFrontD,    //pDx2 half x length at -pDz, +pDy1
		 0.00*degree,//pAlpha1 wrt y-axis from the centre of the side (lower endcap)
		 endCapVacBackD,    //pDy2 half y length at +pDz
		 endCapVacBackD,    //pDx3 half x length at +pDz, -pDy2
		 endCapVacBackD,    //pDx4 half x length at +pDz, +pDy2
		 0.00*degree); //pAlpha2 wrt y-axis from the centre of the side (upper endcap)

  //rectangular part
  G4Box*         endCapVacBox  = new G4Box("CloverendCapVacBox",endCapVacBackD,endCapVacBackD,endCapVacBoxL/2.);
  G4ThreeVector transVacBox(   0.*mm, 0.*mm, (endCapVacTaperL/2.+endCapVacBoxL/2.-0.0001*mm));

  //add them together
  solidVacuum = new G4UnionSolid("CloverVac_Box+Taper",solidTaperVac,endCapVacBox,0,transVacBox);


  //---------------------------------------------------------
  //The Ge crystal...
  G4double GeTaperL    = 36.0*mm;

  G4double smallSquare = (41.0-lappingSize)*mm;
  G4double largeSquare = (45.5-lappingSize)*mm;

  G4double transX = (largeSquare-smallSquare)/2.;
  G4double transY = (largeSquare-smallSquare)/2.;
  fHole_dX = transX + 5.*mm;  //5 mm is a fudge
  fHole_dY = transY + 5.*mm;  //5 mm is a fudge

  //tapered part......
  G4Trap* solidTaper
    = new G4Trap("CloverTaper",
		 GeTaperL/2.,    //Half ? z-length [pDz]
		 5.05*degree,   //Polar angle of line joining centres of the faces @ -/+pDz
		 45.*degree,   //equivalent azimuthal angle  //DOES NOT MAKE SENSE !!
		 smallSquare/2., //pDy1 half y length at -pDz
		 smallSquare/2., //pDx1 half x length at -pDz, -pDy1
		 smallSquare/2., //pDx2 half x length at -pDz, +pDy1
		 0.00*degree,//pAlpha1 wrt y-axis from the centre of the side (lower endcap)
		 largeSquare/2.,    //pDy2 half y length at +pDz
		 largeSquare/2.,    //pDx3 half x length at +pDz, -pDy2
		 largeSquare/2.,    //pDx4 half x length at +pDz, +pDy2
		 0.0*degree); //pAlpha2 wrt y-axis from the centre of the side (upper endcap)

  const G4int numZPlanesGe=4;      // no. polycone planes

  G4double zPlaneGe[numZPlanesGe]; // positions of planes
  zPlaneGe[0] =  0.00*mm;
  zPlaneGe[1] =  2.46*mm;
  zPlaneGe[2] =  5.00*mm;
  zPlaneGe[3] = GeTaperL;

  G4double rInnerGe[numZPlanesGe]; // interior radii
  rInnerGe[0] = rInnerGe[1] = rInnerGe[2] = rInnerGe[3] = 0.0*mm;
  G4double rOuterGe[numZPlanesGe]; // exterior radii
  //rOuterGe[0] = 20.5*mm;  rOuterGe[1] = 23.54*mm;
  rOuterGe[0] = (20.5/25.)*fCrystalR*mm;  rOuterGe[1] = (23.54/25.)*fCrystalR*mm;
  rOuterGe[2] = rOuterGe[3] = fCrystalR;


  G4Polycone* solidCone = new G4Polycone("CloverCone", 0.0*degree, 360.0*degree,
					 numZPlanesGe,
					 zPlaneGe,
					 rInnerGe,
					 rOuterGe);

  G4ThreeVector  transGeCone( -transX/2., -transY/2., -GeTaperL/2.);
  G4IntersectionSolid* taperedCone = new G4IntersectionSolid("CloverTaper+Cone",solidTaper,solidCone,0,transGeCone);

  //back part....
  G4double geBoxL = fTotalGeL - GeTaperL;

  G4Box*    GeBox = new G4Box("CloverGeBox",largeSquare/2.,largeSquare/2.,geBoxL/2.);
  G4Tubs*   GeCyl = new G4Tubs("CloverGeCyl",0.0*mm,fCrystalR,geBoxL/2.,0.*degree,360.*degree);

  G4ThreeVector transGeBox( transX, transY, 0.0*mm);
  G4IntersectionSolid* backPart = new G4IntersectionSolid("CloverBox+Cyl",GeCyl,GeBox,0,transGeBox);

  //add front and back
  G4ThreeVector transBack( -transX/2., -transY/2., (GeTaperL/2.+geBoxL/2.));
  solidGeLeaf = new G4UnionSolid("Clvoergermanium",taperedCone,backPart,0,transBack);

  //z-position of Ge-leaf wrt vacuum
  fGeLeafPosZ = -endCapVacTaperL/2. + GeTaperL/2. + GeGap - endCapFrontThickness;

  //------------------------------------------------------------------
  // Inner bore hole + lithium contact + passivated Ge
  G4double GeDepth      = 15.00 * mm;  //Hole dirilled to this far from face
  G4double passiveThick = fPassiveThick;  //passivated Ge
  G4double contactThick = fContactThick;  //Li contact

  //G4double innerRHole =  0.00*mm;
  G4double holeR      = fHoleR;
  G4double contactR   = holeR + contactThick;
  G4double passiveR   = contactR + passiveThick;
  G4double holeL      = fTotalGeL - GeDepth;
  G4double tubeL      = holeL - holeR;

  //the same translation works for all the following rounded tubes
  G4ThreeVector transSphere(0.01*mm, 0.01*mm, -tubeL/2.-0.1*mm); //if offsets are 0. it does not display !!

  //now add a passivated layer
  G4Sphere* passivatedSphere = new G4Sphere("CloverpassSphere", 0.0*mm, passiveR,           0.*degree, 360.*degree, 0.*degree, 180.*degree);
  G4Tubs*   passivatedTube   = new G4Tubs(  "CloverpassTube",   0.0*mm, passiveR, tubeL/2., 0.*degree, 360.*degree);
  solidPassivated    = new G4UnionSolid("CloverpassivatedGe",passivatedTube,passivatedSphere,0,transSphere);

  //and the Li contact
  G4Sphere* contactSphere  = new G4Sphere("Cloversphere1", 0.0*mm, contactR,           0.*deg, 360.*deg, 0.*deg, 180.*deg);
  G4Tubs*   contactTube    = new G4Tubs(  "Clovertube1",   0.0*mm, contactR, tubeL/2., 0.*deg, 360.*deg);
  solidContact = new G4UnionSolid("CloverliContact",contactTube,contactSphere,0,transSphere);

  //bore out a hole
  G4Sphere* boreSphere  = new G4Sphere("CloverboreSphere", 0.0*mm, holeR,           0.*degree, 360.*degree, 0.*degree, 180.*degree);
  G4Tubs*   boreTube    = new G4Tubs("CloverboreTube",   0.0*mm, holeR, tubeL/2., 0.*degree, 360.*degree);
  solidBoreHole = new G4UnionSolid("CloverboreHole",boreTube,boreSphere,0,transSphere);

  //save this for placement
  fContact_dZ = holeL/2. - contactThick;// - passiveThick;

  //put corners @ (0,0)
  fGeLeaf_dX = largeSquare/2. - transX/2.;


}


//------------------------------------------------------------------
void Clover::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps, G4ThreeVector position, G4RotationMatrix rotation, int i)
{

  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();

  //assign default materials.....
  G4Material* geMaterial = nistManager->FindOrBuildMaterial("G4_Ge");

  G4Material* endCapMaterial = nistManager->FindOrBuildMaterial("DurAl");

  G4Material* vacuumMaterial = nistManager->FindOrBuildMaterial("LabVacuum");

  G4Material* contactMaterial = nistManager->FindOrBuildMaterial("G4_Li");

  //Creation of the logical volumes

  G4LogicalVolume* logicEndCap = new G4LogicalVolume(solidEndCap, endCapMaterial, ("clover_EC"+to_string(i)).c_str(),   0, 0, 0);
  G4LogicalVolume* logicVacuum = new G4LogicalVolume(solidVacuum, vacuumMaterial, ("clover_Vac"+to_string(i)).c_str(),  0, 0, 0);

  G4LogicalVolume* logicGeLeaf     = new G4LogicalVolume(solidGeLeaf, geMaterial, ("CloverLeaf"+to_string(i)).c_str(), 0, 0, 0);
  G4LogicalVolume* logicPassivated = new G4LogicalVolume(solidPassivated, geMaterial, ("CloverpassivatedGe"+to_string(i)).c_str(),  0, 0, 0);
  G4LogicalVolume* logicContact    = new G4LogicalVolume(solidContact, contactMaterial, ("Cloverinner_contact"+to_string(i)).c_str(), 0, 0, 0);
  G4LogicalVolume* logicBoreHole   = new G4LogicalVolume(solidBoreHole, vacuumMaterial, ("Cloverbore_hole"+to_string(i)).c_str(), 0, 0, 0);

  //Creation of the visualisation attributes
  G4VisAttributes* visAttAlCap = new G4VisAttributes( G4Colour(0.9,0.9,0.9,1.0) );
  visAttAlCap->SetVisibility(true);
  visAttAlCap->SetForceWireframe(true);

  G4VisAttributes* visAttGeVac = new G4VisAttributes( G4Colour(0.9,1.0,0.9) );
  visAttGeVac->SetForceWireframe(true);
  visAttGeVac->SetVisibility(false);

  G4VisAttributes* visAttActive = new G4VisAttributes( G4Colour(1.0,1.0,0.0) );
  //visAttActive->SetForceWireframe(true);
  visAttActive->SetVisibility(true);
  G4VisAttributes* visAttActive0 = new G4VisAttributes( G4Colour(1.0,0.2,0.2) );
  //visAttActive->SetForceWireframe(true);
  visAttActive0->SetVisibility(true);
  G4VisAttributes* visAttActive1 = new G4VisAttributes( G4Colour(0.2,1.0,0.2) );
  //visAttActive->SetForceWireframe(true);
  visAttActive1->SetVisibility(true);
  G4VisAttributes* visAttActive2 = new G4VisAttributes( G4Colour(0.2,0.2,1.0) );
  //visAttActive->SetForceWireframe(true);
  visAttActive2->SetVisibility(true);
  G4VisAttributes* visAttActive3 = new G4VisAttributes( G4Colour(0.1,0.1,0.1) );
  //visAttActive->SetForceWireframe(true);
  visAttActive3->SetVisibility(true);

  G4VisAttributes* visAttPassive = new G4VisAttributes(G4Colour(0.0,1.0,1.0) );
  visAttPassive->SetForceWireframe(true);
  visAttPassive->SetVisibility(false);

  G4VisAttributes* visAttLiContact = new G4VisAttributes(G4Colour(1.0,0.0,1.0) );
  //visAttLiContact->SetVisibility(false);
  visAttLiContact->SetVisibility(true);

  G4VisAttributes* visAttHole = new G4VisAttributes( G4Colour(1.0,0.0,1.0) );
  //visAttHole->SetVisibility(false);
  visAttHole->SetVisibility(true);

  logicEndCap->SetVisAttributes(visAttAlCap);
  logicVacuum->SetVisAttributes(visAttGeVac);
  logicGeLeaf->SetVisAttributes(visAttActive0);
  logicPassivated->SetVisAttributes(visAttPassive);
  logicContact->SetVisAttributes(visAttLiContact);
  logicBoreHole->SetVisAttributes(visAttHole);


  // Defining Scoring Volumes
  fScoringVolume.push_back(logicGeLeaf);

  //=================================================================================
  //Do not know why, but the positioning seems to be with respect to the Taper-part :
  //setting the z-position as endCapTaperL/2 puts the front face at z = 0 mm
  //=================================================================================
  G4double vacuum_PosZ = fVacuumPosZ;
  G4double geLeaf_PosZ = fGeLeafPosZ;


  //Physical placement of these solids......
  //physiEndCap = new G4PVPlacement(rotMat,
  //					 position,

  auto physiEndCap = new G4PVPlacement(G4Transform3D(rotation, position),
					 "Clover_EC",       //its name
					 logicEndCap,//its logical volume
					 physiMother,         //its mother
					 true,               //no boolean operat
					 copyNo,             //copy number
					 checkOverlaps);              //overlap check

  auto physiVacuum = new G4PVPlacement(0,                   //rotation
					 G4ThreeVector(0.*mm,0.*mm,vacuum_PosZ),
					 "Clover_Vac",       //its name
					 logicVacuum, //its logical volume
					 physiEndCap, //its mother
					 true,                //no boolean operat
					 copyNo,              //copy number
					 checkOverlaps);               //overlap check

  //Now for the placement of the leaves in each clover......
  G4RotationMatrix* rmC;
  G4double dPos = fGeLeaf_dX + fGapBetweenLeaves/2.;
  G4double leafX;
  G4double leafY;
  G4double leafZ;

  for(G4int l = 0; l < 4; l++)
  {
    //the rotation
    rmC = new G4RotationMatrix;
    rmC->set(0,0,0);
    rmC->rotateZ(90.*degree*(4-l));
    rmC->invert();
    //the x-translation
    if(l < 2) {
      leafX = dPos;
    } else {
      leafX = -dPos;
    }
    //the y-translation
    if(l == 0 || l == 3 ) {
      leafY = dPos;
    } else {
      leafY = -dPos;
    }
    //the z-translation
    leafZ = geLeaf_PosZ;

    G4ThreeVector leafdist;
    leafdist.operator=(position);
    leafdist.operator+=(G4ThreeVector(0.*mm,0.*mm,-vacuum_PosZ));
    leafdist.operator+=(G4ThreeVector(0.*mm,0.*mm,-leafZ));

    auto physiGeLeaf = new G4PVPlacement(rmC,                       //rotation
					      G4ThreeVector(leafX,leafY,leafZ),
					      ("CloverGermanium"+to_string(l)).c_str(),           //its name
					      logicGeLeaf,     //its logical volume
					      physiVacuum,        //its mother
					      true,               //no boolean operat
					      copyNo+l,           //copy number
					      checkOverlaps);     //overlap check

    auto physiPassivated = new G4PVPlacement(0,                   //rotation
						  G4ThreeVector(-fHole_dX, -fHole_dY, fContact_dZ),
						  ("CloverGePassivated"+to_string(l)).c_str(),
						  logicPassivated,
						  physiGeLeaf,
						  false,copyNo+l,checkOverlaps);

    auto physiContact = new G4PVPlacement(0,                   //rotation
					       G4ThreeVector(0.*mm,0.*mm, 0.0*mm),//-fContact_dZ),
					      ("CloverLiContact"+to_string(l)).c_str(),
					       logicContact,
					       physiPassivated,
					       false,copyNo+l,checkOverlaps);


    new G4PVPlacement(0,                   //rotation
						G4ThreeVector(0.*mm,0.*mm, 0.0*mm),//-fContact_dZ),
						("CloverBoreHole"+to_string(l)).c_str(),
						logicBoreHole,
						physiContact,
						false,copyNo+l,checkOverlaps);
  }
}
