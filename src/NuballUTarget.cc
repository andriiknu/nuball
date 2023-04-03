//---------------------------------------------------------------------
// Create the solids defining the target
//---------------------------------------------------------------------
#include "NuballUTarget.hh"
//to
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

Target::Target()
{
  //create the solids.....
  CreateSolids();

}

//Destructor
Target::~Target() { }

void Target::SetPosition(G4ThreeVector thisPos)
{
  position = thisPos*mm;
}

void Target::SetRotation(G4RotationMatrix thisRot) { rotation = thisRot; }

//---------------------------------------------------------------------
// Create the solids defining Phase-II Target
//---------------------------------------------------------------------
void  Target::CreateSolids()
{

  //---I-create-the-troncated-cons---//
  con_radius_a = 1.25*cm;
  con_radius_b = 3.0*cm;
  con_thickness = 7.5*cm;

  con_solid = new G4Cons("ConeTronque",
    0, con_radius_a, 0, con_radius_b, con_thickness/2,
    0, 360.*deg);

  //---I-create-the-cylinder---//
  cylinder_radius = 3.*cm;
  cylinder_thickness = 3.*cm;

  cylinder_solid =  new G4Tubs("Cylindre",
			       0,
			       cylinder_radius,
			       cylinder_thickness/2,
			       0,
			       360.*deg);

  //---I-unify-the-two-parts-of-the-target---//
  target_solid = new G4UnionSolid("Target_solid", cylinder_solid, con_solid, 0, G4ThreeVector(0, 0, -5.25*cm));

  G4Element* U8 = new G4Element("Uranium238", "U8", 92, 238*g/mole);
  G4Material* target_mat = new G4Material("U", 1.05*g/cm3, 1);
  target_mat->AddElement(U8, 1.);

  target_log = new G4LogicalVolume(target_solid,
				   target_mat,
				   "Target_log");


  //-----------------------------------------------------------------------//
  //                    Visual Attributes in the visualisation             //
  //-----------------------------------------------------------------------//

  G4VisAttributes* CylinderVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));   // gray
  //CoverVisAtt->SetForceSolid(true);
  CylinderVisAtt ->SetForceWireframe(true);
  target_log->SetVisAttributes(CylinderVisAtt);


}


//------------------------------------------------------------------
void Target::Placement(G4int copyNo, G4VPhysicalVolume* physWorld, G4bool checkOverlaps)
{
new G4PVPlacement(&rotation,
				   G4ThreeVector(position.x(),position.y(),position.z()+4.5*cm),
				   "Target",  //its name
				   target_log,     //its logical volume
				   physWorld,      //its mother
				   true,           //no boolean operat
				   copyNo,         //copy number
				   0);
}
