//---------------------------------------------------------------------
// Create the solids defining the target
//---------------------------------------------------------------------
#include "NuballThTarget.hh"
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

TargetTh::TargetTh()
{
  //create the solids.....
  CreateSolids();

}

//Destructor
TargetTh::~TargetTh() { }

void TargetTh::SetPosition(G4ThreeVector thisPos)
{
  position = thisPos*mm;
}

void TargetTh::SetRotation(G4RotationMatrix thisRot) { rotation = thisRot; }

//---------------------------------------------------------------------
// Create the solids defining Phase-II Target
//---------------------------------------------------------------------
void  TargetTh::CreateSolids()
{

  //---I-create-the-cylinder---//
  cylinder_radius = 5.5*cm/2;
  cylinder_thickness = 8*cm;

  cylinder_solid =  new G4Tubs("control_volume",
  0,
  cylinder_radius,
  cylinder_thickness/2,
  0,
  360.*deg);

  target1 =  new G4Tubs("target1",
  0,
  1.386*cm,
  1.*mm/2,
  0,
  360.*deg);

  target2 =  new G4Tubs("target2",
  0,
  1.659*cm,
  1.*mm/2,
  0,
  360.*deg);

  target3 =  new G4Tubs("target3",
  0,
  1.932*cm,
  1.*mm/2,
  0,
  360.*deg);

  target4 =  new G4Tubs("target4",
  0,
  2.205*cm,
  1.*mm/2,
  0,
  360.*deg);

  target5 =  new G4Tubs("target5",
  0,
  2.477*cm,
  1.*mm/2,
  0,
  360.*deg);

  target6 =  new G4Tubs("target6",
  0,
  2.75*cm,
  1.*mm/2,
  0,
  360.*deg);

  //---I-unify-the-two-parts-of-the-target---//

  G4Element* Th232 = new G4Element("Thorium232", "232Th", 90, 232*g/mole);
  G4Material* target_mat = new G4Material("Th", 11.72*g/cm3, 1);
  target_mat->AddElement(Th232, 1.);

  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");

  control_log = new G4LogicalVolume(cylinder_solid,
    air,
    "control_log");

    target1_log = new G4LogicalVolume(target1,
      target_mat,
      "target1_log");

      target2_log = new G4LogicalVolume(target2,
        target_mat,
        "target2_log");

        target3_log = new G4LogicalVolume(target3,
          target_mat,
          "target3_log");

          target4_log = new G4LogicalVolume(target4,
            target_mat,
            "target4_log");

            target5_log = new G4LogicalVolume(target5,
              target_mat,
              "target5_log");

              target6_log = new G4LogicalVolume(target6,
                target_mat,
                "target6_log");


                //-----------------------------------------------------------------------//
                //                    Visual Attributes in the visualisation             //
                //-----------------------------------------------------------------------//

                G4VisAttributes* ControlAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));   // gray
                //CoverVisAtt->SetForceSolid(true);
                ControlAtt->SetForceWireframe(true);
                ControlAtt->SetVisibility(false);
                control_log->SetVisAttributes(ControlAtt);

                G4VisAttributes* TargetAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));   // gray
                //CoverVisAtt->SetForceSolid(true);
                TargetAtt->SetForceWireframe(true);
                TargetAtt->SetVisibility(true);
                target1_log->SetVisAttributes(TargetAtt);
                target2_log->SetVisAttributes(TargetAtt);
                target3_log->SetVisAttributes(TargetAtt);
                target4_log->SetVisAttributes(TargetAtt);
                target5_log->SetVisAttributes(TargetAtt);
                target6_log->SetVisAttributes(TargetAtt);
              }


              //------------------------------------------------------------------
              void TargetTh::Placement(G4int copyNo, G4VPhysicalVolume* physWorld, G4bool checkOverlaps)
              {
                auto physControlTarget = new G4PVPlacement(&rotation,
                  G4ThreeVector(position.x(),position.y(),position.z()),
                  "Control_target",  //its name
                  control_log,     //its logical volume
                  physWorld,      //its mother
                  true,           //no boolean operat
                  copyNo,         //copy number
                  0);

                  new G4PVPlacement(&rotation,
                    G4ThreeVector(0, 0, -3.5*cm),
                    "Target_1",  //its name
                    target1_log,     //its logical volume
                    physControlTarget,      //its mother
                    true,           //no boolean operat
                    copyNo,         //copy number
                    0);

                    new G4PVPlacement(&rotation,
                      G4ThreeVector(0, 0, -2.5*cm),
                      "Target_2",  //its name
                      target2_log,     //its logical volume
                      physControlTarget,      //its mother
                      true,           //no boolean operat
                      copyNo,         //copy number
                      0);

                      new G4PVPlacement(&rotation,
                        G4ThreeVector(0, 0, -1.5*cm),
                        "Target_3",  //its name
                        target3_log,     //its logical volume
                        physControlTarget,      //its mother
                        true,           //no boolean operat
                        copyNo,         //copy number
                        0);

                        new G4PVPlacement(&rotation,
                          G4ThreeVector(0, 0, -0.5*cm),
                          "Target_4",  //its name
                          target4_log,     //its logical volume
                          physControlTarget,      //its mother
                          true,           //no boolean operat
                          copyNo,         //copy number
                          0);

                          new G4PVPlacement(&rotation,
                            G4ThreeVector(0, 0, 0.5*cm),
                            "Target_5",  //its name
                            target5_log,     //its logical volume
                            physControlTarget,      //its mother
                            true,           //no boolean operat
                            copyNo,         //copy number
                            0);


                            new G4PVPlacement(&rotation,
                              G4ThreeVector(0, 0, 1.5*cm),
                              "Target_6",  //its name
                              target6_log,     //its logical volume
                              physControlTarget,      //its mother
                              true,           //no boolean operat
                              copyNo,         //copy number
                              0);

                              new G4PVPlacement(&rotation,
                                G4ThreeVector(0, 0, 2.5*cm),
                                "Target_6",  //its name
                                target6_log,     //its logical volume
                                physControlTarget,      //its mother
                                true,           //no boolean operat
                                copyNo,         //copy number
                                0);

                                new G4PVPlacement(&rotation,
                                  G4ThreeVector(0, 0, 3.5*cm),
                                  "Target_6",  //its name
                                  target6_log,     //its logical volume
                                  physControlTarget,      //its mother
                                  true,           //no boolean operat
                                  copyNo,         //copy number
                                  0);

                                }
