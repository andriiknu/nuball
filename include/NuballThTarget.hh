#ifndef NuballThTarget_H
#define NuballThTarget_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Tubs;
class G4Cons;

class TargetTh
{
public:
  TargetTh();
  ~TargetTh();

  void SetPosition( G4ThreeVector );
  void SetRotation( G4RotationMatrix );
  void Placement(G4int, G4VPhysicalVolume*, G4bool);

private:

  // Placement Parameters
  G4ThreeVector        position;
  G4RotationMatrix     rotation;

  // Solid parameters

  //---Cone---//
  G4double con_radius_a;
  G4double con_radius_b;
  G4double con_thickness;

  G4Cons* con_solid;

  //---Cylinder---//
  G4double cylinder_radius;
  G4double cylinder_thickness;

  G4Tubs* cylinder_solid;

  G4Tubs* target1;
  G4Tubs* target2;
  G4Tubs* target3;
  G4Tubs* target4;
  G4Tubs* target5;
  G4Tubs* target6;

  // Logical Volume
  G4LogicalVolume* control_log;

  G4LogicalVolume* target1_log;
  G4LogicalVolume* target2_log;
  G4LogicalVolume* target3_log;
  G4LogicalVolume* target4_log;
  G4LogicalVolume* target5_log;
  G4LogicalVolume* target6_log;

private:
  void CreateSolids();
};

#endif
