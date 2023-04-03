#ifndef NuballUTarget_H
#define NuballUTarget_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Tubs;
class G4Cons;

class Target
{    
public:
  Target();
  ~Target();
    
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
  
  //---Target---//
  G4UnionSolid* target_solid;
    
  // Logical Volume
  G4LogicalVolume* target_log;
    
private:
  void CreateSolids();
};

#endif
