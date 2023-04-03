#ifndef LICORNE_Chamber_H
#define LICORNE_Chamber_H 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;


class LICORNE_Chamber {

public:
    LICORNE_Chamber();
    ~LICORNE_Chamber();

public:
    void SetPosition( G4ThreeVector );
    void SetRotation( G4RotationMatrix );
    void Placement(G4int, G4VPhysicalVolume*, G4bool);

public:
    inline G4double GetCellThickness(){return gascell_length;};
    inline G4double GetChamberThickness(){return chamber_length;};
    inline G4double GetChamberRadius(){return chamber_diameter/2.;};

private:

    G4Material* ChamberMaterial;
    G4Material* innerChamberMaterial;
    G4Material* cellMaterial;
    G4Material* innercellMaterial;
    G4Material* ControlMaterial;
    G4Material* innershieldingMaterial;
    G4Material* shieldingMaterial;

    // Placement Parameters
    G4ThreeVector        position;
    G4RotationMatrix     rotation;

    // Solid parameters
    G4double shielding_length;
    G4double shielding_diameter;
    G4double shielding_innerdiameter;
    G4double chamber_diameter;
    G4double chamber_length;
    G4double gascell_diameter;
    G4double gascell_length;
    G4double control_diameter;
    G4double control_length;

    // Shapes
    G4Tubs* controlchamber_solid;
    G4Tubs* chamber_solid;
    G4Tubs* innerchamber_solid;
    G4Tubs* cell_solid;
    G4Tubs* innercell_solid;
    G4Tubs* shielding_solid;
    G4Tubs* innershielding_solid;

    // Logical Volume
    G4LogicalVolume* controlchamber_log;
    G4LogicalVolume* chamber_log;
    G4LogicalVolume* innerchamber_log;
    G4LogicalVolume* cell_log;
    G4LogicalVolume* innercell_log;
    G4LogicalVolume* shielding_log;
    G4LogicalVolume* innershielding_log;

    // Physical Volume
    G4VPhysicalVolume* controlchamber_phys;
    G4VPhysicalVolume* chamber_phys;
    G4VPhysicalVolume* innerchamber_phys;
    G4VPhysicalVolume* cell_phys;
    G4VPhysicalVolume* innercell_phys;
    G4VPhysicalVolume* shielding_phys;
    G4VPhysicalVolume* innershielding_phys;

private:
    void CreateSolids();

};

#endif
