#ifndef OUPS_Chamber_H
#define OUPS_Chamber_H 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

//class G4Box;
//class G4Tubs;
//class G4Polycone;
class G4UnionSolid;
//class G4SubtractionSolid;
//class G4IntersectionSolid;
//class G4Polyhedra;
class G4LogicalVolume;
class G4VPhysicalVolume;
class MyDetectorMessenger;
class MyMaterials;


class OUPS_Chamber {
    
public:
    OUPS_Chamber();
    ~OUPS_Chamber();
    
public:
    void SetPosition( G4ThreeVector );
    void SetRotation( G4RotationMatrix );
    void Placement(G4int, G4VPhysicalVolume*, G4bool);
  
private:
    //General materials....
    MyMaterials*   fMat;
    
    G4Material* plastic;
    G4Material* matCapsule;
    
    // Placement Parameters
    G4ThreeVector        position;
    G4RotationMatrix     rotation;
    
    // Solid parameters
    G4int cp1;
    G4double *z;
    G4double *router;
    G4double *rinner;
    G4int cp2;
    G4double *z2;
    G4double *router2;
    G4double *rinner2;
    G4int cp3;
    G4double *z3;
    G4double *router3;
    G4double *rinner3;
   
    // Shapes
    G4Polycone *solchamber1;
    G4Tubs *subtub1;
    G4SubtractionSolid *sub1chamber;
    G4Polycone *solchamber2;
    G4Tubs *subtub2;
    G4SubtractionSolid *sub2chamber;
    G4UnionSolid *union1chamber;
    G4Polycone *solplastic;
    
    // Logical Volume
    G4LogicalVolume *logchamber1;
    G4LogicalVolume *logplastic;
    
    // Physical Volume
    G4VPhysicalVolume *physchamber1;
    G4VPhysicalVolume *physplastic;
    
    
private:
    void CreateSolids();
    void MakeMaterials();
    
};

#endif
