#ifndef NuballClover_H
#define NuballClover_H 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include<vector>

class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

class Clover {

public:
    Clover();
    ~Clover();

public:
    void SetPosition( G4ThreeVector );
    void SetRotation( G4RotationMatrix );
    void Placement(G4int, G4VPhysicalVolume*, G4bool, G4ThreeVector, G4RotationMatrix, int);

    inline G4double GetTaperedEndCapL() {return fEndCapTaperL/mm;};

    std::vector<G4LogicalVolume*> GetScoringVolumes() const { return fScoringVolume; }
    G4LogicalVolume* GetScoringVolume(int i) const { return fScoringVolume.at(i); }

private:

    G4double             lappingSize;

    G4double             fCrystalR;
    G4double             fTotalGeL;
    G4double             fHoleR;
    G4double             fContactThick;
    G4double             fPassiveThick;
    G4double             fEndCapTaperL;
    G4double             fEndCapThickness;
    G4double             fEndCap2Ge;
    G4double             fFudge;
    G4double             fVacuumPosZ;
    G4double             fContact_dZ;
    G4double             fGeLeafPosZ;
    G4double             fGapBetweenLeaves;
    G4double             fGeLeaf_dX;
    G4double             fHole_dX;
    G4double             fHole_dY;

    G4UnionSolid*        solidEndCap;
    G4UnionSolid*        solidVacuum;
    G4UnionSolid*        solidGeLeaf;
    G4UnionSolid*        solidPassivated;
    G4UnionSolid*        solidContact;
    G4UnionSolid*        solidBoreHole;

private:
    void CreateSolids();

protected:
    std::vector<G4LogicalVolume*>  fScoringVolume;
};

#endif
