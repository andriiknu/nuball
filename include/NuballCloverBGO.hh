#ifndef NuballCloverBGO_H
#define NuballCloverBGO_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include<vector>

class G4UnionSolid;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;


class CloverBGO {

public:
    CloverBGO();
    ~CloverBGO();

public:
    void SetPosition( G4ThreeVector );
    void SetRotation( G4RotationMatrix );
    void Placement(G4int, G4VPhysicalVolume*, G4bool, G4ThreeVector, G4RotationMatrix, int);
    void SetAbsorber( G4bool );
    void SetHeavyMet( G4bool );

    std::vector<G4LogicalVolume*> GetScoringVolumes() const { return fScoringVolume; }
    G4LogicalVolume* GetScoringVolume(int i) const { return fScoringVolume.at(i); }

private:

    G4double BGOHevimetThickness;
    G4double BGOHevimetDist;
    G4double BGOShieldDist;
    G4double BGOCrystalDist;
    G4double BGOCrystalH;
    G4double BGOCrystalAngle;
    G4double BGOCrystalBoxD;
    G4double BGOCrystalBoxD2;
    G4double BGOCrystalTrapL;
    G4double BGOCrystalW;
    G4double BGOCrystalW2;
    G4double BGOCrystalW3;
    G4double abs1_thickness;
    G4double abs2_thickness;
    G4double AbsorberWidth;

    G4bool   useAbsorber;
    G4bool   useCollimator;

    G4UnionSolid*       solidBGOCrystal; // 1 side
    G4VSolid*           solidBGOCrystal2; // upper crystal, sliced
    G4VSolid*           solidBGOCrystal3; // lower crystal, sliced (mirrored)
    G4VSolid*           solidBGOHevimet;
    G4VSolid*           solidBGOShield;
    G4VSolid*           CloverAbsorber1;
    G4VSolid*           CloverAbsorber2;


private:
    void CreateSolids();

protected:
    std::vector<G4LogicalVolume*>  fScoringVolume;
};

#endif
