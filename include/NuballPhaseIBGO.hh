#ifndef NuballPhaseIBGO_H
#define NuballPhaseIBGO_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include <string>
#include<vector>

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;


class PhaseIBGO {

public:
    PhaseIBGO();
    ~PhaseIBGO();

public:
    //void SetPosition( G4ThreeVector );
    //void SetRotation( G4RotationMatrix );
    void Placement(G4int, G4VPhysicalVolume*, G4bool, G4ThreeVector, G4RotationMatrix, int);
    void SetAbsorber( G4bool );
    void SetHeavyMet( G4bool );

    //inline G4double GetShieldDistance {return BGOShieldDist;};
    std::vector<G4LogicalVolume*> GetScoringVolumes() const { return fScoringVolume; }
    G4LogicalVolume* GetScoringVolume(int i) const {return fScoringVolume.at(i);}

private:
    //General materials....

    G4bool   useAbsorber;
    G4bool   useCollimator;

    //G4ThreeVector        position;
    //G4RotationMatrix     rotation;

    G4double BGOShieldLength;
    G4double BGOHevimetThickness;
    G4double BGOHevimetDist;
    G4double BGOShieldDist;
    G4double BGOCrystalDist;
    G4double BGOShieldOuterRadius;
    G4double BGOShieldInnerRadius;
    G4double subInnerConeL;
    G4double abs1_thickness;
    G4double abs2_thickness;

    G4VSolid*           solidBGOCrystal;
    G4VSolid*           solidBGOHevimet; //G4SubtractionSolid*
    G4VSolid*           solidBGOShield;
    G4VSolid*           Phase1Absorber1;
    G4VSolid*           Phase1Absorber2;

private:
    void CreateSolids();

protected:
    std::vector<G4LogicalVolume*>  fScoringVolume;
};

#endif
