#ifndef PhaseI_H
#define PhaseI_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include<vector>

class G4Polycone;
class G4UnionSolid;

class G4LogicalVolume;
class G4VPhysicalVolume;

class PhaseI {

public:
    PhaseI();
    ~PhaseI();

public:
    //void SetPosition( G4ThreeVector );
    //void SetRotation( G4RotationMatrix );
    void Placement(G4int, G4VPhysicalVolume*, G4bool, G4ThreeVector, G4RotationMatrix, int);

    void SetGeRingZPos   (G4double);

    std::vector<G4LogicalVolume*> GetScoringVolumes() const { return fScoringVolume; }
    G4LogicalVolume* GetScoringVolume(int i) const {return fScoringVolume.at(i);}


private:

    //G4ThreeVector        position;
    //G4RotationMatrix     rotation;

    //-------------------------------------------
    // A Phase-I tapered Ge
    G4double            fGeTaperL_PhaseI;
    G4double            fTotalGeL_PhaseI;
    G4double            fGeOuterD_PhaseI;
    G4double            fGeInnerD_PhaseI;
    G4double            fAlCap2Ge_PhaseI;
    G4double            fEndCapThickness_PhaseI;
    G4double            fEndCapTaperL_PhaseI;
    G4double            fEndCapTubeL_PhaseI;
    G4double            fEndCapFrontR_PhaseI;
    G4double            fEndCapBackR_PhaseI;
    //
    G4Polycone*         solidAlCap_PhaseI;  //Al end-cap
    G4Polycone*         solidVacuum_PhaseI; //vacuum
    G4Polycone*         solidGe_PhaseI; //Ge
    G4double           fContact_dZ_PhaseI;  //dZ to position Passivated Ge
    G4UnionSolid*      solidPassivated_PhaseI;
    G4UnionSolid*      solidContact_PhaseI; //inner Li contact
    G4UnionSolid*      solidBoreHole_PhaseI; //inner Li contact
private:

    void CreatePhaseISolids();

    std::vector<G4LogicalVolume*> fScoringVolume;

};

#endif
