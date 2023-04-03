#ifndef NuballLaBr3_H
#define NuballLaBr3_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include<vector>


class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;


class LaBr3 {

public:
  LaBr3();
  ~LaBr3();

public:
  void Placement(G4int, G4VPhysicalVolume*, G4bool, G4ThreeVector, G4RotationMatrix, int);

public:
  inline G4double GetPMRadius() {return fPMouterR;};
  inline G4double GetCrystalRadius() {return fCrystalR;};
  inline G4double GetCrystalThickness() {return fCrystalThickness;};
  inline G4double GetCapThickness() {return fCapThickness;};

  std::vector<G4LogicalVolume*> GetLogicalVolumes() const {return fLogical;}
  G4LogicalVolume* GetLogicalVolume(int i) const {return fLogical.at(i);}

private:

  G4double             lappingSize;
  G4double             LGlappingSize;
  G4double             PMlappingSize;
  G4double             alulappingsize;

  // Solid parameters
  G4double             fBoxsizez;
  G4double             fBoxsizexy;
  G4double             fCapR;
  G4double             fCapD;
  G4double             fCapThickness;
  G4double             fCrystalR;
  G4double             fCrystalD;
  G4double             fCrystalThickness;
  G4double             fPMinnerR;
  G4double             fPMouterR;
  G4double             fPMinnerD;
  G4double             fPMouterD;
  G4double             fPMThickness;
  G4double             flightguideR;
  G4double             flightguideThickness;

  // Shapes
  G4Box* LaBr_solid;
  G4Tubs* LaBr3_solid;
  G4Tubs* Alcover_solid;
  G4Tubs* lightguide_solid;
  G4Tubs* PM_solid;

protected:
  std::vector<G4LogicalVolume*> fLogical;

private:
  void CreateSolids();
};

#endif
