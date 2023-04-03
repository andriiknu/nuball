#ifndef NuballLaBr3_valencia_H
#define NuballLaBr3_valencia_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4VSolid.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include<vector>

class G4UnionSolid;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;

class LaBr3_valencia {

public:
    LaBr3_valencia();
    ~LaBr3_valencia();

public:
    void Placement(G4int, G4VPhysicalVolume*, G4bool, G4ThreeVector, G4RotationMatrix, int);

public:
    inline G4double GetPMRadius(){return fPMR;};
    inline G4double GetCrystalRadius(){return fCrystalR;};
    inline G4double GetCrystalThickness(){return fCrystalThickness;};
    inline G4double GetCapThickness(){return fCapThickness;};
    inline G4double GetControlThickness(){return fcontrolthickness;};

    std::vector<G4LogicalVolume*> GetLogicalVolumes() const {return fLogical;}
    G4LogicalVolume* GetLogicalVolume(int i) const {return fLogical.at(i);}

private:

    G4double             lappingSize;
    G4double             Refllappingsize;
    G4double             LGlappingSize;
    G4double             PMlappingSize;
    G4double             alulappingsize;

    // Solid parameters
    //Definition of Lead Shielding
    G4double fshieldsmallerdiameter;
    G4double fshieldinnersmallerdiameter;
    G4double fshieldgreaterdiameter;
    G4double fshieldinnergreaterdiameter;
    G4double fshieldconethickness;

    //Definition of Aluminium Cap
    G4double fcapsmallerdiameter;
    G4double fcapgreaterdiameter;
    G4double fcapconethickness;
    G4double ffrontthickness;
    G4double fCapThickness,fCapD,fCapR;

    // Definition of Reflector parameters
    G4double fRsmallerdiameter;
    G4double fRgreaterdiameter;
    G4double fRconethickness;

    // Definition of crystal parameters
    G4double fsmallerdiameter;
    G4double fgreaterdiameter;
    G4double fconethickness;
    G4double fCrystalD,fCrystalR,fCrystalThickness;


    // Definition of PM coat parameters
    G4double fPMcoatD;
    G4double fPMcoatThickness;

    // Definition of PM parameters
    G4double fPMinnerD,fPMouterD,fPMR;
    G4double fPMThickness;

    // Definition of light guide parameters
    G4double flightguideThickness;
    G4double flightguideDiameter;


    // Definition of Control Volume
    G4double fcontroldiameter;
    G4double fcontrolthickness;


    // Shapes
    G4Tubs* LaBr_solid;
    G4VSolid* Pbcover_cone;
    G4VSolid* Alcover_cone;
    G4VSolid* Alcover_front;
    G4UnionSolid* Alcover_solid;
    G4VSolid* Refl_solid;
    G4VSolid* LaBr3_solid;
    G4Tubs* lightguide_solid;
    G4Tubs* PMCap_solid;
    G4Tubs* PM_solid;

protected:
    std::vector<G4LogicalVolume*> fLogical;
private:
    void CreateSolids();

};

#endif
