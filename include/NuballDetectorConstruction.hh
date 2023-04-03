//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef NuballDetectorConstruction_h
#define NuballDetectorConstruction_h 1


#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Sphere;
class G4Cons;

/// Detector construction class to define materials and geometry.

class NuballDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  NuballDetectorConstruction();
  virtual ~NuballDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();


private:
  // Booleans
  G4bool fMakeLICORNE;
  G4bool fMakeOUPS;
  G4bool fMakeNEUTRONS;
  G4bool fMakePhaseI;
  G4bool fMakeClovers;
  G4bool fMakeLaBr3;
  G4bool fMakeFATIMA;
  G4bool fMakeTarget;

  G4double cell_length;
  G4double chamber_length;
  G4double translationdistance;

  // Geometry Variables
  G4double world_size;
  G4double chamberdistance;

  // Shapes
  G4Sphere* experimentalHall_solid;
  G4Cons* neutron_cone_solid;

  // Logical volumes
  G4LogicalVolume* experimentalHall_log;
  G4LogicalVolume* neutron_cone_log;


  // Physical Volumes
  G4VPhysicalVolume* experimentalHall_phys;
  G4VPhysicalVolume* neutron_cone_phys;
  // StringToRotationMatrix() converts a string "X90,Y45" into a
  // G4RotationMatrix.
  // This is an active rotation, in that the object is first rotated
  // around the parent's X axis by 90 degrees, then the object is
  // further rotated around the parent's Y axis by 45 degrees.
  // The return value points to a G4RotationMatrix on the heap, so
  // it is persistent. Angles are in degrees, can have decimals,
  // and can be negative. Axes are X, Y, Z.

  static G4RotationMatrix StringToRotationMatrix(G4String rotation);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
