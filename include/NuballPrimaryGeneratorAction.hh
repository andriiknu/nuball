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
// $Id: NuballPrimaryGeneratorAction.hh 90623 2015-06-05 09:24:30Z gcosmo $
//
/// \file NuballPrimaryGeneratorAction.hh
/// \brief Definition of the NuballPrimaryGeneratorAction class

#ifndef NuballPrimaryGeneratorAction_h
#define NuballPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "TTree.h"
#include "TFile.h"
#include<vector>
class G4ParticleGun;
class G4Event;
class G4Box;

class NuballPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  NuballPrimaryGeneratorAction();
  virtual ~NuballPrimaryGeneratorAction();

  // method from the base class
  virtual void GeneratePrimaries(G4Event*);

  // method to access particle gun
  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  void SetMonoenenergeticGamma(G4double energy, G4double time);
  void SetSourceDecay(G4int A, G4int Z);
  void SetSourceDecay(G4int A, G4int Z, G4double E, G4int J);
  void SetThoriumDecay();
  void SetFifrelinDecay(G4Event* anEvent);
  void SetMultipleGamma(std::function<int(G4double*)> GetRandomMultiplicity, std::function<double(/* int, double, double */)>, G4Event* anEvent);
  G4bool IsSourceConfined(G4ThreeVector);
  void SetCo60decay(G4Event* anEvent);

private:
  G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
  G4Box* fWorldBox;
  TFile* fifrelin_file;
  TTree* fifrelin_tree;
	Long_t fifrelin_index = -1;

	std::vector<Double_t>* E_gamma;
	std::vector<Double_t>* E_neutron;
  G4double *Cumu;
  G4double *MultiplicityCumu;

  G4ThreeVector* interaction;

  Long_t compteur = 0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
