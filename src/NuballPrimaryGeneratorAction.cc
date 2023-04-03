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
// $Id: NuballPrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file NuballPrimaryGeneratorAction.cc
/// \brief Implementation of the NuballPrimaryGeneratorAction class

#include "NuballPrimaryGeneratorAction.hh"

//#include "G4LogicalVolumeStore.hh"
//#include "G4LogicalVolume.hh"
//#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "NuballAnalysis.hh"
#include "matrix.hh"

#include <cmath>
#include<vector>

// #define NDEBUG

Double_t ScissionSpinFunction(double x, double avgspin)
{
	double sig=avgspin/1.15;
	double J=x;
	double p=(((2*J)+1)/(2*sig*sig))*exp(-((J+0.5)*(J+0.5)/(2.0*sig*sig)));
	return p;
}

double* GetMultiplicityCumuFunc() {
  double avgspin = 8;
  double tot=0;
  int n = 20;
  double *ScissionSpinDist=new double[20];
  for (int i=0; i < n; i++) 
	{
    ScissionSpinDist[i]=ScissionSpinFunction(i, avgspin);
    assert(ScissionSpinDist[i] != 0);
    tot+=ScissionSpinDist[i];

	}
  
  for (int i=0; i < n; i++) {ScissionSpinDist[i]/=tot;}
  double *cumu=new double[20];
  cumu[0]=ScissionSpinDist[0];
  assert(cumu[0] != 0);
  for (int i=1; i < n; i++) {
    cumu[i]=cumu[i-1]+ScissionSpinDist[i];
    
  }


  
  delete [] ScissionSpinDist;
  return cumu;
  
}

G4int GetMultiplicityValue(G4double *cumu) {
  int n = 20;

#ifndef NDEBUG
  for (int i = 0; i < 20; ++i) assert(cumu[i] != 0);
#endif



  int M = 0;
  double r1=drand48()*(1 - cumu[0]) + cumu[0];
  for (int j=0; j < n; j++) {if (r1 < cumu[j]) {M=j; break;}}
  assert(M != 0);
  return M;
}

G4double GetFLATdistrEnergy() {
  int _integer = 18*G4UniformRand();
  G4double particleEnergy = 0.2*MeV + ((double)_integer/10)*MeV;
  return particleEnergy;
}

int n = 180;

double w (double x) {
    double a0 = 1.;
    double a2 = 1./8.;
    double a4 = 1./24.;
    return a0+a2*pow(cos(x), 2) + a4*pow(cos(x), 4);
}

double *GetCumuletiveFunction () {
    
    double *p = new double[n];
    double Sum = 0;
    for (int i = 0; i < n; i++) {
        double x = i * deg;
        p[i] = w(x);
        Sum += p[i];
        assert(p[i] > 0);
    }
    for (int i = 0; i < n; ++i) p[i] /= Sum;
    

    double *cumu = new double[n];
    cumu[0] = p[0];
    for (int i = 1; i < n; ++i) cumu[i] = cumu[i-1] + p[i]; 
    delete [] p;   
    return cumu;
}

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuballPrimaryGeneratorAction::NuballPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0),
fWorldBox(0)
{
  G4int n_particle = 1; //Number of particle to be generated
  fParticleGun  = new G4ParticleGun(n_particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  interaction = new G4ThreeVector();
  Cumu = GetCumuletiveFunction();
  MultiplicityCumu = GetMultiplicityCumuFunc();
  G4Random::setTheSeed(time(0));

#ifndef NDEBUG
std::cout << "DEBUG = ON\n";
  for (int i = 0; i < 20; ++i) assert(MultiplicityCumu[i] != 0);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuballPrimaryGeneratorAction::~NuballPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete Cumu;
  delete MultiplicityCumu;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NuballPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // SetCo60decay(anEvent);
  SetMultipleGamma(&GetMultiplicityValue, &GetFLATdistrEnergy, anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NuballPrimaryGeneratorAction::SetMonoenenergeticGamma(G4double energy, G4double time)
{
  fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
  G4double theta = acos(1-2*G4UniformRand());
  G4double phi = 2*M_PI*G4UniformRand();

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi)
  ,sin(theta)*sin(phi)
  ,cos(theta)));
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->SetParticleTime(time);
}

void NuballPrimaryGeneratorAction::SetMultipleGamma(std::function<G4int(G4double*)> GetRandomMultiplicity, std::function<double(/* int, double, double */)> GetRandomEnergy, G4Event* anEvent) {
  fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
  G4double Energy = 0;
  G4int Multiplicity = GetRandomMultiplicity(MultiplicityCumu);
  if (Multiplicity == 0) std::cout << "Multiplicity == 0\n";

  for (int i  = 0; i < Multiplicity; ++i) {
    G4double particleEnergy = GetRandomEnergy();
    fParticleGun->SetParticleEnergy(particleEnergy);
    G4double theta = acos(1 - 2*G4UniformRand());
    G4double phi = 2*M_PI*G4UniformRand();
    fParticleGun->SetParticleMomentumDirection (G4ThreeVector(sin(theta)*cos(phi)
                                                ,sin(theta)*sin(phi)
                                                ,cos(theta)));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    Energy +=  particleEnergy; 
  }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(3, 0, Energy);
  analysisManager->FillNtupleIColumn(3, 1, Multiplicity);
  analysisManager->FillNtupleIColumn(3, 2, anEvent->GetEventID());
  analysisManager->AddNtupleRow(3);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NuballPrimaryGeneratorAction::SetCo60decay(G4Event* anEvent) {

  fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
  G4double theta = acos(1-2*G4UniformRand());
  G4double phi = 2*M_PI*G4UniformRand();

  G4double alpha0 = 0;
  srand48(time(0)*drand48());
  double u = drand48();//*(1.0 - cumu[0]) + cumu[0];
  u = G4UniformRand();
  for (int i = 0; i < n; ++i) {
    if (u <= Cumu[i]) {
      alpha0 = i*deg;
      break;
    }
    
  }
  alpha0 += 0.1 * deg;
  assert(1);
  G4double phi_ = 2*M_PI*G4UniformRand();
  TMatrixD g1_ket_ = GetKetVector(alpha0/2, phi_);
  double g2_c[] = {-g1_ket_.GetMatrixArray()[0], -g1_ket_.GetMatrixArray()[1], g1_ket_.GetMatrixArray()[2]};
  double *g1_c = g1_ket_.GetMatrixArray();
  TMatrixD g2_ket_(3, 1, g2_c);
  assert(abs(TVector3(g1_c).Angle(TVector3(g2_c)) - alpha0) < 0.000001);
  assert(abs(TVector3(g1_ket_.GetMatrixArray()).Angle(TVector3(g2_ket_.GetMatrixArray())) - alpha0) < 0.000001);
  
  TMatrixD trM = GetTransMatrix(theta, phi);
  trM.Transpose(trM);
  auto g1_ket = trM * g1_ket_;
  auto g2_ket = trM * g2_ket_;
  TVector3 g1_v(g1_ket.GetMatrixArray());
  TVector3 g2_v(g2_ket.GetMatrixArray());
  double alpha_ = TVector3(g1_ket_.GetMatrixArray()).Angle(TVector3(g2_ket_.GetMatrixArray()));
  double d1 = abs(alpha_ - alpha0);
  double d2 = abs(g1_v.Angle(g2_v) - alpha0);
  double eps = 0.01;
  assert(d1 < eps);
  assert(d2 < eps);

  auto analysisManager = G4AnalysisManager::Instance();

  double x = g1_v.X(); double y = g1_v.Y(); double z = g1_v.Z();
  G4double Energy = 1.1732;
  fParticleGun->SetParticleEnergy(Energy);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x, y, z));
  fParticleGun->GeneratePrimaryVertex(anEvent);


  analysisManager->FillNtupleDColumn(1, 0, x);
  analysisManager->FillNtupleDColumn(1, 1, y);
  analysisManager->FillNtupleDColumn(1, 2, z);
  analysisManager->FillNtupleIColumn(1, 3, anEvent->GetEventID());
  analysisManager->AddNtupleRow(1);

  x = g2_v.X(); y = g2_v.Y(); z = g2_v.Z();
  Energy = 1.3325;
  fParticleGun->SetParticleEnergy(Energy);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x, y, z));
  fParticleGun->GeneratePrimaryVertex(anEvent);



  analysisManager->FillNtupleDColumn(2, 0, x);
  analysisManager->FillNtupleDColumn(2, 1, y);
  analysisManager->FillNtupleDColumn(2, 2, z);
  analysisManager->FillNtupleIColumn(2, 3, anEvent->GetEventID());
  analysisManager->AddNtupleRow(2);




}
void NuballPrimaryGeneratorAction::SetSourceDecay(G4int A, G4int Z)
{
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z, A, 1*keV, 0);
  if(ion)
  {
    
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(0);
    fParticleGun->SetParticleEnergy(0*eV); //<== It is necessary to set the energy
    fParticleGun->SetParticleTime(0); // <== To be studied
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,0));
  }
}

void NuballPrimaryGeneratorAction::SetSourceDecay(G4int A, G4int Z, G4double E, G4int J)
{
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z, A, E*keV, J);
  if(ion)
  {
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(0);
    fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
    fParticleGun->SetParticleTime(0); // <== To be studied
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,0));
  }
}

void NuballPrimaryGeneratorAction::SetThoriumDecay()
{

  G4double rand = G4UniformRand();

  if(rand < 0.1)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(90, 232, 1*keV, 0);
    if(ion)
    {
      //G4cout << rand << G4endl;
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.2)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(90, 228, 1*keV, 0);
    if(ion)
    {
      //G4cout << rand << G4endl;
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.3)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(88, 228, 1*keV, 0);
    if(ion)
    {
      //G4cout << rand << G4endl;
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.4)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(89, 228, 1*keV, 0);
    if(ion)
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.5)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(88, 224, 1*keV, 0);
    if(ion)
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.6)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(86, 220, 1*keV, 0);
    if(ion)
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.7)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(84, 216, 1*keV, 0);
    if(ion)
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.8)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(82, 212, 1*keV, 0);
    if(ion)
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.9)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(83, 212, 1*keV, 0);
    if(ion)
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else if(rand < 0.9666666)
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(84, 212, 1*keV, 0);
    if(ion)
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
  else
  {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(81, 208, 1*keV, 0);
    if(ion)
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(0);
      fParticleGun->SetParticleEnergy(0*keV); //<== It is necessary to set the energy
      fParticleGun->SetParticleTime(0); // <== To be studied
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
  }
}

void NuballPrimaryGeneratorAction::SetFifrelinDecay(G4Event* anEvent)
{
  if(fifrelin_index>-1)
  {
    fifrelin_tree->GetEntry(fifrelin_index);
    fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("neutron"));

    for(unsigned int i = 0; i < E_neutron->size(); i++)
    {
      G4double theta = acos(1-2*G4UniformRand());
      G4double phi = 2*M_PI*G4UniformRand();

      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi)
      ,sin(theta)*sin(phi)
      ,cos(theta)));
      fParticleGun->SetParticleEnergy(E_neutron->at(i)*MeV);
      fParticleGun->SetParticleTime(0);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
    for(unsigned int i = 0; i < E_gamma->size(); i++)
    {
      G4double theta = acos(1-2*G4UniformRand());
      G4double phi = 2*M_PI*G4UniformRand();

      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi)
      ,sin(theta)*sin(phi)
      ,cos(theta)));
      fParticleGun->SetParticleEnergy(E_gamma->at(i)*MeV);
      fParticleGun->SetParticleTime(0);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fifrelin_index++;
  }
  else
  {
    fifrelin_file = TFile::Open("fifrelin30ns.root");
    if(fifrelin_file->IsOpen())
    {
      fifrelin_tree = static_cast<TTree*>(fifrelin_file->Get("fifrelin30ns"));

      E_gamma = new vector<Double_t>;
      E_neutron = new vector<Double_t>;

      fifrelin_tree->SetBranchAddress("E_gamma", &E_gamma);
      fifrelin_tree->SetBranchAddress("E_neutron", &E_neutron);

      fifrelin_index = 0;

       fifrelin_tree->GetEntry(fifrelin_index);
    fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("neutron"));

    for(unsigned int i = 0; i < E_neutron->size(); i++)
    {
      G4double theta = acos(1-2*G4UniformRand());
      G4double phi = 2*M_PI*G4UniformRand();

      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi)
      ,sin(theta)*sin(phi)
      ,cos(theta)));
      fParticleGun->SetParticleEnergy(E_neutron->at(i)*MeV);
      fParticleGun->SetParticleTime(0);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
    for(unsigned int i = 0; i < E_gamma->size(); i++)
    {
      G4double theta = acos(1-2*G4UniformRand());
      G4double phi = 2*M_PI*G4UniformRand();

      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi)
      ,sin(theta)*sin(phi)
      ,cos(theta)));
      fParticleGun->SetParticleEnergy(E_gamma->at(i)*MeV);
      fParticleGun->SetParticleTime(0);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fifrelin_index++;

    }
    else
    {
      G4cerr << "Fifrelin file not found ... " << G4endl;
      exit(-1);
    }
  }
}

G4bool NuballPrimaryGeneratorAction::IsSourceConfined(G4ThreeVector particle_position)
{
  // Method to check point is within the volume specified
  G4Navigator* gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  G4ThreeVector null(0.,0.,0.);
  G4ThreeVector *ptr;
  ptr = &null;

  // Check particle_position is within VolName, if so true,
  // else false
  G4VPhysicalVolume *theVolume;
  theVolume=gNavigator->LocateGlobalPointAndSetup(particle_position,ptr,true);
  G4String theVolName = theVolume->GetName();
   //G4cout << "VolName : " << theVolName << G4endl;
  if(theVolName == "Target_1" || theVolName == "Target_2" || theVolName == "Target_3" || theVolName == "Target_4" || theVolName == "Target_5" || theVolName == "Target_6")
    {
      return(true);
    }
  else
    return(false);
}
