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
// $Id: NuballRunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file NuballRunAction.cc
/// \brief Implementation of the NuballRunAction class

#include "NuballRunAction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"

#include "NuballAnalysis.hh"

//#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuballRunAction::NuballRunAction()
: G4UserRunAction()
{
    G4cout << "Initializing Ntuple" << G4endl;

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetNtupleMerging(true);
//0
    analysisManager->CreateNtuple("Nuball", "All Hits in Nuball");
    analysisManager->CreateNtupleIColumn("Label");
    analysisManager->CreateNtupleDColumn("Energy");
    analysisManager->CreateNtupleDColumn("Time");
    analysisManager->CreateNtupleIColumn("Type");
    analysisManager->CreateNtupleIColumn("EventNumber");
    analysisManager->FinishNtuple();
//1
    analysisManager->CreateNtuple("p1", "Momentum direction 1");
    analysisManager->CreateNtupleDColumn("x");
    analysisManager->CreateNtupleDColumn("y");
    analysisManager->CreateNtupleDColumn("z");
    analysisManager->CreateNtupleIColumn("EvtID");
    analysisManager->FinishNtuple();
//2
    analysisManager->CreateNtuple("p2", "Momentum direction 2");
    analysisManager->CreateNtupleDColumn("x");
    analysisManager->CreateNtupleDColumn("y");
    analysisManager->CreateNtupleDColumn("z");
    analysisManager->CreateNtupleIColumn("EvtID");
    analysisManager->FinishNtuple();
//3
    analysisManager->CreateNtuple("gun", "Energy & Multiplicity distribution");
    analysisManager->CreateNtupleDColumn("Energy");
    analysisManager->CreateNtupleIColumn("Multiplicity");
    analysisManager->CreateNtupleIColumn("EvtID");
    analysisManager->FinishNtuple();

    analysisManager->SetNtupleActivation(false);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuballRunAction::~NuballRunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NuballRunAction::BeginOfRunAction(const G4Run*)
{

  G4cout << "Begin of Run Action" << G4endl;
  // inform the runManager to save random number seed

  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  

  // Open an output file
  //
  G4String fileName = "Nuball";
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile(fileName);

  if(IsMaster()) // Only for the master thread (in case of multithread mode)
  {
    //Initialize the seed for the random generator
    long seeds[2];
    auto systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    G4Random::setTheSeeds(seeds);

    G4cout << G4endl
    << " -------------------Beginning--of--Run--------------------------"
    << G4endl;
  }
  //----------------------------------------------------------------------------//

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NuballRunAction::EndOfRunAction(const G4Run*)
{
  // save histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
