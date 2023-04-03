#include "NuballCloverSD.hh"
#include "NuballAnalysis.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4VTouchable.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <G4INCLRandom.hh>
#include "TRandom3.h"
#include "G4MTRunManager.hh"

NuballCloverSD::NuballCloverSD(const G4String& name,
                         const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name),
   fHitsCollection(0)
{
  randgen = new TRandom3();
  collectionName.insert(hitsCollectionName);
}

NuballCloverSD::~NuballCloverSD()
{}

void NuballCloverSD::Initialize(G4HCofThisEvent* hce)
{
	fHitsCollection = new CloverHitsCollection(SensitiveDetectorName, collectionName[0]);
	G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	hce->AddHitsCollection( hcID, fHitsCollection);

  for (G4int i=0; i<144; ++i) {
    NuballCloverHit* newHit = new NuballCloverHit();
    newHit->SetLabel(23+i);
    fHitsCollection->insert(newHit);
  }
}

G4bool NuballCloverSD::ProcessHits(G4Step* step, G4TouchableHistory* /*history*/)
{
  G4StepPoint* preStepPoint = step->GetPreStepPoint();

  //Label
  G4int label = preStepPoint->GetTouchable()->GetCopyNumber();

  // Get hit accounting data for this layer
  NuballCloverHit* hit = (*fHitsCollection)[label-23];
  if ( ! hit ) {
    G4cerr << "Cannot access hit of label " << label << G4endl;
    exit(1);
  }

	//Time
  if(hit->GetTime() == 0)
  {
  	G4double time = preStepPoint->GetGlobalTime();
  	hit->SetTime(time);
  }
	//Energy
	G4double energy = step->GetTotalEnergyDeposit();
	hit->SetEnergy(hit->GetEnergy()+energy);

	return false;
}


void NuballCloverSD::EndOfEvent(G4HCofThisEvent* /*hce*/)
{

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  for(int i = 0; i <144; i++)
  {
    NuballCloverHit* hit = (NuballCloverHit*)(*fHitsCollection)[i];
    if(hit->GetEnergy()>0)
    {
      if(i%6 == 1)
      {
        //hit->SetEnergy(randgen->Gaus(hit->GetEnergy(), 0.1148 * hit->GetEnergy()));
        analysisManager->FillNtupleIColumn(3, 0);
      }
      else if(i%6 > 1)
      {
        //hit->SetEnergy(randgen->Gaus(hit->GetEnergy(), 4.742e-4 * hit->GetEnergy() + 2.08e-4));
        analysisManager->FillNtupleIColumn(3, 1);
      }
      analysisManager->FillNtupleIColumn(0, hit->GetLabel());
      analysisManager->FillNtupleIColumn(4, G4MTRunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
      analysisManager->FillNtupleDColumn(1, hit->GetEnergy());
      analysisManager->FillNtupleDColumn(2, hit->GetTime());

      analysisManager->AddNtupleRow();
    }
  }
}
