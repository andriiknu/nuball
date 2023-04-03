#include "NuballPhaseISD.hh"
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


NuballPhaseISD::NuballPhaseISD(const G4String& name,
                         const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name),
   fHitsCollection(0)
{
  randgen = new TRandom3();
  collectionName.insert(hitsCollectionName);
}

NuballPhaseISD::~NuballPhaseISD()
{}

void NuballPhaseISD::Initialize(G4HCofThisEvent* hce)
{
fHitsCollection = new PhaseIHitsCollection(SensitiveDetectorName, collectionName[0]);
G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
hce->AddHitsCollection( hcID, fHitsCollection);

  for (G4int i=0; i<20; ++i)
  {
    NuballPhaseIHit* newHit = new NuballPhaseIHit();
    newHit->SetLabel(i+1);
    fHitsCollection->insert(newHit);
  }
}

G4bool NuballPhaseISD::ProcessHits(G4Step* step, G4TouchableHistory* /*history*/)
{
  G4StepPoint* preStepPoint = step->GetPreStepPoint();

  //Label
  G4int label = preStepPoint->GetTouchable()->GetCopyNumber();

  // Get hit accounting data for this layer
  NuballPhaseIHit* hit = (*fHitsCollection)[label-1];
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


void NuballPhaseISD::EndOfEvent(G4HCofThisEvent* /*hce*/)
{

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  for(int i = 0; i <20; i++)
  {
    NuballPhaseIHit* hit = static_cast<NuballPhaseIHit*>((*fHitsCollection)[i]);
    if(hit->GetEnergy()>0)
    {
      if(hit->GetLabel()%2 == 0)
      {
        //hit->SetEnergy(randgen->Gaus(hit->GetEnergy(), 6.683e-4 * hit->GetEnergy() + 0.000488599));
        analysisManager->FillNtupleIColumn(3, 1);
      }
      else
      {
        //hit->SetEnergy(randgen->Gaus(hit->GetEnergy(), 0.1148 * hit->GetEnergy()));
        analysisManager->FillNtupleIColumn(3, 0);
      }

      analysisManager->FillNtupleIColumn(0, hit->GetLabel());
      analysisManager->FillNtupleDColumn(1, hit->GetEnergy());
      analysisManager->FillNtupleDColumn(2, hit->GetTime());
      analysisManager->FillNtupleIColumn(4, G4MTRunManager::GetRunManager()->GetCurrentEvent()->GetEventID());

      analysisManager->AddNtupleRow();
    }
  }

}
