#ifndef NuballPhaseISD_H
#define NuballPhaseISD_H

#include "G4VSensitiveDetector.hh"
#include "NuballPhaseIHit.hh"
#include "Randomize.hh"
#include "TRandom3.h"
#include <G4INCLRandom.hh>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class NuballPhaseISD : public G4VSensitiveDetector
{
	public:
		NuballPhaseISD(const G4String& name, const G4String& hitsCollectionName);
		virtual ~NuballPhaseISD();

		virtual void Initialize(G4HCofThisEvent* hce);
		virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
		virtual void EndOfEvent(G4HCofThisEvent* hce);

	private:
		PhaseIHitsCollection* fHitsCollection;
		TRandom3* randgen;

};



#endif //NuballPhaseISD_H
