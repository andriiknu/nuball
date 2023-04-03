#ifndef NuballCloverSD_H
#define NuballCloverSD_H

#include "G4VSensitiveDetector.hh"
#include "NuballCloverHit.hh"
#include "Randomize.hh"
#include "TRandom3.h"
#include <G4INCLRandom.hh>


class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class NuballCloverSD : public G4VSensitiveDetector
{
	public:
		NuballCloverSD(const G4String& name, const G4String& hitsCollectionName);
		virtual ~NuballCloverSD();

		virtual void Initialize(G4HCofThisEvent* hce);
		virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
		virtual void EndOfEvent(G4HCofThisEvent* hce);

	private:
		CloverHitsCollection* fHitsCollection;
		TRandom3* randgen;

};



#endif //NuballCloverSD_H
