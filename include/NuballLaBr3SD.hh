#ifndef NuballLaBr3SD_H
#define NuballLaBr3SD_H

#include "G4VSensitiveDetector.hh"
#include "NuballLaBr3Hit.hh"
#include "Randomize.hh"
#include "TRandom3.h"
#include <G4INCLRandom.hh>


class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class NuballLaBr3SD : public G4VSensitiveDetector
{
	public:
		NuballLaBr3SD(const G4String& name, const G4String& hitsCollectionName);
		virtual ~NuballLaBr3SD();

		virtual void Initialize(G4HCofThisEvent* hce);
		virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
		virtual void EndOfEvent(G4HCofThisEvent* hce);

	private:
		LaBr3HitsCollection* fHitsCollection;
		TRandom3* randgen;
};



#endif //NuballLaBr3SD_H
