#ifndef NuballPhaseIHit_H
#define NuballPhaseIHit_H

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class NuballPhaseIHit : public G4VHit
{
	public:
		NuballPhaseIHit();
		virtual ~NuballPhaseIHit();
		NuballPhaseIHit(const NuballPhaseIHit& right);
		const NuballPhaseIHit& operator=(const NuballPhaseIHit& right);
		int operator==(const NuballPhaseIHit &right) const;

		inline void* operator new(size_t);
		inline void operator delete(void* hit);

		virtual void Print();

		//setter methods
		void SetTime(G4double time) {fTime = time;}
		void SetEnergy(G4double energy) {fEnergy = energy;}
		void SetLabel(G4int label) {fLabel = label;}

		//getter methods
		G4double GetTime() const {return fTime;}
		G4double GetEnergy() const {return fEnergy;}
		G4int GetLabel() const {return fLabel;}

	private:
		//data members
		G4double fTime;
		G4double fEnergy;
		G4int fLabel;

};

typedef G4THitsCollection<NuballPhaseIHit> PhaseIHitsCollection;

extern G4ThreadLocal G4Allocator<NuballPhaseIHit>* PhaseIHitAllocator;

inline void* NuballPhaseIHit::operator new(size_t)
{
	if(! PhaseIHitAllocator) PhaseIHitAllocator = new G4Allocator<NuballPhaseIHit>;
	return (void*)PhaseIHitAllocator->MallocSingle();
}

inline void NuballPhaseIHit::operator delete(void* hit)
{
	PhaseIHitAllocator->FreeSingle((NuballPhaseIHit*) hit);
}

#endif
