#ifndef NuballCloverHit_H
#define NuballCloverHit_H

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class NuballCloverHit : public G4VHit
{
	public:
		NuballCloverHit();
		virtual ~NuballCloverHit();
		NuballCloverHit(const NuballCloverHit& right);
		const NuballCloverHit& operator=(const NuballCloverHit& right);
		int operator==(const NuballCloverHit &right) const;

		inline void* operator new(size_t);
		inline void  operator delete(void* hit);

		virtual void Print();

		//setter methods
		void SetTime(G4double time) {fTime = time;}
		void SetEnergy(G4double energy) {fEnergy = energy;}
		void SetLabel(G4int label) {fLabel = label;}

		//getter methods
		G4double GetTime() const   {return fTime;}
		G4double GetEnergy() const {return fEnergy;}
		G4int    GetLabel() const  {return fLabel;}

	private:
		//data members
		G4double fTime;
		G4double fEnergy;
		G4int fLabel;

};

typedef G4THitsCollection<NuballCloverHit> CloverHitsCollection;

extern G4ThreadLocal G4Allocator<NuballCloverHit>* CloverHitAllocator;

inline void* NuballCloverHit::operator new(size_t)
{
	if(! CloverHitAllocator) CloverHitAllocator = new G4Allocator<NuballCloverHit>;
	return (void*)CloverHitAllocator->MallocSingle();
}

inline void NuballCloverHit::operator delete(void* hit)
{
	CloverHitAllocator->FreeSingle((NuballCloverHit*) hit);
}

#endif
