#ifndef NuballLaBr3Hit_H
#define NuballLaBr3Hit_H

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class NuballLaBr3Hit : public G4VHit
{
	public:
		NuballLaBr3Hit();
		virtual ~NuballLaBr3Hit();
		NuballLaBr3Hit(const NuballLaBr3Hit& right);
		const NuballLaBr3Hit& operator=(const NuballLaBr3Hit& right);
		int operator==(const NuballLaBr3Hit &right) const;

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

typedef G4THitsCollection<NuballLaBr3Hit> LaBr3HitsCollection;

extern G4ThreadLocal G4Allocator<NuballLaBr3Hit>* LaBr3HitAllocator;

inline void* NuballLaBr3Hit::operator new(size_t)
{
	if(! LaBr3HitAllocator) LaBr3HitAllocator = new G4Allocator<NuballLaBr3Hit>;
	return (void*)LaBr3HitAllocator->MallocSingle();
}

inline void NuballLaBr3Hit::operator delete(void* hit)
{
	LaBr3HitAllocator->FreeSingle((NuballLaBr3Hit*) hit);
}

#endif
