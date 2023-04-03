#include "NuballPhaseIHit.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4Allocator<NuballPhaseIHit>* PhaseIHitAllocator = 0;

NuballPhaseIHit::NuballPhaseIHit()
		: G4VHit(),
			fTime(0.),
			fEnergy(0.),
			fLabel(0)
{}

NuballPhaseIHit::~NuballPhaseIHit()
{}

NuballPhaseIHit::NuballPhaseIHit(const NuballPhaseIHit& /*right*/)
		: G4VHit()
{}

const NuballPhaseIHit& NuballPhaseIHit::operator=(const NuballPhaseIHit& /*right*/)
{
	return *this;
}

int NuballPhaseIHit::operator==(const NuballPhaseIHit& right) const
{
	if(fTime == right.GetTime() && fEnergy == right.GetEnergy() && fLabel == right.GetLabel())
	{
		return true;
	}
	return false;
}

void NuballPhaseIHit::Print()
{
	G4cout << "Label: " << fLabel << "; Time: " << fTime/s << "; Energy: " << fEnergy/keV << G4endl;
}
