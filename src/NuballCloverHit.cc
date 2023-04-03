#include "NuballCloverHit.hh"

#include "G4SystemOfUnits.hh"

G4ThreadLocal G4Allocator<NuballCloverHit>* CloverHitAllocator = 0;

NuballCloverHit::NuballCloverHit()
		: G4VHit(),
			fTime(0.),
			fEnergy(0.),
			fLabel(0)
{}

NuballCloverHit::~NuballCloverHit()
{}

NuballCloverHit::NuballCloverHit(const NuballCloverHit& /*right*/)
		: G4VHit()
{}

const NuballCloverHit& NuballCloverHit::operator=(const NuballCloverHit& /*right*/)
{
	return *this;
}

int NuballCloverHit::operator==(const NuballCloverHit& right) const
{
	if(fTime == right.GetTime() && fEnergy == right.GetEnergy() && fLabel == right.GetLabel())
	{
		return true;
	}
	return false;
}

void NuballCloverHit::Print()
{
	G4cout << "Label: " << fLabel << "; Time: " << fTime/s << "; Energy: " << fEnergy/keV << G4endl;
}
