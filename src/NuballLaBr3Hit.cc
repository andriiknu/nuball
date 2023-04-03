#include "NuballLaBr3Hit.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4Allocator<NuballLaBr3Hit>* LaBr3HitAllocator = 0;

NuballLaBr3Hit::NuballLaBr3Hit()
		: G4VHit(),
			fTime(0.),
			fEnergy(0.),
			fLabel(0)
{}

NuballLaBr3Hit::~NuballLaBr3Hit()
{}

NuballLaBr3Hit::NuballLaBr3Hit(const NuballLaBr3Hit& /*right*/)
		: G4VHit()
{}

const NuballLaBr3Hit& NuballLaBr3Hit::operator=(const NuballLaBr3Hit& /*right*/)
{
	return *this;
}

int NuballLaBr3Hit::operator==(const NuballLaBr3Hit& right) const
{
	if(fTime == right.GetTime() && fEnergy == right.GetEnergy() && fLabel == right.GetLabel())
	{
		return true;
	}
	return false;
}

void NuballLaBr3Hit::Print()
{
	G4cout << "Label: " << fLabel << "; Time: " << fTime/s << "; Energy: " << fEnergy/keV << G4endl;
}
