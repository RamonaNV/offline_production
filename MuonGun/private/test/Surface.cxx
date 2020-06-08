
#include <I3Test.h>

#include "MuonGun/Cylinder.h"

TEST_GROUP(Surface);

TEST(Equality)
{
	using namespace I3MuonGun;
	
	Cylinder cylinder(1600, 800);
	Cylinder offset_cylinder(1600, 800, I3Position(3,2,1));
	// Sphere sphere(6e3, 900);
	
	ENSURE(cylinder == cylinder);
	ENSURE(offset_cylinder == offset_cylinder);
	// ENSURE(sphere == sphere);
	
	ENSURE(!(cylinder == offset_cylinder));
	// ENSURE(!(cylinder == sphere));
}
