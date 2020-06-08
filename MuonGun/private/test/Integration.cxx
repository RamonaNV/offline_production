
#include <I3Test.h>

#include "common.h"
#include "MuonGun/Cylinder.h"

#include <boost/make_shared.hpp>

TEST_GROUP(Integration);

using namespace I3MuonGun;

namespace {

double one(double depth, double cos_theta) {
	return 1.;
}

}

TEST(Constant)
{
	Cylinder surface(1000, 500);
	
	ENSURE_DISTANCE(surface.IntegrateFlux(one, 0., 1.),
	    surface.GetAcceptance(0., 1.), surface.GetAcceptance(0., 1.)/1e4,
	    "Numerical integration of a constant is accurate to 1e-4");
}
