
#include <I3Test.h>
#include "MuonGun/SplineTable.h"

#include "common.h"
#include "MuonGun/Flux.h"
#include "MuonGun/RadialDistribution.h"
#include "MuonGun/EnergyDistribution.h"
#include <boost/make_shared.hpp>

namespace I3MuonGun {

std::string
get_tabledir()
{
	const char *I3_BUILD = getenv("I3_BUILD");
	ENSURE(I3_BUILD, "I3_BUILD is set");
	
	return std::string(I3_BUILD) + "/MuonGun/resources/tables/";
}

BundleModel
load_model(const std::string &base)
{
	using boost::make_shared;
	return BundleModel(boost::make_shared<SplineFlux>(get_tabledir() + base + ".single_flux.fits", get_tabledir() + base + ".bundle_flux.fits"),
	    boost::make_shared<SplineRadialDistribution>(get_tabledir() + base + ".radius.fits"),
	    boost::make_shared<SplineEnergyDistribution>(get_tabledir() + base + ".single_energy.fits", get_tabledir() + base + ".bundle_energy.fits"));
}

}


TEST_GROUP(SplineTable);

TEST(Equality)
{
	using namespace I3MuonGun;
	
	const SplineTable t1(get_tabledir() + "Hoerandel5_atmod12_SIBYLL.single_flux.fits");
	const SplineTable t2(get_tabledir() + "GaisserH4a_atmod12_SIBYLL.single_flux.fits");
	
	ENSURE(t1 == t1);
	ENSURE(t2 == t2);
	ENSURE(!(t1 == t2));
	ENSURE(!(t2 == t1));
}
