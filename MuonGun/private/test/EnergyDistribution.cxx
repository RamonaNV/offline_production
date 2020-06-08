
#include <I3Test.h>

#include "MuonGun/EnergyDistribution.h"
#include "MuonGun/RadialDistribution.h"
#include "phys-services/I3GSLRandomService.h"
#include "common.h"

TEST_GROUP(EnergyDistribution);

TEST(Equality)
{
	using namespace I3MuonGun;
	
	{
		OffsetPowerLaw p1(2, 500, 1, 1e10);
		OffsetPowerLaw p2(2, 500, 2, 1e10);
	
		ENSURE(p1 == p1);
		ENSURE(p2 == p2);
		ENSURE(!(p1 == p2));
	}
	
	{
		BundleModel p1 = load_model("Hoerandel5_atmod12_SIBYLL");
		BundleModel p2 = load_model("GaisserH4a_atmod12_SIBYLL");
		ENSURE(*p1.energy == *p1.energy);
		ENSURE(*p2.energy == *p2.energy);
		ENSURE(!(*p1.energy == *p2.energy));
	}
}

TEST(Sampling)
{
	using namespace I3MuonGun;
	typedef std::pair<double, double> pair;
	
	double depth = 2.5;
	double ct = 0.01;
	unsigned m = 100;
	
	I3GSLRandomService rng(1);
	// GaisserH4a_atmod12_SIBYLL
	// GaisserH4a_atmod12_DPMJET
	BundleModel model = load_model("GaisserH4a_atmod12_SIBYLL");
	
	{
		BMSSEnergyDistribution edist;
		std::vector<pair> vals = edist.Generate(rng, depth, ct, m, 1000u);
		BOOST_FOREACH(const pair &val, vals) {
			ENSURE(BMSSRadialDistribution()(depth, ct, m, val.first) > 0);
			double pdf = edist(depth, ct, m, val.first, val.second);
			ENSURE(pdf > 0);
		}
		
	}
	
	{
		std::vector<pair> vals = model.energy->Generate(rng, depth, ct, m, 1000);
	}
}
