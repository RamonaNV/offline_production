
#include <I3Test.h>

#include "common.h"
#include "MuonGun/Generator.h"
#include "MuonGun/CORSIKAGenerationProbability.h"
#include "MuonGun/Cylinder.h"

#include <boost/make_shared.hpp>

TEST_GROUP(Generator);

TEST(Addition)
{
	using namespace I3MuonGun;
	using boost::make_shared;
	
	BundleModel soft = load_model("Hoerandel5_atmod12_SIBYLL");
	BundleModel hard = load_model("GaisserH4a_atmod12_SIBYLL");
	
	boost::shared_ptr<CORSIKAGenerationProbability> soft_g(
	    new CORSIKAGenerationProbability(make_shared<Cylinder>(1600, 800), soft.flux, soft.radius, soft.energy));
	boost::shared_ptr<CORSIKAGenerationProbability> hard_g(
	    new CORSIKAGenerationProbability(make_shared<Cylinder>(1600, 800), hard.flux, hard.radius, hard.energy));
	
	ENSURE(soft_g->IsCompatible(soft_g));
	ENSURE_EQUAL((soft_g+soft_g)->GetTotalEvents(), 2.);
	ENSURE((bool)boost::dynamic_pointer_cast<const CORSIKAGenerationProbability>(soft_g+soft_g),
	    "Sum of compatible generators is a generator of the same type");
	ENSURE((bool)boost::dynamic_pointer_cast<const CORSIKAGenerationProbability>(soft_g+soft_g+soft_g),
	    "Sum of compatible generators is a generator of the same type");
	
	ENSURE(!soft_g->IsCompatible(hard_g));
	ENSURE_EQUAL((soft_g+hard_g)->GetTotalEvents(), 1.);
	ENSURE((bool)boost::dynamic_pointer_cast<const GenerationProbabilityCollection>(soft_g+hard_g),
	    "Sum of incompatible generators is a collection");
}
