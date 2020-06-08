#include <I3Test.h>

#include "MuonGun/EnsembleSampler.h"
#include "phys-services/I3GSLRandomService.h"

TEST_GROUP(EnsembleSampler);

namespace {

double gaussian2(double x, double y) {
	return -0.5*(x*x + y*y) - std::log(2*M_PI);
}

};

TEST(Initialize)
{
	typedef double (Signature)(double, double);
	typedef I3MuonGun::EnsembleSampler<Signature> Sampler;
	
	std::vector<Sampler::array_type> ensemble(8);
	for (unsigned i=0; i < ensemble.size(); i++) {
		double phi = (i*2*M_PI)/ensemble.size();
		Sampler::array_type point = {{std::cos(phi), std::sin(phi)}};
		ensemble[i] = point;
	}
	Sampler sampler(gaussian2, ensemble);
}

TEST(Gaussian2d)
{
	using namespace I3MuonGun;
	I3GSLRandomService rng(0);
	
	typedef double (Signature)(double, double);
	typedef I3MuonGun::EnsembleSampler<Signature> Sampler;
	const unsigned walkers = 64;
	const unsigned thin = 10;
	const unsigned nsamples = 100000;
	
	std::vector<Sampler::array_type> ensemble(walkers);
	for (unsigned i=0; i < ensemble.size(); i++) {
		double phi = (i*2*M_PI)/ensemble.size();
		Sampler::array_type point = {{std::cos(phi), std::sin(phi)}};
		ensemble[i] = point;
	}
	Sampler sampler(gaussian2, ensemble);
	
	for (unsigned i=0; i < thin; i++)
		sampler.Sample(rng);
	sampler.Reset();
	
	std::vector<Sampler::array_type > samples;
	samples.reserve(nsamples);
	const unsigned iterations = (thin*nsamples + walkers - 1)/walkers;
	ENSURE(iterations > thin);
	for (unsigned i=0; i < iterations; i++) {
		const std::vector<Sampler::sample> &ensemble = sampler.Sample(rng);
		if (i % thin == 0) {
			unsigned todo = std::min(walkers, nsamples - unsigned(samples.size()));
			for (unsigned j=0; j < todo; j++) {
				samples.push_back(ensemble[j].point);
			}
		}
	}
	ENSURE_EQUAL(samples.size(), nsamples);
	
	for (unsigned dim=0; dim < 2; dim++) {
		double sum = 0;
		double sum2 = 0;
		BOOST_FOREACH(const Sampler::array_type &point, samples) {
			sum += point[dim];
			sum2 += point[dim]*point[dim];
		}
		double mean = sum/samples.size();
		double var = (sum2 - (mean*mean))/(samples.size() - 1);
		ENSURE_DISTANCE(mean, 0, 1e-2, "mean is 0");
		ENSURE_DISTANCE(var, 1, 2e-2, "variance is 1");
		
	}
	ENSURE(sampler.GetAcceptanceRate() > 0);
}
