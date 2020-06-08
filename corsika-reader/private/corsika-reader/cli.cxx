
#include "corsika-reader/I3CORSIKAService.h"
#include "I3CORSIKAReaderUtils.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>

#include <dataclasses/physics/I3MCTreeUtils.h>
#include <phys-services/surfaces/Cylinder.h>
#include <phys-services/I3GSLRandomService.h>

#include <neutrino-generator/legacy/I3NuGIncrementalEventGenerator.h>
#include <neutrino-generator/Steering.h>
#include <neutrino-generator/Steering.h>

#include <MuonGun/Cylinder.h>

#include <earthmodel-service/EarthModelService.h>

#include <weighting/fluxes.h>
#include <boost/range/irange.hpp>

double muon_energy_threshold(double max_range)
{
	double a=0.212/1.2, b=0.251e-3/1.2;
	return (std::exp(max_range*b) - 1)*a/b;
}

void hoboplopia(I3MCTree &tree, I3MuonGun::SamplingSurface &surface, I3RandomService &rng)
{
#if 1
	auto generator = GaisserH3a(6e2, 1e11);
	
	double t0 = tree.get_heads().at(0).GetTime();
	double window = 40*I3Units::microsecond;
	double rate = surface.GetAcceptance(0,1)*generator.total()*I3Units::hertz;
	int n = rng.Poisson(window*rate);
	
	for (auto i : boost::irange(0,n,1)) {
		I3Particle::ParticleType type;
		double energy;
		I3Position pos;
		I3Direction dir;
		std::tie(type, energy) = generator(rng);
		surface.SampleImpactRay(pos, dir, rng, 0, 1);

		I3Particle p;
		p.SetType(type);
		if (energy/nucleonNumber(type) < muon_energy_threshold(I3CORSIKAReaderUtils::GetSlantDepth(dir, pos)))
		{
			continue;
		}
		p.SetEnergy(energy);
		p.SetTime(t0 + rng.Uniform(-window/2, window/2));
		p.SetDir(dir);
		p.SetPos(pos);
		
		I3MCTreeUtils::AddPrimary(tree, p);
	}
	// log_info_stream("skipped " << n+1-I3MCTreeUtils::GetPrimaries(tree).size() << " of " << n << " coincident showers");
#endif
}

void run_single_event(const std::string &corsika_executable, int nshowers=1, bool neutrinos=true)
{
	auto rng = boost::make_shared<I3GSLRandomService>(0);

	// Accept InIce particles (neutrino secondaries), plus any muons that can
	// reach the sampling surface. Other particles are dropped.
	I3MuonGun::Cylinder surface(1600,800);
	auto inice_filter = [&surface](const I3Particle &p) {
		if (p.GetLocationType() == I3Particle::InIce)
			return true;
		else if (p.GetType() != I3Particle::MuMinus && p.GetType() != I3Particle::MuPlus)
			return false;
		double a=0.212/1.2, b=0.251e-3/1.2;
		double max_range = std::log(1 + p.GetEnergy()*b/a)/b;
		// i3_assert(surface.GetIntersection(p.GetPos(), p.GetDir()).first < max_range);
		return surface.GetIntersection(p.GetPos(), p.GetDir()).first < max_range;
	};
	
	// Set up CORSIKA
	boost::shared_ptr<CorsikaService> caas(new CorsikaService(corsika_executable));
	
	// Set up NeutrinoGenerator
	boost::shared_ptr<I3CosmicEventGenerator> generator;
	if (neutrinos) {
		auto earthmodel = boost::make_shared<earthmodel::EarthModelService>();
		auto steering = boost::make_shared<nugen::Steering>(earthmodel);
		auto interaction_info = boost::make_shared<I3NuGInteractionInfo>(rng, steering, "csms");
		interaction_info->Initialize();
		auto neutrino_propagator = boost::make_shared<I3NuGIncrementalEventGenerator>(rng, steering, interaction_info);
		auto neutrino_selector = boost::make_shared<I3NuGSelector>(steering, rng, 10);
	
		// Mix them all together in a big jug
		generator.reset(new I3CosmicEventGenerator(caas, inice_filter, neutrino_selector, neutrino_propagator));
	} else {
		generator.reset(new I3CosmicEventGenerator(caas, inice_filter));
	}

	I3MCTree source_tree;
	{
		I3Frame frame;
		I3MCTree tree;
	
		// Make a primary
		I3Particle p;
		p.SetType(I3Particle::PPlus);
		p.SetEnergy(1e2);
		p.SetTime(0);
		p.SetDir(I3Direction(70*I3Units::degree, 37*I3Units::degree));
		p.SetPos(I3Position(-8,7,3));
		I3MCTreeUtils::AddPrimary(tree, p);
	
		// Add coincident showers
		hoboplopia(tree, surface, *rng);
		source_tree = tree;
	}

	typedef I3Map<I3ParticleID, std::tuple<CorsikaClient::BiasParticleType, double, double>> ShowerBiasMap;

	typedef std::chrono::high_resolution_clock clock;
	auto start = clock::now();
	for (int i=0; i < nshowers; i++) {
		try {
#if 1
			I3Frame frame;
			I3MCTree tree;
		
			// Make a primary
			I3Particle p;
			p.SetType(I3Particle::PPlus);
			p.SetEnergy(1e2);
			p.SetTime(0);
			p.SetDir(I3Direction(70*I3Units::degree, 37*I3Units::degree));
			p.SetPos(I3Position(-8,7,3));
			I3MCTreeUtils::AddPrimary(tree, p);
		
			// Add coincident showers
			// hoboplopia(tree, surface, *rng);
			
			auto biases = boost::make_shared<ShowerBiasMap>();
			biases->emplace(p.GetID(), ShowerBiasMap::mapped_type(CorsikaClient::Mu,1e-3,1));
			frame.Put("ShowerBias", biases);
#else
			I3Frame frame;
			I3MCTree tree = source_tree;
#endif
			// Callback to feed particles down the simulation chain as they are 
			// produced. This could enqueue particles into CLSim, or simulate the
			// IceTop tank response (or both!)
			auto emitParticle = [](I3Particle &p) {
				// std::cout << p.GetTypeString() << " " << p.GetEnergy() << std::endl;
			};
			// Simulate showers from all primaries in the tree, be they cosmic rays
			// or neutrinos
			generator->Generate(tree, frame, emitParticle);
		
			for (auto &primary : tree.get_heads()) {
				auto node = tree.find(primary);
				if (I3MCTree::sibling_iterator(node) == tree.end_sibling())
					tree.erase(node);
			}
			// std::cout << tree << std::endl;
		
		} catch (std::exception &e) {
			log_error_stream(e.what());
		}
		if (i % std::max(nshowers/100, 10000) == 0) {
			auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now()-start);
			log_info_stream(i+1 << " of " << nshowers << ": " << elapsed.count()/double(i+1) << " ms/shower");
		}
	}
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now()-start);
	std::cout << elapsed.count()/double(nshowers) << " ms/shower" << std::endl;
}

// see: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
class welford_gap_timer {

public:
	welford_gap_timer(bool storeit=false) : t0(clock::now()), count(0), mean(0), mean2(0), store(storeit) {}
	
	void update()
	{
		auto t1 = clock::now();
		int value = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
		t0 = t1;
		count++;
		double delta = value - mean;
		mean += delta / count;
		double delta2 = value - mean;
		mean2 += delta*delta2;
		if (store)
			storage.push_back(value);
	}
	
	void report()
	{
		double variance = mean2/(count-1);
		if (count >= 2)
			log_info_stream(mean << " +- " << std::sqrt(variance) << " microseconds between events ("<<count<<")");
		
	}
	
	void dump(const std::string &fname)
	{
		std::fstream fs(fname.c_str(), std::ios::out | std::ios::binary);
		fs.write((char*)(storage.data()), sizeof(int)*storage.size());
		fs.close();
	}
	std::vector<int> storage;

private:
	typedef std::chrono::high_resolution_clock clock;
	decltype(clock::now()) t0;
	size_t count;
	double mean, mean2;
	bool store;
};

void run_powerlaw(const std::string &corsika_executable, int nshowers, double emin, double emax, double gamma, bool skip_subthreshold=false)
{
	auto rng = boost::make_shared<I3GSLRandomService>(0);

	auto spectrum = PowerLaw(1, emin, emax, gamma);

	// Accept InIce particles (neutrino secondaries), plus any muons that can
	// reach the sampling surface. Other particles are dropped.
	I3MuonGun::Cylinder surface(1600,800);
	auto inice_filter = [&surface](const I3Particle &p) {
		if (p.GetLocationType() == I3Particle::InIce)
			return true;
		else if (p.GetType() != I3Particle::MuMinus && p.GetType() != I3Particle::MuPlus)
			return false;
		double a=0.212/1.2, b=0.251e-3/1.2;
		double max_range = std::log(1 + p.GetEnergy()*b/a)/b;
		// i3_assert(surface.GetIntersection(p.GetPos(), p.GetDir()).first < max_range);
		return surface.GetIntersection(p.GetPos(), p.GetDir()).first < max_range;
	};
	
	// Set up CORSIKA
	boost::shared_ptr<CorsikaService> caas(new CorsikaService(corsika_executable));
	
	// Set up NeutrinoGenerator
	boost::shared_ptr<I3CosmicEventGenerator> generator;
	generator.reset(new I3CosmicEventGenerator(caas, inice_filter));

	typedef std::chrono::high_resolution_clock clock;
	auto start = clock::now();
	welford_gap_timer gap_timer(true);
	welford_gap_timer event_timer(true);
	
	for (int i=0; i < nshowers; i++) {
		try {
			I3Frame frame;
			I3MCTree tree;
		
			// Make a primary
			I3Particle p;
			p.SetType(I3Particle::PPlus);

			p.SetEnergy(spectrum(*rng));
			p.SetTime(0);
			I3Direction dir;
			I3Position pos;
			surface.SampleImpactRay(pos, dir, *rng, 0, 1);
			p.SetDir(dir);
			p.SetPos(pos);
			if (p.GetEnergy()/nucleonNumber(p.GetType()) < muon_energy_threshold(I3CORSIKAReaderUtils::GetSlantDepth(dir, pos)))
			{
				if (skip_subthreshold)
					continue;
			}
			I3MCTreeUtils::AddPrimary(tree, p);
		
			// Callback to feed particles down the simulation chain as they are 
			// produced. This could enqueue particles into CLSim, or simulate the
			// IceTop tank response (or both!)
			auto emitParticle = [&gap_timer](I3Particle &p) {
				if (std::abs(p.GetType()) == std::abs(I3Particle::MuMinus))
					gap_timer.update();
				// std::cout << p.GetTypeString() << " " << p.GetEnergy() << std::endl;
			};
			// Simulate showers from all primaries in the tree, be they cosmic rays
			// or neutrinos
			generator->Generate(tree, frame, emitParticle);
		
			for (auto &primary : tree.get_heads()) {
				auto node = tree.find(primary);
				if (I3MCTree::sibling_iterator(node) == tree.end_sibling())
					tree.erase(node);
			}
			event_timer.update();
			// std::cout << tree << std::endl;
		
		} catch (std::exception &e) {
			log_error_stream(e.what());
		}
		if (i % std::max(nshowers/20, 10) == 0) {
			auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now()-start);
			log_info_stream(i+1 << " of " << nshowers << ": " << elapsed.count()/double(i+1) << " ms/shower");
			gap_timer.report();
		}
	}
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now()-start);
	std::cout << elapsed.count()/double(nshowers) << " ms/shower" << std::endl;
	gap_timer.report();
	gap_timer.dump("corsika-remotecontrol-gaps.dat");
	event_timer.dump("corsika-remotecontrol-event-gaps.dat");
	
}


void run_single_shower(const std::string &corsika_executable)
{
	CorsikaClient::config_type config;
	config.rng[0].seed = 42;
	config.rng[0].ncalls = 17;
	config.rng[1].seed = 538;
	config.atmosphere = 11;
	// config.elcuts[2] = 1;
	// config.elcuts[3] = 1;
	CorsikaClient caas(corsika_executable, config);

	typedef std::chrono::high_resolution_clock clock;
	auto start = clock::now();
	int nshowers = 10000;
	double total_weight = 0;
	for (int i=0; i < nshowers; i++) {
		std::array<double,4> elcuts = {{320.,320.,1e21,1e21}};
		caas.StartShower(14, 1e6, 70*M_PI/180, 0, CorsikaClient::NuE, 1e-3, elcuts);
		bool finished = false;
		int muons = 0;
		while (true) {
			std::vector<double> particle = caas.NextParticle();
			if (particle.empty())
				break;
		
			int ptype = particle[0];
			double gamma = particle[1];
			double ct = particle[2];
			double height = particle[5];
			double x = particle[7];
			double y = particle[8];
		
#if 0
			std::cout << "obslevel  "
				<< " type " << ptype
				<< " gamma " << gamma
				<< " ct " << ct
				<< " height " << height
				<< " x " << x
				<< " y " << y
				<< " fields " << particle.size()
			
				<< std::endl;
#endif
		
			if (ptype == 5 || ptype == 6) {
				muons++;
				if (muons % 100 == 0) {
					std::cout << "--- " << muons << std::endl;
				}
			}

		}
		auto header = caas.GetEventHeader();
		auto trailer = caas.GetEventEnd();
		if (i % 100 == 0)
			std::cout << i << ": " << muons << " muons, bias: " << header[219] << ", weight: " << trailer[266] << ", max_x: " << trailer[267] << std::endl;
		total_weight += trailer[266];
	}
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now()-start);
	std::cout << "total weight: " << total_weight << " ("<<nshowers<<" showers)" << std::endl;
	std::cout << elapsed.count()/double(nshowers) << " ms/shower" << std::endl;
	
}

int main (int argc, char const *argv[])
{
	if (argc < 2) {
		std::cerr << "need at least 1 arguments (got "<<argc-1<<")" << std::endl;
		exit(1);
	}
	
	// run_powerlaw(argv[1], 10000,  1e6, 1e9, -2, true);
	// run_powerlaw(argv[1], 1000000,  6e2, 1e5, -2.6, true);
	// run_powerlaw(argv[1], 100000,  3e4, 1e6, -2.6, true);
	
	// run_powerlaw(argv[1], 10,  1e4, 1e5, -2, true);
	
	run_single_event(argv[1], 100, false);
	// run_single_shower(argv[1]);
	
	/* code */
	return 0;
}
