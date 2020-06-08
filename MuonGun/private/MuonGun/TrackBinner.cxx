/** $Id: TrackBinner.cxx 128962 2015-02-12 06:57:51Z olivas $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 128962 $
 * $Date: 2015-02-11 23:57:51 -0700 (Wed, 11 Feb 2015) $
 */

#include <MuonGun/TrackBinner.h>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <icetray/I3Units.h>

namespace I3MuonGun {

TrackBinner::TrackBinner(double mindepth, double maxdepth, unsigned steps)
{	
	std::vector<double> multbins;
	{
		using namespace boost::assign;
		multbins += 0.5, 1.5, 2.5, 3.5, 9.5, 19.5, 39.5, 99.5;
	}
	
	using namespace histogram::binning;
	histogram::histogram<5>::bin_specification mspecs;
	mspecs[0] = uniform<cosine>::create(0, M_PI/2, 11);
	double dh = (maxdepth-mindepth)/(2*(steps-1));
	mspecs[1] = uniform<>::create(mindepth-dh, maxdepth-dh, steps);
	mspecs[2] = multbins;
	mspecs[3] = uniform<power<2> >::create(0, 250, 31);
	mspecs[4] = uniform<histogram::binning::log10>::create(1, 1e7, 101);
	
	{
		histogram::histogram<2>::bin_specification specs = {{mspecs[0], uniform<histogram::binning::log10>::create(1e2, 1e11, 101)}};
		primary_ = boost::make_shared<histogram::histogram<2> >(specs);
	}
	
	{
		histogram::histogram<3>::bin_specification specs = {{
		    uniform<cosine>::create(0, M_PI/2, 101), mspecs[1], uniform<>::create(0.5, 99.5, 100)}};
		multiplicity_ = boost::make_shared<histogram::histogram<3> >(specs);
	}
	
	{
		histogram::histogram<4>::bin_specification specs = {{mspecs[0], mspecs[1], mspecs[2],
			uniform<power<2> >::create(0, 250, 101)}};
		radius_ = boost::make_shared<histogram::histogram<4> >(specs);
	}
	
	energy_ = boost::make_shared<histogram::histogram<5> >(mspecs);
	
}

void TrackBinner::Consume(boost::shared_ptr<const TrackBundle> tracks,
    double energy, double zenith, double weight)
{
	{
		boost::array<double, 2> values = {{zenith, energy}};
		primary_->fill(values, weight);

	}
	for (std::map<double, std::vector<CompactTrack> >::const_iterator
	    i = tracks->begin(); i != tracks->end(); i++) {

		double mult = i->second.size();
		boost::array<double, 3> mult_values = {{zenith, i->first/I3Units::km, mult}};
		multiplicity_->fill(mult_values, weight);
		boost::array<double, 4> radius_values = {{zenith, i->first/I3Units::km, mult, 0}};
		boost::array<double, 5> energy_values = {{zenith, i->first/I3Units::km, mult, 0, 0}};
		BOOST_FOREACH(const CompactTrack &track, i->second) {
			radius_values[3] = energy_values[3] = track.GetRadius();
			energy_values[4] = track.GetEnergy();
			radius_->fill(radius_values, weight/mult);
			energy_->fill(energy_values, weight/mult);
		}
	}
}

NeutrinoBinner::NeutrinoBinner()
{
	using namespace histogram::binning;
	// Zenith angle, bundle energy, neutrino energy, parent type
	histogram::histogram<3>::bin_specification mspecs;
	mspecs[0] = uniform<cosine>::create(0, M_PI/2, 101);
	mspecs[1] = uniform<histogram::binning::log10>::create(1, 1e8, 101);
	mspecs[2] = uniform<histogram::binning::log10>::create(1e2, 1e8, 101);
	
	boost::array<I3Particle::ParticleType, 3> nutypes = {{I3Particle::NuE, I3Particle::NuMu, I3Particle::NuTau}};
	BOOST_FOREACH(I3Particle::ParticleType ptype, nutypes) {
		histmap::mapped_type &vec = histograms_[ptype];
		for (int i=0; i < 4; i++)
			vec.push_back(boost::make_shared<hist>(mspecs));
	}
}

void NeutrinoBinner::Consume(boost::shared_ptr<const TrackBundle> tracks,
    I3MCTreeConstPtr tree, I3MapI3ParticleIDDoubleConstPtr particle_weights, double flux_weight)

{
	if (tree->size() == 0)
		return;
	double zenith = tree->begin()->GetZenith();
	double total_energy = 0;
	if (tracks->size() > 0)
		BOOST_FOREACH(const CompactTrack &track, tracks->begin()->second)
			total_energy += track.GetEnergy();
	else
		total_energy = 0;
	
	boost::array<double, 3> values = {{zenith, total_energy, 0.}};
	for (I3MCTree::const_iterator it = tree->begin(); it != tree->end(); it++) {
		const I3Particle &p = *it;
		if (p.IsNeutrino()) {
			double weight = flux_weight;
			if (particle_weights != NULL) {
				I3MapI3ParticleIDDouble::const_iterator w = particle_weights->find(p.GetID());
				if (w != particle_weights->end())
					weight *= w->second;
			}
			values[2] = p.GetEnergy();
			
			// Find neutrino type
			histmap::iterator ptype = histograms_.find(std::abs(p.GetPdgEncoding()));
			if (ptype == histograms_.end())
				log_fatal("Unknown neutrino type %s!", p.GetTypeString().c_str());
			
			histmap::mapped_type::iterator target;
			
			// Find parent type
			I3MCTree::const_iterator parent = tree->parent(it);
			int pcode = parent != tree->end() ? std::abs(parent->GetPdgEncoding()) : 0;
			// Classify mesons by heaviest consituent quark
			if (pcode/1000 == 0 && pcode/100 > 0 /* mesons */) {
				// charm
				if (pcode/100 == 4)
					target = ptype->second.begin()+2;
				// strange (including K_long)
				else if (pcode/100 == 3 || pcode == 130)
					target = ptype->second.begin()+1;
				// all other mesons (almost entirely pions)
				else
					target = ptype->second.begin();
			// Put everything else in a big jug
			} else {
				target = ptype->second.begin()+3;
			}
			
			(*target)->fill(values, weight);
		}
	}
}


}
