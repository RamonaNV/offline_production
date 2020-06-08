/** $Id: EnsembleSampler.h 160287 2018-01-05 20:41:31Z juancarlos $
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision: 160287 $
 * $Date: 2018-01-05 13:41:31 -0700 (Fri, 05 Jan 2018) $
 */

#ifndef MUONGUN_ENSEMBLESAMPLER_H_INCLUDED
#define MUONGUN_ENSEMBLESAMPLER_H_INCLUDED

#include <MuonGun/I3MuonGun.h>
#include <icetray/I3Logging.h>

#include <vector>
#include <cmath>
#include <boost/foreach.hpp>

#include <phys-services/I3RandomService.h>

namespace I3MuonGun {

template <typename T, size_t N>
std::ostream& operator<<(std::ostream &o, const boost::array<T,N> &vals) {
	o << "[ ";
	BOOST_FOREACH(const T& v, vals)
		o << v << " ";
	o << "]";
	return o;
}

/// @brief An affine invariant Markov chain Monte Carlo (MCMC) sampler.
///
/// Goodman & Weare, Ensemble Samplers With Affine Invariance
///   Comm. App. Math. Comp. Sci., Vol. 5 (2010), No. 1, 65–80
///
/// This implementation is a simplified C++ port of emcee's EnsembleSampler:
/// http://dan.iel.fm/emcee/current/
///
template <typename Signature>
class EnsembleSampler {
public:
	typedef typename detail::traits<Signature>::array_type array_type;
	
	struct sample {
		sample(const array_type &p, double logprob) : point(p), log_probability(logprob) {}
		array_type point;
		double log_probability;
	};
	
	EnsembleSampler(boost::function<Signature> log_posterior,
	    const std::vector<array_type>& initial_ensemble) : log_posterior_(log_posterior), stretch_scale_(2.) {
		// Check each dimension for >1 unique value. Dimensions with only one
		// unique value reduce the effective dimensionality of the ensemble.
		unsigned effective_dimensions = 0;
		for (unsigned dim = 0; dim < detail::traits<Signature>::arity; dim++) {
			for (unsigned i = 1; i<initial_ensemble.size(); i++) {
				if (initial_ensemble[i][dim] != initial_ensemble[0][dim]) {
					effective_dimensions++;
					break;
				}
			}
		}
		if (effective_dimensions == 0) {
			log_fatal("Initial ensemble has only one unique point. Can't use this to propose a stretch move.");
		} 
		assert(effective_dimensions <= detail::traits<Signature>::arity);
		// See Goodman & Weare, Eq 9, 3rd line
		dimension_scale_ = effective_dimensions - 1.;
		if (initial_ensemble.size() % 2 != 0 || initial_ensemble.size() < 2*effective_dimensions) {
			log_fatal("Ensemble must be at least twice the dimensionality of the sampled space (and even)");
		}
		BOOST_FOREACH(const array_type &point, initial_ensemble) {
			ensemble_.push_back(sample(point, LogProbability(point)));
			if (!std::isfinite(ensemble_.back().log_probability))
				log_fatal("Initial ensemble point has non-finite probability");
		}
		half_size_ = ensemble_.size()/2;
		Reset();
	}
	
	void Reset() {
		total_samples_ = 0;
		accepted_samples_ = 0;
	}
	
	const std::vector<sample>& Sample(I3RandomService &rng) {
		// Update one half of the ensemble based on the positions of the other
		// half, then vice versa. These sub-loops can each be parallelized.
		for (unsigned i=0; i < half_size_; i++) {
			this->ProposeStretch(rng, i, half_size_);
		}
		for (unsigned i=0; i < half_size_; i++) {
			this->ProposeStretch(rng, i + half_size_, 0);
		}
		
		return ensemble_;
	}
	
	double GetAcceptanceRate() const {
		return double(accepted_samples_)/double(total_samples_);
	}
private:
	boost::function<Signature> log_posterior_;
	std::vector<sample> ensemble_;
	double stretch_scale_;
	double dimension_scale_;
	unsigned half_size_;
	unsigned total_samples_;
	unsigned accepted_samples_;
	
	double LogProbability(const array_type &point) {
		double value = detail::call(&log_posterior_, &point.front());
		return value;
	}
	
	// Step one point in the ensemble toward or away from a random point in the
	// other half of the ensemble
	void ProposeStretch(I3RandomService &rng, unsigned pos, unsigned offset) {
		const sample &p0 = ensemble_.at(pos);
		const sample &p1 = ensemble_.at(offset + rng.Integer(half_size_));
		
		double z = this->Stretch(rng);
		array_type q;
		for (unsigned i=0; i < detail::traits<Signature>::arity; i++) {
			q[i] = p1.point[i] - z*(p1.point[i] - p0.point[i]);
		}
		double log_probability = this->LogProbability(q);
		double log_ratio = dimension_scale_*std::log(z)
		    + log_probability - p0.log_probability;
		if (log_ratio > std::log(rng.Uniform())) {
			// accept new point
			ensemble_.at(pos) = sample(q, log_probability);
			accepted_samples_++;
		}
		total_samples_++;
	}
	
	double Stretch(I3RandomService &rng) const {
		return std::pow((stretch_scale_ - 1.)*rng.Uniform() + 1, 2)/stretch_scale_;
	}
};

}

#endif // MUONGUN_ENSEMBLESAMPLER_H_INCLUDED
