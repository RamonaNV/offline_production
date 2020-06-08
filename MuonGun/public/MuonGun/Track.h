/** $Id: Track.h 128654 2015-02-04 18:34:51Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 128654 $
 * $Date: 2015-02-04 11:34:51 -0700 (Wed, 04 Feb 2015) $
 */

#ifndef I3MUONGUN_TRACK_H_INCLUDED
#define I3MUONGUN_TRACK_H_INCLUDED

#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTree.h>
#include <boost/tuple/tuple.hpp>

I3_FORWARD_DECLARATION(I3MMCTrack);
typedef I3Vector<I3MMCTrack> I3MMCTrackList;

namespace I3MuonGun {

/**
 * @brief A particle that moves in a straight line, and loses energy as it does so.
 */
class Track : public I3Particle {
public:
	Track() {};
	Track(const I3MMCTrack &mmctrack,
	    const I3MCTree::sibling_const_iterator &secondaries_begin,
	    const I3MCTree::sibling_const_iterator &secondaries_end);
	
	/**
	 * Get the energy of the particle at the given down-track distance,
	 * assuming that the "continuous" contribution to the energy loss
	 * is constant between measurement points.
	 *
	 * @param[in] length distance from the track origin
	 * @returns an energy if 0 <= length < range; otherwise 0
	 */
	double GetEnergy(double length) const;
	I3Position GetPos(double length) const;
	double GetTime(double length) const;
	
	// Un-hide overridden base class methods
	using I3Particle::GetEnergy;
	using I3Particle::GetPos;
	using I3Particle::GetTime;
	
	/**
	 * @brief Extract energy losses from frame objects.
	 *
	 * Find the stochastic energy losses associated with each I3MMCTrack
	 * in the I3MCTree, and store them together in a Track.
	 */
	static std::list<Track> Harvest(const I3MCTree &, const I3MMCTrackList &);

	/**
	 * @brief A point at which the absolute energy of the particle is known
	 */
	struct Checkpoint {
		double length, energy;
		ptrdiff_t offset;
		Checkpoint(double l, double e=0., size_t o=0)
		    : length(l), energy(e), offset(o) {}
	};
	/**
	 * @brief The sum of stochastic energy losses since the last checkpoint
	 */
	struct LossSum {
		double length, energy;
		LossSum(double l, double e=0.)
		    : length(l), energy(e) {}
	};
	
	std::vector<Checkpoint> GetCheckpoints() const { return checkpoints_; }
	std::vector<LossSum> GetLosses() const { return losses_; }
	
private:
	std::vector<Checkpoint> checkpoints_;
	std::vector<LossSum> losses_;
};

I3_POINTER_TYPEDEFS(Track);

}

#endif
