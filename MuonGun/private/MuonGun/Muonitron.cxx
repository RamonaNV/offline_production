/** $Id: Muonitron.cxx 126481 2014-12-02 22:45:21Z david.schultz $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 126481 $
 * $Date: 2014-12-02 15:45:21 -0700 (Tue, 02 Dec 2014) $
 */

#include <MuonGun/Muonitron.h>
#include <MuonGun/CompactTrack.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <phys-services/I3Calculator.h>
#include <simclasses/I3MMCTrack.h>

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

I3_MODULE(Muonitron);

Muonitron::Muonitron(const I3Context &ctx) : I3Module(ctx), iceWorld_(I3Constants::SurfaceElev - I3Constants::OriginElev, 6374134)
{
	AddParameter("Depths", "Propagate muons to these vertical depths (in meters)", depths_);
	AddParameter("Crust", "Air and firn layers surrounding the bulk medium", crust_);
	AddParameter("Propagator", "MuonPropagator instance", bulkPropagator_);
	
	AddOutBox("OutBox");
}

void
Muonitron::Configure()
{
	GetParameter("Depths", depths_);
	GetParameter("Crust", crust_);
	GetParameter("Propagator", bulkPropagator_);
	
	if (depths_.size() == 0)
		log_fatal("You must specify at least one vertical depth!");
	if (!crust_)
		log_fatal("No crust layers configured!");
	if (!bulkPropagator_)
		log_fatal("No MMC propagator configured!");
}

// Convert to a coordinate system where the zenith is given by the given direction
// rather than I3Direction(0,0,-1)
I3Direction
Muonitron::RotateToZenith(const I3Direction &direction, const I3Direction &dir)
{
	I3Direction p(dir);
	p.RotateZ(-direction.GetAzimuth());
	p.RotateY(-direction.GetZenith());
	p.RotateZ(direction.GetAzimuth()); //get x-y orientation right
	return p;
}

I3Position
Muonitron::RotateToZenith(const I3Direction &direction, const I3Position &pos)
{
	I3Position p(pos);
	p.RotateZ(-direction.GetAzimuth());
	p.RotateY(-direction.GetZenith());
	p.RotateZ(direction.GetAzimuth()); //get x-y orientation right
	return p;
}

I3Position
Muonitron::Impact(const I3Particle &p)
{
	// Vector pointing from the origin to the anchor point,
	// projected along the track
	double l = (p.GetPos().GetX()*p.GetDir().GetX()
	          + p.GetPos().GetY()*p.GetDir().GetY()
	          + p.GetPos().GetZ()*p.GetDir().GetZ());
	// Shift anchor point down the track to the closest approach
	// to the origin
	return I3Position(p.GetPos().GetX() - l*p.GetDir().GetX(),
	                  p.GetPos().GetY() - l*p.GetDir().GetY(),
	                  p.GetPos().GetZ() - l*p.GetDir().GetZ());
}

I3Particle
Muonitron::RotateToZenith(const I3Particle &reference, const I3Particle &part)
{
	I3Particle p(part);
	p.SetDir(RotateToZenith(reference.GetDir(), p.GetDir()));
	// Force reference axis to pass through the origin
	I3Position impact = Impact(reference);
	I3Position anchor = I3Position(p.GetPos().GetX()-impact.GetX(),
	                               p.GetPos().GetY()-impact.GetY(),
	                               p.GetPos().GetZ()-impact.GetZ());
	p.SetPos(RotateToZenith(reference.GetDir(), anchor));
	return p;
}

// Find the distance to the surface from a point at depth d
double
Muonitron::GetOverburden(double zenith, double d, double r)
{
	double ct = cos(zenith);
	return sqrt(2*r*d + ct*ct*(r-d)*(r-d) - d*d) - (r-d)*ct;
}

// Transform a detector-centered zenith angle to an earth-centered zenith angle
double
Muonitron::GetGeocentricZenith(double zenith, double d, double r)
{
	double p = GetOverburden(zenith, d, r);
	return atan2(p*sin(zenith), p*cos(zenith) + (r-d));
}


double
Muonitron::GetSurfaceZenith(double zenith, double d, double r)
{
	return zenith - GetGeocentricZenith(zenith, d, r);
}

bool
Muonitron::PropagateTrack(I3Particle &target, double slant_depth)
{	
	double l = target.GetLength();
	target = bulkPropagator_->propagate(target, slant_depth);
	target.SetLength(l + std::min(slant_depth, target.GetLength()));
	return target.GetEnergy() > 0;
}

void
Muonitron::DAQ(I3FramePtr frame)
{
	I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>();
	if (!mctree)
		log_fatal("No MCTree!");
	
	I3MCTree::const_iterator it = mctree->begin();
	const I3Particle &primary = *it;
	
	std::list<I3Particle> tracks;
	for (it++; it != mctree->end(); it++)
		if (it->GetType() == I3Particle::MuPlus || it->GetType() == I3Particle::MuMinus) {
			// Transport from the atmosphere through the firn layer into the bulk ice
			I3Particle track = crust_->Ingest(*it);
			if (track.GetEnergy() > 0)
				tracks.push_back(track);
		}	
	
	I3MuonGun::TrackBundlePtr bundle = boost::make_shared<I3MuonGun::TrackBundle>();
	BOOST_FOREACH(double vdepth, depths_) {
		
		// Find the slant depth from the surface of the glacier to
		// a point vdepth meters below grid center.
		double dx = GetOverburden(primary.GetDir().GetZenith(), vdepth, 6374134);
		std::vector<I3MuonGun::CompactTrack> deep_tracks;
		for (std::list<I3Particle>::iterator pit = tracks.begin(); pit != tracks.end(); ) {
			
			// Propagate the muon the remaining distance to reach the desired depth.
			if (PropagateTrack(*pit, dx-pit->GetLength())) {
				deep_tracks.push_back(I3MuonGun::CompactTrack(RotateToZenith(primary, *pit)));
				pit++;
			} else {
				pit = tracks.erase(pit);
			}
		}
		
		(*bundle)[vdepth].swap(deep_tracks); 
	}
	
	frame->Put("MCPrimary", boost::make_shared<I3Particle>(primary));
	frame->Put("Tracks", bundle);
	PushFrame(frame);
}
