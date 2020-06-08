/** $Id: Muonitron.h 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 12:24:47 -0600 (Mon, 31 Aug 2015) $
 */

#ifndef MUONGUN_MUONITRON_H_INCLUDED
#define MUONGUN_MUONITRON_H_INCLUDED

#include <icetray/I3Module.h>
#include <dataclasses/I3Constants.h>
#include <dataclasses/physics/I3Particle.h>
#include <MuonGun/MuonPropagator.h>
#include <phys-services/surfaces/Sphere.h>

namespace {
	const double EarthRadius = 637131500*I3Units::cm; // as in CORSIKA
	const double SurfaceRadius = EarthRadius+I3Constants::SurfaceElev;
	const double OriginDepth = I3Constants::SurfaceElev-I3Constants::OriginElev;
	const double IceDensity = 0.917*1.005; // as in MMC
}

/**
 * @brief a utility module for recording the state of a muon bundle as it propagates to depth
 */
class Muonitron : public I3Module {
public:
	Muonitron(const I3Context &);
	void Configure();
	void DAQ(I3FramePtr);
		
	static double GetOverburden(double z, double d=OriginDepth, double r=SurfaceRadius);
	static double GetSurfaceZenith(double z, double d=OriginDepth, double r=SurfaceRadius);
	static double GetGeocentricZenith(double z, double d=OriginDepth, double r=SurfaceRadius);
	static double GetDetectorZenith(double z, double d=OriginDepth, double r=SurfaceRadius);
	
	static I3Direction RotateToZenith(const I3Direction&, const I3Direction&);
	static I3Position RotateToZenith(const I3Direction&, const I3Position&);
	static I3Particle RotateToZenith(const I3Particle&, const I3Particle&);
	static I3Position Impact(const I3Particle &);
private:
	bool PropagateTrack(I3Particle &target, double depth);
	
	boost::shared_ptr<I3MuonGun::MuonPropagator> bulkPropagator_;
	boost::shared_ptr<I3MuonGun::Crust> crust_;
	I3Surfaces::Sphere iceWorld_;
	
	std::vector<double> depths_;
};

#endif // MUONGUN_MUONITRON_H_INCLUDED
