/**
 *  $Id: CORSIKAResampler.cxx 150715 2016-10-12 20:23:40Z hpandya $
 *  
 *  Copyright (C) 2012
 *  The IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 */

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3Module.h>
#include <dataio/I3FileStager.h>

#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/I3Double.h>
#include <simclasses/I3CorsikaInfo.h>
#include <simclasses/I3CorsikaShowerInfo.h>
#include <phys-services/I3RandomService.h>
#include <phys-services/surfaces/Cylinder.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/assign.hpp>
#include <boost/python.hpp>
#include <boost/scoped_ptr.hpp>

#include <cmath>

const double EarthRadius = 637131500*I3Units::cm; // as in CORSIKA

typedef I3Map<I3ParticleID, double> I3ParticleIDMap;
typedef boost::shared_ptr<I3ParticleIDMap> I3ParticleIDMapPtr;

class CORSIKAResampler : public I3Module
{
  public:
	CORSIKAResampler(const I3Context& context);
	void Configure();
	void Simulation(I3FramePtr frame);
	void DAQ(I3FramePtr frame);

  private:
	SET_LOGGER("I3CorsikaReader");

	I3MCTreePtr ResampleEvent(const I3MCTree, const I3CorsikaShowerInfo);

	I3RandomServicePtr rng_;
	std::string weightdictname_; 
	std::string orig_mctree_name_; 
	std::string mctree_name_; 
	int oversampling_;
	I3Surfaces::Cylinder surface_;
	bool volumeCorr_;

};

I3_MODULE(CORSIKAResampler);

CORSIKAResampler::CORSIKAResampler(const I3Context& context)
    : I3Module(context), surface_(1600*I3Units::m, 800*I3Units::m)
{

	AddParameter("OverSampling", "number of times to oversample showers ",1);
	AddParameter("CylinderHeight", "Height of detector cylinder",
	    surface_.GetLength());
	AddParameter("CylinderRadius", "Radius of detector cylinder",
	    surface_.GetRadius());
	AddParameter("ZenithBias", "True if CORSIKA was compiled with the "
	    "VOLUMECORR option, False if the VOLUMEDET option was used. The "
	    "default zenith bias (proportional to sin(theta) for surface "
	    "detectors) is not supported.", true);
	AddParameter("WeightMapName", "Name of weights map to store in frame",
	    "CorsikaWeightMap");
	AddParameter("InputMCTreeName", "Name of mctree in frame", "I3MCTree_preSampling");
	AddParameter("OutputMCTreeName", "Name of mctree in frame", "I3MCTree");

	AddOutBox("OutBox");
}

void CORSIKAResampler::Configure()
{
	std::string filename;
	
	GetParameter("OverSampling", oversampling_);
	// Oversampling by reading input files multiple times
	double length, radius;
	GetParameter("CylinderHeight", length);
	GetParameter("CylinderRadius", radius);
	surface_ = I3Surfaces::Cylinder(length, radius);
	GetParameter("ZenithBias", volumeCorr_);
	GetParameter("WeightMapName", weightdictname_ );
	GetParameter("InputMCTreeName", orig_mctree_name_); 
	GetParameter("OutputMCTreeName", mctree_name_); 

	rng_ = context_.Get<I3RandomServicePtr>();
	if (rng_ == NULL)
		log_fatal("No I3RandomService in context!");
	
}

void CORSIKAResampler::Simulation(I3FramePtr frame)
{
  const I3CorsikaInfo old_info = frame->Get<I3CorsikaInfo>();
  I3CorsikaInfoPtr new_info(new I3CorsikaInfo(old_info));
  new_info->cylinder_height = surface_.GetLength();  
  new_info->cylinder_radius = surface_.GetRadius();
  new_info->oversampling = oversampling_;
  frame->Delete("I3CorsikaInfo");
  frame->Put(new_info);
  log_info_stream (GetName() << " " << __PRETTY_FUNCTION__ << " " << *new_info);
  PushFrame(frame);
}

void CORSIKAResampler::DAQ(I3FramePtr frame)
{
	log_debug("%s: %s", GetName().c_str(), __PRETTY_FUNCTION__);
	log_trace("%s: %s", GetName().c_str(), __PRETTY_FUNCTION__);
	if (!frame->Has(orig_mctree_name_) )
		log_fatal("No tree %s found in frame!",orig_mctree_name_.c_str());
	if (frame->Has(mctree_name_) )
		log_fatal("Found %s: This event has already been sampled!",mctree_name_.c_str());
	const I3MCTree mctree = frame->Get<I3MCTree>(orig_mctree_name_);
	const I3MapStringDouble weightMap = frame->Get<I3MapStringDouble>(weightdictname_);
	const I3CorsikaShowerInfo showerInfo = frame->Get<I3CorsikaShowerInfo>("ShowerInfo");

	for ( int f_i = 0; f_i < oversampling_ ; f_i++)
	{ 
		I3FramePtr newframe(new I3Frame(I3Frame::DAQ));
		I3MCTreePtr o = I3MCTreePtr(new I3MCTree(mctree));
		newframe->Put(orig_mctree_name_,o);

		I3MCTreePtr t = ResampleEvent(mctree,showerInfo);
		const I3Particle primary = I3MCTreeUtils::GetPrimaries(*t)[0];
		newframe->Put(mctree_name_,t);

		// Update weights
		I3MapStringDoublePtr weightdict(new I3MapStringDouble(weightMap));

		double thetaMin = (*weightdict)["ThetaMin"];
		double thetaMax = (*weightdict)["ThetaMax"];

		(*weightdict)["OverSampling"] = (double)oversampling_;
		(*weightdict)["AreaSum"] = surface_.GetAcceptance(cos(thetaMax), cos(thetaMin));
		(*weightdict)["CylinderRadius"] = surface_.GetRadius();
		(*weightdict)["CylinderLength"] = surface_.GetLength();

		if (!volumeCorr_) { 
			(*weightdict)["AreaSum"] = 2*M_PI*
		    surface_.GetArea(primary.GetDir())
		    /(cos(thetaMin)-cos(thetaMax)); 
		}
		newframe->Put(weightdictname_, weightdict);

		I3DoublePtr height(new I3Double(showerInfo.firstIntHeight));
		newframe->Put("CorsikaInteractionHeight", height);

		PushFrame(newframe);
	}
}

// Distance to the surface
double
GetSlantDepth(const I3Direction &dir, const I3Position &pos)
{
	double d = I3Constants::SurfaceElev-I3Constants::OriginElev-pos.GetZ();
	double r = EarthRadius;
	double ct = cos(dir.GetZenith());
	return sqrt(2*r*d + ct*ct*(r-d)*(r-d) - d*d) - (r-d)*ct;
}



I3MCTreePtr
CORSIKAResampler::ResampleEvent(I3MCTree mctree,I3CorsikaShowerInfo shower_info)
{
	I3MCTreePtr newtree = I3MCTreePtr(new I3MCTree()); 
	const I3Particle orig_primary = I3MCTreeUtils::GetPrimaries(mctree)[0];
	I3Particle primary = orig_primary.Clone();

	// Randomly shift shower core around
	I3Position shower_intersection =
	    surface_.SampleImpactPosition(primary.GetDir(), *rng_);
	
	double xspeed, yspeed, zspeed;
	xspeed = primary.GetSpeed()*sin(primary.GetZenith() / I3Units::radian)*
	    cos(primary.GetAzimuth() / I3Units::radian);
	yspeed = primary.GetSpeed()*sin(primary.GetZenith() / I3Units::radian)*
	    sin(primary.GetAzimuth() / I3Units::radian);
	zspeed = primary.GetSpeed()*cos(primary.GetZenith() / I3Units::radian);


	// Time 0 is when the particle intersects the observation level
	double time_to_intersection;
	if (shower_info.curved) {
		double entry_altitude = EarthRadius +
		    shower_info.entryHeight;
		double detection_altitude = EarthRadius + 
			shower_info.obsLevelHeight + I3Constants::OriginElev;
		double detection_angle = primary.GetZenith();

		double flight_distance =
		    -detection_altitude*cos(detection_angle) +
		    sqrt(pow(entry_altitude, 2) -
		    pow(detection_altitude*sin(detection_angle), 2));
		primary.SetTime(-flight_distance/primary.GetSpeed());
		time_to_intersection = primary.GetTime() - GetSlantDepth(
		    primary.GetDir(), shower_intersection)/primary.GetSpeed();
	} else {
		primary.SetTime(-(shower_info.entryHeight -
		    I3Constants::SurfaceElev)/zspeed);
		time_to_intersection = primary.GetTime() -
		    (I3Constants::SurfaceElev - I3Constants::OriginElev
		    - shower_intersection.GetZ())/zspeed;
	}

	double time_to_interaction = 
		(shower_info.entryHeight - shower_info.firstIntHeight)/zspeed;
	primary.SetLength(primary.GetSpeed()*time_to_interaction);
	primary.SetPos(
             shower_intersection.GetX() - time_to_intersection*xspeed,
             shower_intersection.GetY() - time_to_intersection*yspeed,
             shower_intersection.GetZ() - time_to_intersection*zspeed
             ); 

        double distance = primary.GetTime()*primary.GetSpeed();
	I3Position surface_pos = primary.GetPos() - distance*primary.GetDir();

	I3MCTreeUtils::AddPrimary(*newtree, primary);

	//I3MCTree::iterator iter;
	std::vector<I3Particle>::iterator iter; 
	std::vector<I3Particle> daughters = I3MCTreeUtils::GetDaughters(mctree,orig_primary);
	std::vector<I3Particle> new_daughters;

	for (iter = daughters.begin(); iter != daughters.end();iter++) { 
		I3Particle part = iter->Clone();
		I3Position temp_pos = iter->GetPos(); 
		if (shower_info.curvedObs) {
			double theta, phi, D;
			// Normal particles are recorded at the observation
			// level, extended-history entries at production level
			D = EarthRadius + I3Constants::OriginElev + 
				shower_info.obsLevelHeight; 
			theta = hypot(
			  temp_pos.GetX() + surface_pos.GetX(),
			  temp_pos.GetY() + surface_pos.GetY()
			  )/D;
			phi = atan2(temp_pos.GetY() + surface_pos.GetY(),
			  temp_pos.GetX() + surface_pos.GetX());

			part.SetPos(
				    D*sin(theta)*cos(phi),
				    D*sin(theta)*sin(phi),
				    D*cos(theta) - EarthRadius -
				    I3Constants::OriginElev
				);
		} else {
			part.SetPos(temp_pos.GetX() + surface_pos.GetX(),
			  temp_pos.GetY() + surface_pos.GetY(),
			  shower_info.obsLevelHeight);
		}
		part.SetTime(iter->GetTime()+primary.GetTime());

		new_daughters.push_back(part);
	} 
	newtree->append_children(primary, new_daughters);

	return newtree;
}

