#include "simclasses/I3ShowerBias.h"
#include "simclasses/I3CorsikaWeight.h"
#include "corsika-reader/I3CORSIKAService.h"
#include "corsika-reader/I3CORSIKAReaderUtils.h"

double muon_energy_threshold(double max_range)
{
	double a=0.212/1.2, b=0.251e-3/1.2;
	return (std::exp(max_range*b) - 1)*a/b;
}

void CorsikaService::StartShower(I3Particle &primary, const I3Frame &frame)
{
	using I3CORSIKAReaderUtils::EarthRadius;
	using I3Constants::OriginElev;
	
	i3_assert(std::isfinite(primary.GetTime()));
	i3_assert(std::isfinite(primary.GetEnergy()));
	i3_assert(std::isfinite(primary.GetDir().GetZenith()));
	i3_assert(std::isfinite(primary.GetDir().GetAzimuth()));
	i3_assert(std::isfinite(primary.GetPos().GetX()));
	i3_assert(std::isfinite(primary.GetPos().GetY()));
	i3_assert(std::isfinite(primary.GetPos().GetZ()));
	
	// Store the desired shower axis and time
	I3Position target_pos = primary.GetPos();
	I3Direction target_dir = primary.GetDir();
	double target_time = primary.GetTime();
	
	// Find the position where the shower core will hit the surface
	double overburden = I3CORSIKAReaderUtils::GetSlantDepth(target_dir, target_pos);
	I3Position core = target_pos - overburden*target_dir;
	// The local zenith is a line from the core impact position to the
	// center of the Earth. This is used to convert between the IceCube
	// coordinate frame (centered at 0,0, 1948.07 m below the surface) and
	// the CORSIKA coordinate frame (centered on the shower axis at the
	// surface).
	I3Direction local_zenith(core.GetX(),core.GetY(),core.GetZ()+EarthRadius+OriginElev);
	
	// Set energy threshold such that muons that will not reach the depth of
	// the requested core position are dropped.
	// FIXME: hadron energy threshold needs to be independent of overburden if
	// neutrino biasing is requested
	std::array<double,4> elcuts = {{
	    muon_energy_threshold(overburden),
	    muon_energy_threshold(overburden),
	    1e21,
	    1e21
	}};

        save_weight_=false;
	I3ShowerBias bias(I3ShowerBias::Mu,1.0);
	auto biases = frame.Get<I3ShowerBiasMapConstPtr>();
	if (biases) {
		auto it = biases->find(primary.GetID());
		if (it != biases->end()){
                  save_weight_=true;                  
                  bias = it->second;
                }
	}
	// Some conversions between IceTray and CORSIKA conventions:
	// - Particle encoding: PDG -> CORSIKA
	// - Zenith angle: relative to z axis -> relative to zenith angle at
	//   the center of the observation level (i.e. the shower axis)
	// - Azimuth angle: origin relative to x axis -> direction relative to
	//   magnetic north
	std::vector<double> block = CorsikaClient::StartShower(
	    I3CORSIKAReaderUtils::PDGToCorsika(primary.GetPdgEncoding()),
	    primary.GetEnergy(),
	    std::acos(-local_zenith*primary.GetDir()),
	    (-primary.GetDir()).GetAzimuth() - magnetic_north_,
            bias.type,bias.target,elcuts);
	
	std::vector<float> header = CorsikaClient::GetEventHeader();
	i3_assert(header.size() == 273);
	// Check that our assumptions are valid.
	// NB: indices here are zero based, i.e. one smaller than those in the
	//     CORSIKA Users Guide.
	if (header[78] == 0)
		log_fatal("CORSIKA must be built with the CURVED option");
	if (header[167] == 0)
		log_fatal("Observation level is not curved");
	if (header[92] != 0)
		log_fatal("Coordinate system is not aligned with magnetic north!");

	// Check that bias specification came back the same
        i3_assert(header[219] == float(bias.target));
        i3_assert(header[220] == float(bias.type));
	
	// Store core position as a rotation around the center of the Earth
	core_position_.SetOrientation(local_zenith, I3Direction(0,1,0).Cross(local_zenith));
	
	// Reset time shift
	time_shift_ = 0;
	
	FillParticle(block, primary, false);
	primary.SetShape(I3Particle::Primary);
	{
		double interaction_height = -header[6]*I3Units::cm;
		if (interaction_height < 0)
			log_fatal("Tracking should start at the top of the atmosphere (TSTART = true)");
		primary.SetLength(-I3CORSIKAReaderUtils::GetSlantDepth(
		    primary.GetDir(), primary.GetPos(), interaction_height));
	}
	
	// Verify coordinate frame conversions

        if (std::acos(primary.GetDir()*target_dir) > 1e-1*I3Units::degree){
          log_warn("target direction %g deg from primary",std::acos(primary.GetDir()*target_dir));
        }
	I3Position actual_core = primary.GetPos() - I3CORSIKAReaderUtils::GetSlantDepth(primary.GetDir(), primary.GetPos())*primary.GetDir();
	i3_assert((actual_core-core).Magnitude() < I3Units::m);
	
	// Shift times so that the shower front arrives at the target position
	// at the desired time
	time_shift_ = primary.GetTime() - target_time + (primary.GetPos()-target_pos).Magnitude()/primary.GetSpeed();
	primary.SetTime(primary.GetTime() - time_shift_);

	// save primary information so we can add it to the weight information in EndEvent()
	primary_=primary;
}

bool CorsikaService::NextParticle(I3Particle &particle)
{
	std::vector<double> block = CorsikaClient::NextParticle();
	if (block.empty())
		return false;
	
	FillParticle(block, particle, true);
	particle.SetShape(I3Particle::StartingTrack);
	particle.SetLocationType(I3Particle::IceTop);

	return true;
}

void CorsikaService::EndEvent(I3Frame &frame)
{

  if (save_weight_){
    auto header = GetEventHeader();
    auto trailer = GetEventEnd();
    I3ShowerBias bias(I3ShowerBias::BiasParticleType(header[220]),
                      header[219]);
    I3CorsikaWeightPtr weight_obj=boost::make_shared<I3CorsikaWeight>();
    weight_obj->primary=primary_;
    weight_obj->bias = bias;
    weight_obj->weight = trailer[266];
    weight_obj->max_x = trailer[267];
    
    frame.Put("I3CorsikaWeight",weight_obj);
  }
}
	
void CorsikaService::FillParticle(const std::vector<double> &block, I3Particle &particle, bool use_elevation)
{
	particle.SetPdgEncoding(I3CORSIKAReaderUtils::CorsikaToPDG(block[0]));
	double mass = particle.GetMass()/I3Units::GeV;
	if (mass > 0) {
		particle.SetEnergy(mass*block[1]);
		double beta = std::sqrt((1-1/block[1])*(1+1/block[1]));
		i3_assert(beta <= 1);
		particle.SetSpeed(beta*I3Constants::c);
	} else {
		particle.SetEnergy(block[1]*I3Units::GeV);
		particle.SetSpeed(I3Constants::c);
	}
	particle.SetTime(block[6]*I3Units::second - time_shift_);
	{
		// Convert local zenith angle to zenith angle relative to obs level
		double zenith = I3CORSIKAReaderUtils::LocalZenith(std::acos(block[2]), block[5]*I3Units::cm, I3Constants::SurfaceElev);
		double azimuth = std::atan2(-block[4], -block[3]) + magnetic_north_;
		// Rotate into a coordinate system where the shower core passes
		// through the target point
		particle.SetDir(core_position_.RotateOut(I3Direction(zenith, azimuth)));
	}
	{
		using I3CORSIKAReaderUtils::EarthRadius;
		double x = block[7]*I3Units::cm;
		double y = block[8]*I3Units::cm;
		// for primary, x and y are measured at sea level
		// for observation level, they're at the observation level
		double D = EarthRadius + (use_elevation ? block[5]*I3Units::cm : 0);
	
		double phi = std::atan2(y,x) + magnetic_north_;
		double theta = std::sqrt(std::pow(x,2)+std::pow(y,2))/D;
		if (!use_elevation)
			D += (block[5]*I3Units::cm);
		I3Position pos(D*std::cos(phi)*std::sin(theta),
		               D*std::sin(phi)*std::sin(theta),
		               D*std::cos(theta));
		// Rotate into a coordinate system where the shower core passes
		// through the target point
		pos = core_position_.RotateOut(pos);
		// Convert from Earth-centered to IceCube-centered frame
		pos.SetZ(pos.GetZ() - EarthRadius - I3Constants::OriginElev);
		particle.SetPos(pos);
	}

}

const double CorsikaService::magnetic_north_ = -119*I3Units::degree;

//I3_SERIALIZABLE(ShowerBiasMap);
