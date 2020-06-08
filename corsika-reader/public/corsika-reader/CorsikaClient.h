
#ifndef CORSIKACLIENT_H_INCLUDED
#define CORSIKACLIENT_H_INCLUDED

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <array>
#include "simclasses/I3ShowerBias.h"

class CorsikaClientImpl;

/// @brief Basic client for RemoteControl CORSIKA
///
/// This client is stateful, and requires calls in the following sequence:
/// 1. Call StartShower()
/// 2. Call NextParticle() until it returns an empty vector
/// 3. Call GetEventHeader()/GetEventTrailer()
///
class CorsikaClient {
public:
	struct config_type {
		int atmosphere;
		std::array<double,4> elcuts;
		struct rng_stream {
			uint32_t seed;
			uint64_t ncalls;
		};
		std::vector<rng_stream> rng;
		config_type();
	};
	
	/// param[in] corsika_executable path to CORSIKA binary
	CorsikaClient(const std::string &corsika_executable, config_type=config_type());
  virtual ~CorsikaClient();
	
	/// Start a new shower simulation. Blocks until an event header and primary
	/// are received from CORSIKA
	/// @param[in] particle_id Particle number in CORSIKA convention
	/// @param[in] energy      total energy of primary particle in GeV
	/// @param[in] theta       zenith angle in rad as measured from the
	///                        observation level
	/// @param[in] azimuth     azimuth angle in rad
	///
	/// @returns A particle block containing the primary at the top of the
	///          atmosphere.
	std::vector<double> StartShower(
          uint32_t particle_id, double energy, double theta, double phi,
          I3ShowerBias::BiasParticleType I3ShowerBiasBiasbias_target=I3ShowerBias::Mu,
          double bias_factor=1, std::array<double,4> elcuts={{0.3, 0.3, 0.003, 0.003}});
	
	/// Get the next particle at the observation level. Blocks until the next
	/// particle or the event end block is recieved.
	///
	/// May be called after StartShower().
	///
	/// @returns A particle block containing the particle at the lowest
	///          observation level, or an empty vector if there are no more
	///          particles.
	std::vector<double> NextParticle();
	
	/// Get the EVTH block from the current shower
	///
	/// May be called after StartShower()
	const std::vector<float>& GetEventHeader() const;
	/// Get the EVTE block from the current shower
	///
	/// May be called after NextParticle() has returned an empty vector
	const std::vector<float>& GetEventEnd() const;
	
private:
	std::unique_ptr<CorsikaClientImpl> impl_;

	bool more_particles_available_;
	
	std::vector<float> event_header_, event_trailer_;
	
};

#endif // CORSIKACLIENT_H_INCLUDED
