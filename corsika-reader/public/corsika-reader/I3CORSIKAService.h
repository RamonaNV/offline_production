
#ifndef I3CORSIKASERVICE_H_INCLUDED
#define I3CORSIKASERVICE_H_INCLUDED

#include <sim-services/I3CosmicEventGenerator.h>
#include <dataclasses/I3Orientation.h>
#include "corsika-reader/CorsikaClient.h"

/// @brief IceTray client for RemoteControl CORSIKA
///
/// This client is stateful; see also CorsikaClient
class CorsikaService : public I3IncrementalEventGeneratorService, private CorsikaClient {
public:
	CorsikaService(const std::string &corsika_executable) :
          CorsikaClient(corsika_executable), save_weight_(false)
	{}

  virtual ~CorsikaService(){;}
	
	/// @brief Start simulation for an air shower
	///
	/// @param[inout] primary The primary particle for the air shower. The
	///                       particle should be initialized with the type,
	///                       total energy, and direction (in IceCube frame) of
	///                       the primary, as well as a position and time. On
	///                       return, contains the primary particle at the top
	///                       of the atmosphere.
	void StartShower(I3Particle &primary, const I3Frame &frame);
	
	/// @brief Fetch the next particle at the surface
	///
	/// @param[out] particle particle to fill
	/// @returns false if the shower is finished, true otherwise
	bool NextParticle(I3Particle &particle);
	
	/// @brief Emit weighting info
	void EndEvent(I3Frame &frame);

private:
	I3Orientation core_position_;
	double time_shift_;
	// Azimuth angle of magnetic North
	const static double magnetic_north_;
	
	void FillParticle(const std::vector<double> &block, I3Particle &particle, bool use_elevation);

	I3Particle primary_;
  bool save_weight_;
};

#endif // I3CORSIKASERVICE_H_INCLUDED
