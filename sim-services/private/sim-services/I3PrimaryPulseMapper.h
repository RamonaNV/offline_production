#ifndef I3PRIMARYPULSEMAPPER_H
#define I3PRIMARYPULSEMAPPER_H

#include <icetray/I3ConditionalModule.h>
#include <dataclasses/physics/I3MCTree.h>

//forward declare test harness
template<typename ModuleType>
class SingleModuleTestSetup;

/**
 * \brief Converts mapping information describing which particles produced each
 *        MCPE to a mapping to primary particles which were the parents of the
 *        light emitting particles.
 */
class I3PrimaryPulseMapper: public I3Module {
public:
  I3PrimaryPulseMapper(const I3Context& ctx);
  ~I3PrimaryPulseMapper(){};
  
  void Configure();
  void DAQ(I3FramePtr frame);
  void Finish(){};
  
  using PrimaryIDMap = std::unordered_map<I3ParticleID,I3ParticleID,i3hash<I3ParticleID>>;
  ///Construct a map of particle IDs to those particles' parent primary particles
  static PrimaryIDMap buildPrimaryMap(const I3MCTree& tree);
private:
  std::string inputParticleIDMapName_;
  std::string outputParticleIDMapName_;
  std::string mcTreeName_;
  
  SET_LOGGER("I3PrimaryPulseMapper");
  friend class SingleModuleTestSetup<I3PrimaryPulseMapper>;
};

#endif //I3PRIMARYPULSEMAPPER_H
