#ifndef SIM_SERVICES_I3COMBINEMCPE_h
#define SIM_SERVICES_I3COMBINEMCPE_h

#include <icetray/I3ConditionalModule.h>
#include <simclasses/I3MCPE.h>
#include <simclasses/I3ParticleIDMap.hpp>

/**
 * @brief Combines several I3MCPEHitSeriesMaps into one.
 */
class I3CombineMCPE : public I3ConditionalModule
{
public:
  I3CombineMCPE(const I3Context& ctx);
  ~I3CombineMCPE(){};
  
  void Configure();
  void DAQ(I3FramePtr frame);
  void Finish(){};
  
private:
  std::vector<std::string> inputResponses_;
  std::string outputResponse_;
  
private:
  SET_LOGGER("I3CombineMCPE");
};

/**
  * MergeHits - Takes two I3MCPESeriesMaps and merges their hits.
  * The times of the hits in the second map are moved by offsetTime.
  *
  * @param map1 'original' map to which we will add the new hits
  * @param map2 second map hit times will be set within window of first map
  * @param offsetTime time difference between first hit in map2 relative to map1
  * @param map1Info any external particle parentage info for map1. 
  *                 Will be created if necessary.
  * @param map2Info any external particle parentage info for map2
  * @return an updated version of map1Info, or nothing if that object continues
  *         to be unnecessary.
  */
I3ParticleIDMapPtr
MergeMCPEs(I3MCPESeriesMapPtr map1, I3MCPESeriesMapConstPtr map2,
           float offsetTime,
           I3ParticleIDMapPtr map1Info,
           I3ParticleIDMapConstPtr map2Info);

#endif
