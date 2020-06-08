#include <sim-services/I3MCPEConverters.h>
#include <dataclasses/physics/I3MCHit.h>
#include <simclasses/I3MCPE.h>
#include <boost/foreach.hpp>

I3MCHit PEConversions::PEToHit(const I3MCPE& pe){
  I3MCHit hit(pe.ID.majorID, pe.ID.minorID);
  hit.SetTime(pe.time);
  hit.SetNPE(pe.npe);

  // PEs are only either SPEs or noise so we don't 
  // carry the enum HitSource.  If the I3ParticleID
  // is (0,0) then it's noise.  Otherwise it's an SPE
  I3MCHit::HitSource source = pe.ID.majorID == 0 && pe.ID.minorID == 0 ?
    I3MCHit::RANDOM : I3MCHit::SPE;
  hit.SetHitSource(source);
  return hit;
}

// not all I3MCHits can be converted to I3MCPEs
// only SPE and RANDOM can be converted
// the user has to check if the returned ptr is
// not NULL before dereferencing.
I3MCPEPtr PEConversions::HitToPE(const I3MCHit& hit){
  I3MCPEPtr pe;
  if(hit.GetHitSource() == I3MCHit::SPE ||
     hit.GetHitSource() == I3MCHit::RANDOM){
    uint64_t majorid(0);
    int32_t minorid(0);
    if(hit.GetHitSource() == I3MCHit::SPE){
      majorid = hit.GetParticleMajorID();
      minorid = hit.GetParticleMinorID();
    }

    pe = I3MCPEPtr(new I3MCPE);
    pe->ID.majorID = majorid;
    pe->ID.minorID = minorid;
    pe->time = hit.GetTime();
    pe->npe = hit.GetNPE();
  }else{
    log_debug("conversion from I3MCHit to I3MCPE failed.");
    log_debug("don't dereference the returned pointer.");
  }
  return pe;
}

void PEConversions::PEToHit(const I3MCPESeries& pes, 
                            I3MCHitSeries& mchits){
  BOOST_FOREACH(const I3MCPE& pe, pes) 
    mchits.push_back(PEToHit(pe));
}

void PEConversions::HitToPE(const I3MCHitSeries& mchits, 
                            I3MCPESeries& mcpes){
  BOOST_FOREACH(const I3MCHit& hit, mchits){
    I3MCPEPtr pe = HitToPE(hit);
    if(pe) mcpes.push_back(*pe);      
  }
}

void PEConversions::PEToHit(const I3MCPESeriesMap& pemap, 
                            I3MCHitSeriesMap& hsmap){
  BOOST_FOREACH(I3MCPESeriesMap::const_reference map_pair, pemap){
    I3MCHitSeries hits;
    PEToHit(map_pair.second, hits);
    hsmap[map_pair.first] = hits;
  }
}

void PEConversions::HitToPE(const I3MCHitSeriesMap& hsmap, 
                            I3MCPESeriesMap& pemap){
  BOOST_FOREACH(I3MCHitSeriesMap::const_reference map_pair, hsmap){
    I3MCPESeries pes;
    HitToPE(map_pair.second, pes);
    pemap[map_pair.first] = pes;
  }
}
