/**
 * Copyright (c) 2013
 * Juan Carlos Diaz-Velez <juancarlos@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id: I3PrimaryPulseMapper.cxx 157457 2017-07-31 21:02:29Z cweaver $
 *
 * @file I3PrimaryPulseMapper.cxx
 * @version $Revision: 157457 $
 * @date $Date: 2017-07-31 15:02:29 -0600 (Mon, 31 Jul 2017) $
 * @author Juan Carlos Diaz-Velez
 */

#include "I3PrimaryPulseMapper.h"
#include <sim-services/MCPEMCPulseTools.hpp>

I3PrimaryPulseMapper::PrimaryIDMap
I3PrimaryPulseMapper::buildPrimaryMap(const I3MCTree& tree)
{
  PrimaryIDMap map;
  std::vector<I3Particle> primaries = tree.get_heads();
  auto end=tree.end();
  //collect the children of each primary
  for(I3ParticleID primary : primaries)
  {
    auto iter = tree.find(primary);
    //used to detect when we would jump to the next subtree
    auto next_primary=tree.next_sibling(iter);
    bool has_next=(next_primary!=end);
    while(iter!=end && !(has_next && iter==next_primary))
    {
      map.insert(std::make_pair(iter->GetID(),primary));
      iter++;
    }
  }
  return map;
}

I3PrimaryPulseMapper::I3PrimaryPulseMapper(const I3Context& ctx) :
  I3Module(ctx),
  inputParticleIDMapName_("I3MCPulseSeriesMapParticleIDMap"),
  outputParticleIDMapName_("I3MCPulseSeriesMapPrimaryIDMap"),
  mcTreeName_("I3MCTree")
{
  AddParameter("InputParticleIDMapname",
               "Name of I3ParticleIDMap to read",
               inputParticleIDMapName_);

  AddParameter("OutputParticleIDMapname",
               "Name of I3ParticleIDMap to write",
               outputParticleIDMapName_);

  AddParameter("I3MCTreeName",
               "Name of I3MCTree to get particle relations from",
               mcTreeName_);

  AddOutBox("OutBox");
}

void I3PrimaryPulseMapper::Configure()
{
    GetParameter("InputParticleIDMapname", inputParticleIDMapName_);
    GetParameter("OutputParticleIDMapname", outputParticleIDMapName_);
    GetParameter("I3MCTreeName", mcTreeName_);
}

void I3PrimaryPulseMapper::DAQ(I3FramePtr frame)
{
  if(!frame->Has(mcTreeName_))
    log_fatal_stream("I3MCTree '" << mcTreeName_ << "' doesn't exist in the frame.");
  I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>(mcTreeName_);

  if(!frame->Has(inputParticleIDMapName_))
    log_fatal_stream("ParticleIDMap'" << inputParticleIDMapName_ << "' doesn't exist in the frame.");
  I3ParticleIDMapConstPtr pIDMap = frame->Get<I3ParticleIDMapConstPtr>(inputParticleIDMapName_);

  PrimaryIDMap primary_map = buildPrimaryMap(*mctree);
  I3ParticleIDMapPtr primaryPIDMap(new I3ParticleIDMap());
  //iterate over DOMs
  for(const auto& dom_entry : *pIDMap)
  {
	ParticlePulseIndexMap primary_pulse_map;
    //iterate over particles which produced light on this DOM
    for(const auto& particle_entry : dom_entry.second)
    {
      I3ParticleID source=particle_entry.first; //the source particle
      auto primary_it=primary_map.find(source);
      if(primary_it!=primary_map.end()) //does it have a known primary?
        source=primary_it->second;
      //append these p.e. indices to any that we may already have for this primary
      primary_pulse_map[source].insert(primary_pulse_map[source].end(),
        particle_entry.second.begin(),particle_entry.second.end());
    }
    //sort and deduplicate the per-primary index lists
    for(auto& particle_entry : primary_pulse_map)
    {
      std::sort(particle_entry.second.begin(),particle_entry.second.end());
      std::vector<uint32_t> temp;
      std::unique_copy(particle_entry.second.begin(),particle_entry.second.end(),std::back_inserter(temp));
      particle_entry.second.swap(temp);
    }
    
    primaryPIDMap->insert(std::make_pair(dom_entry.first,primary_pulse_map));
  }
  frame->Put(outputParticleIDMapName_,primaryPIDMap);
  PushFrame(frame);
}

I3_MODULE(I3PrimaryPulseMapper);
