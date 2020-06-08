/**
 * Copyright (c) 2013
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
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
 * $Id: I3CombineMCPE.cxx 159050 2017-10-29 01:10:08Z cweaver $
 *
 * @file I3CombineMCPE.cxx
 * @version $Revision: 159050 $
 * @date $Date: 2017-10-28 19:10:08 -0600 (Sat, 28 Oct 2017) $
 * @author Claudio Kopper
 */

#include <sim-services/I3CombineMCPE.h>
#include "sim-services/MCPEMCPulseTools.hpp"

namespace{
  bool time_ordered(const I3MCPE& lhs, const I3MCPE& rhs){
    return (lhs.time < rhs.time);
  }
}

I3CombineMCPE::I3CombineMCPE(const I3Context& ctx) :
I3ConditionalModule(ctx),
outputResponse_("CombinedMCPEs")
{
  AddParameter("InputResponses",
               "Vector of names of the input response serieses to combine (no DEFAULT)",
               inputResponses_);

  AddParameter("OutputResponse",
               "Name of the output response series",
               outputResponse_);

  AddOutBox("OutBox");
}

void I3CombineMCPE::Configure()
{
    GetParameter("InputResponses", inputResponses_);
    GetParameter("OutputResponse", outputResponse_);

    if(inputResponses_.empty())
      log_fatal("No input responses set.");
}

void I3CombineMCPE::DAQ(I3FramePtr frame)
{
    I3MCPESeriesMapPtr output(new I3MCPESeriesMap);
    I3ParticleIDMapPtr outputInfo;

    for(const std::string& inputName : inputResponses_)
    {
        I3MCPESeriesMapConstPtr input = frame->Get<I3MCPESeriesMapConstPtr>(inputName);

        if(!input)
        {
          log_warn_stream("Frame is missing an input response " << inputName << ". Combining the ones I can find.");
            continue;
        }
        
        I3ParticleIDMapConstPtr inputInfo = frame->Get<I3ParticleIDMapConstPtr>(inputName+"ParticleIDMap");
        
        outputInfo = MergeMCPEs(output,input,0./*no time offset*/,outputInfo,inputInfo);
    }

    frame->Put(outputResponse_, output);
    PushFrame(frame,"OutBox");
}

I3_MODULE(I3CombineMCPE);

I3ParticleIDMapPtr
MergeMCPEs(I3MCPESeriesMapPtr output, I3MCPESeriesMapConstPtr input,
           float offsetTime,
           I3ParticleIDMapPtr outputInfo,
           I3ParticleIDMapConstPtr inputInfo)
{
    //Deal with all possible combinations of 'inline' vs. 'external' particle
    //info. Note that the four case below are not written exclusively (no
    //`else`s) because some cases are handled by touching up the data and
    //then falling into another case.
    if (!inputInfo && !outputInfo)
    { //easy case; no cumbersome external particle tracking
        for (const auto& om_entry : *input)
        {
            const OMKey& omkey = om_entry.first;
            I3MCPESeries& dest=(*output)[omkey];
            for (const I3MCPE& pe : om_entry.second)
            {
                dest.push_back(pe);
                dest.back().time+=offsetTime;
            }
        }
    }
    if (!inputInfo && outputInfo)
    { //need to 'upgrade' input to store particle information externally
        //A little ugly; we need to duplicate the input to make it mutable
        I3MCPESeriesMapPtr tmp(new I3MCPESeriesMap(*input));
        inputInfo = MCHitMerging::extractPIDInfo(*tmp);
        input = tmp;
    }
    if (inputInfo && !outputInfo)
    { //need to 'upgrade' output to store particle information externally
        outputInfo = MCHitMerging::extractPIDInfo(*output);
    }
    if (inputInfo && outputInfo)
    { //annoying case: need to merge the external tables
        for (const auto& om_entry : *input)
        {
            const OMKey& omkey = om_entry.first;
            //copy over the hits themselves
            I3MCPESeries& dest = (*output)[omkey];
            uint32_t countBefore=dest.size();
            for (const I3MCPE& pe : om_entry.second)
            {
                dest.push_back(pe);
                dest.back().time+=offsetTime;
            }
            //copy the parent information, updating the indices to account
            //for however much data was already there
            for (const auto& particle_entry : inputInfo->find(omkey)->second)
            {
                std::vector<uint32_t>& pdest = (*outputInfo)[omkey][particle_entry.first];
                for (uint32_t index : particle_entry.second)
                    pdest.push_back(index+countBefore);
            }
        }
    }
    
    // Make sure the resultant reponses are in time order
    for (auto& om_entry : *output)
    {
        I3MCPESeries::iterator beginning = om_entry.second.begin();
        I3MCPESeries::iterator ending = om_entry.second.end();
        if(!outputInfo)
            std::sort(beginning,ending,time_ordered);
        else
            MCHitMerging::sortMCHits(beginning,ending,(*outputInfo)[om_entry.first]);
    }
    
    return(outputInfo);
}
