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
 * $Id: I3MCPEtoI3MCHitConverter.cxx 137416 2015-09-09 20:15:44Z olivas $
 *
 * @file I3MCPEtoI3MCHitConverter.cxx
 * @version $Revision: 137416 $
 * @date $Date: 2015-09-09 14:15:44 -0600 (Wed, 09 Sep 2015) $
 * @author Claudio Kopper
 */

#include <string>
#include <vector>

#include "icetray/I3ConditionalModule.h"

#include "simclasses/I3MCPE.h"
#include "dataclasses/physics/I3MCHit.h"

void pe_to_hit(const I3MCPE&, I3MCHit&);

/**
 * @brief Converts an I3MCPESeriesMap to an I3MCHitSeriesMap.
 */
class I3MCPEtoI3MCHitConverter: public I3ConditionalModule
{
public:

  I3MCPEtoI3MCHitConverter(const I3Context& ctx);
  ~I3MCPEtoI3MCHitConverter(){};

  void Configure();
  void DAQ(I3FramePtr frame);
  void Finish(){};

private:
  std::string inputResponse_;
  std::string outputResponse_;

private:
  SET_LOGGER("I3MCPEtoI3MCHitSeriesConverter");
};


I3MCPEtoI3MCHitConverter::I3MCPEtoI3MCHitConverter(const I3Context& ctx) :
I3ConditionalModule(ctx),
inputResponse_("I3MCPESeriesMap"),
outputResponse_("I3MCHitSeriesMap")
{
  AddParameter("InputResponse",
               "Vector of names of the input response serieses to combine.",
               inputResponse_);

  AddParameter("OutputResponse",
               "Name of the output response series",
               outputResponse_);

  AddOutBox("OutBox");
}

void I3MCPEtoI3MCHitConverter::Configure()
{
  GetParameter("InputResponse", inputResponse_);
  GetParameter("OutputResponse", outputResponse_);
}

void I3MCPEtoI3MCHitConverter::DAQ(I3FramePtr frame)
{
  I3MCHitSeriesMapPtr output(new I3MCHitSeriesMap);

  I3MCPESeriesMapConstPtr input = frame->Get<I3MCPESeriesMapConstPtr>(inputResponse_);

  if(!input){
    log_fatal("Frame is missing an input response.");
  }

  for (I3MCPESeriesMap::const_iterator map_iter = input->begin();
       map_iter != input->end(); map_iter++)
    {
      const OMKey& omkey = map_iter->first;
      const I3MCPESeries &pe_series = map_iter->second;

      for (I3MCPESeries::const_iterator series_iter = pe_series.begin();
           series_iter != pe_series.end(); ++series_iter)
        {
          // this will insert an empty vector automatically in case there is none

          I3MCHitPtr hit(new I3MCHit());
          pe_to_hit(*series_iter,*hit);
          (*output)[omkey].push_back(*hit);
        }

    }

    frame->Put(outputResponse_, output);
    PushFrame(frame,"OutBox");
}

I3_MODULE(I3MCPEtoI3MCHitConverter);

void pe_to_hit(const I3MCPE& pe, I3MCHit& hit){
    hit.SetParticleID(pe.ID.majorID,pe.ID.minorID);
    hit.SetTime(pe.time);
    hit.SetNPE(pe.npe);
    hit.SetHitSource(I3MCHit::SPE);
}
