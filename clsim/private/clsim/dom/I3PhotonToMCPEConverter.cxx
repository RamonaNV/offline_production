/**
 * Copyright (c) 2011, 2012
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
 * $Id: I3PhotonToMCPEConverter.cxx 179373 2020-03-10 21:25:48Z jvansanten $
 *
 * @file I3PhotonToMCPEConverter.cxx
 * @version $Revision: 179373 $
 * @date $Date: 2020-03-10 15:25:48 -0600 (Tue, 10 Mar 2020) $
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <algorithm>
#include <cmath>

#include "clsim/dom/I3PhotonToMCPEConverter.h"

#include <boost/foreach.hpp>

#include "simclasses/I3Photon.h"
#include "simclasses/I3CompressedPhoton.h"

#include "phys-services/I3SummaryService.h"

#include "dataclasses/geometry/I3OMGeo.h"
#include "dataclasses/geometry/I3ModuleGeo.h"

#include "simclasses/I3MCPE.h"
#include "dataclasses/physics/I3ParticleID.h"
#include <sim-services/MCPEMCPulseTools.hpp>

#include "dataclasses/I3Constants.h"

// The module
I3_MODULE(I3PhotonToMCPEConverter);

I3PhotonToMCPEConverter::I3PhotonToMCPEConverter(const I3Context& context) 
: I3ConditionalModule(context)
{
    AddParameter("MCPEGenerator",
                 "Instance of I3CLSimPhotonToMCPEConverter",
                 mcpeGenerator_);

    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object.",
                 inputPhotonSeriesMapName_);

    outputMCPESeriesMapName_="MCPESeriesMap";
    AddParameter("OutputMCPESeriesMapName",
                 "Name of the output I3MCPESeries frame object.",
                 outputMCPESeriesMapName_);

    mergeHits_=false;
    AddParameter("MergeHits",
                 "Merge photoelectrons which are very close in time in the output.",
                 mergeHits_);

    // add an outbox
    AddOutBox("OutBox");
    
    numGeneratedHits_=0;

}

I3PhotonToMCPEConverter::~I3PhotonToMCPEConverter()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3PhotonToMCPEConverter::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("MCPEGenerator", mcpeGenerator_);
    if (!mcpeGenerator_)
        log_fatal("You must set an MCPEGenerator");

    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputMCPESeriesMapName", outputMCPESeriesMapName_);

    GetParameter("MergeHits",mergeHits_);
}

namespace {

    template <typename Photon>
    void CheckSanity(const Photon &photon)
    {}
    
    template <>
    void CheckSanity<I3Photon>(const I3Photon &photon) {
        // sanity check for unscattered photons: is their direction ok
        // w.r.t. the vector from emission to detection?
        if (photon.GetNumScattered()==0)
        {
            I3Position pp(photon.GetPos() - photon.GetStartPos());
            const double ppl = pp.Magnitude();
            const double cosang = (pp*photon.GetDir())/ppl;
            
            if ((cosang < 0.9) && (ppl>1.*I3Units::m)) {
                log_fatal("unscattered photon direction is inconsistent: cos(ang)==%f, d=(%f,%f,%f), pp=(%f,%f,%f) pp_l=%f",
                          cosang,
                          photon.GetDir().GetX(), photon.GetDir().GetY(), photon.GetDir().GetZ(),
                          pp.GetX(), pp.GetY(), pp.GetZ(),
                          ppl
                          );
            }
        }
    }
}

template <typename PhotonMapType>
std::pair<I3MCPESeriesMapPtr,I3ParticleIDMapPtr>
I3PhotonToMCPEConverter::Convert(I3FramePtr frame)
{
    boost::shared_ptr<const PhotonMapType> inputPhotonSeriesMap = frame->Get<boost::shared_ptr<const PhotonMapType> >(inputPhotonSeriesMapName_);
    
    // allocate the output hitSeriesMap
    I3MCPESeriesMapPtr outputMCPESeriesMap(new I3MCPESeriesMap());
    I3ParticleIDMapPtr outputParentInfo;
    typedef std::map<OMKey, MCHitMerging::MCPEStream> MCPEStreamMap;
    boost::optional<MCPEStreamMap> mcpeStreams;
    if (mergeHits_) //only create this if it is needed
        mcpeStreams.emplace(MCPEStreamMap());

    for (const auto &it : *inputPhotonSeriesMap)
    {
        const ModuleKey &module_key = it.first;
        const auto &photons = it.second;

        for (const auto &photon : photons)
        {
            CheckSanity(photon);

            auto hit = mcpeGenerator_->Convert(module_key, photon);
            if (!hit)
                continue;

            // add a new hit
            if (mergeHits_) {
                (*mcpeStreams)[std::get<0>(*hit)].insert(std::get<1>(*hit));
            } else {
                (*outputMCPESeriesMap)[std::get<0>(*hit)].push_back(std::get<1>(*hit));
            }
            numGeneratedHits_++;
        }
    }

    if (mergeHits_) {
        outputParentInfo.reset(new I3ParticleIDMap());
        for (auto &pair : *mcpeStreams) {
            std::tie(
                (*outputMCPESeriesMap)[pair.first],
                (*outputParentInfo)[pair.first]
            ) = pair.second.extractMCPEsWithPIDInfo();
        }
    } else {
        for (auto &pair : *outputMCPESeriesMap) {
            std::sort(
                pair.second.begin(),
                pair.second.end(),
                [](const I3MCPE &a, const I3MCPE &b) { return a.time < b.time; }
            );
        }
    }
    return std::make_pair(outputMCPESeriesMap,outputParentInfo);
}

void I3PhotonToMCPEConverter::DAQ(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3MCPESeriesMapPtr outputMCPESeriesMap;
    I3ParticleIDMapPtr outputParentInfo;
    if (frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_)) {
        std::tie(outputMCPESeriesMap,outputParentInfo) = Convert<I3PhotonSeriesMap>(frame);
    } else if (frame->Get<I3CompressedPhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_)) {
        std::tie(outputMCPESeriesMap,outputParentInfo) = Convert<I3CompressedPhotonSeriesMap>(frame);
    } else {
        log_debug("Frame does not contain an I3PhotonSeriesMap named \"%s\".",
                  inputPhotonSeriesMapName_.c_str());
        
        // do nothing if there is no input data
        PushFrame(frame);
        return;
    }
    
    // store the output I3MCPESeriesMap
    frame->Put(outputMCPESeriesMapName_, outputMCPESeriesMap);
    if (outputParentInfo) {
        frame->Put(outputMCPESeriesMapName_+"ParticleIDMap", outputParentInfo);
    }
    
    // that's it!
    PushFrame(frame);
}

void I3PhotonToMCPEConverter::Finish()
{
    // add some summary information to a potential I3SummaryService
    I3MapStringDoublePtr summary = context_.Get<I3MapStringDoublePtr>("I3SummaryService");
    if (summary) {
        const std::string prefix = "I3PhotonToMCPEConverter_" + GetName() + "_";
        
        (*summary)[prefix+"NumGeneratedHits"] = numGeneratedHits_;
    }
    
}

I3CLSimPhotonToMCPEConverter::~I3CLSimPhotonToMCPEConverter() {}

I3CLSimPhotonToMCPEConverterForDOMs::I3CLSimPhotonToMCPEConverterForDOMs(I3RandomServicePtr random,
    boost::shared_ptr<const std::map<OMKey, I3CLSimFunctionConstPtr>> wavelengthAcceptance, I3CLSimFunctionConstPtr angularAcceptance)
    : randomService_(random), wavelengthAcceptance_(wavelengthAcceptance), angularAcceptance_(angularAcceptance)
{}

I3CLSimPhotonToMCPEConverterForDOMs::~I3CLSimPhotonToMCPEConverterForDOMs() {}

boost::optional<std::tuple<OMKey,I3MCPE>>
I3CLSimPhotonToMCPEConverterForDOMs::Convert(const ModuleKey &mkey, const I3CompressedPhoton &photon) const
{
    boost::optional<std::tuple<OMKey,I3MCPE>> hit;
    double hitProbability = photon.GetWeight();
    if (hitProbability < 0.) log_fatal("Photon with negative weight found.");
    if (hitProbability == 0.) return hit;
    
    // Only treat DOM-sized DOMs for now
    const double domRadius = 165.1*I3Units::mm;
    const double distFromDOMCenter = photon.GetPos().Magnitude();
    if (std::abs(distFromDOMCenter - domRadius) > 3*I3Units::cm) {
        log_fatal("distance not %fmm.. it is %fmm (diff=%gmm) (OMKey=(%i,%u) (photon @ pos=(%g,%g,%g)mm))",
                  domRadius/I3Units::mm,
                  distFromDOMCenter/I3Units::mm,
                  (distFromDOMCenter-domRadius)/I3Units::mm,
                  mkey.GetString(), mkey.GetOM(),
                  photon.GetPos().GetX()/I3Units::mm,
                  photon.GetPos().GetY()/I3Units::mm,
                  photon.GetPos().GetZ()/I3Units::mm
                  );
    }
    
    double photonCosAngle = photon.GetDir().GetZ();
    photonCosAngle = std::max(-1., std::min(1., photonCosAngle));
    
    OMKey omkey(mkey.GetString(),mkey.GetOM());
    auto domAcceptance = wavelengthAcceptance_->find(omkey);
    if (domAcceptance == wavelengthAcceptance_->end()) {
        log_debug_stream("No wavelength acceptance configured for "<<omkey);
        return hit;
    }
    hitProbability *= domAcceptance->second->GetValue(photon.GetWavelength());
    log_trace("After wlen acceptance: prob=%g (wlen acceptance is %f)",
             hitProbability, domAcceptance->second->GetValue(photon.GetWavelength()));

    hitProbability *= angularAcceptance_->GetValue(photonCosAngle);
    log_trace("After wlen&angular acceptance: prob=%g (angular acceptance is %f)",
              hitProbability, angularAcceptance_->GetValue(photonCosAngle));

    if (hitProbability > 1.+1e-4) {
        log_warn("hitProbability==%f > 1: your hit weights are too high. (hitProbability-1=%g)", hitProbability, hitProbability-1.);

        double hitProbability = photon.GetWeight();

        const double photonAngle = std::acos(photonCosAngle);
        log_warn("Photon (lambda=%fnm, angle=%fdeg, dist=%fm, module=(%d,%d)) has weight %g, 1/weight %g",
                 photon.GetWavelength()/I3Units::nanometer,
                 photonAngle/I3Units::deg,
                 distFromDOMCenter/I3Units::m,
                 mkey.GetString(),
                 mkey.GetOM(),
                 hitProbability,
                 1./hitProbability);

        hitProbability *= domAcceptance->second->GetValue(photon.GetWavelength());
        log_warn("After wlen acceptance: prob=%g (wlen acceptance is %f)",
                 hitProbability, domAcceptance->second->GetValue(photon.GetWavelength()));

        hitProbability *= angularAcceptance_->GetValue(photonCosAngle);
        log_warn("After wlen&angular acceptance: prob=%g (angular acceptance is %f)",
                  hitProbability, angularAcceptance_->GetValue(photonCosAngle));
        
        log_fatal("cannot continue.");
    }
    
    // does it survive?
    if (hitProbability <= randomService_->Uniform()) return hit;
    
    // FIXME: roll arrival time correction into kernel
    hit.emplace(omkey,I3MCPE(photon.GetParticleID(), 1, photon.GetTime()));
    return hit;
}
