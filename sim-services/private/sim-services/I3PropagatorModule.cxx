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
 * $Id: I3PropagatorModule.cxx 176200 2019-10-17 12:37:12Z jvansanten $
 *
 * @file I3PropagatorModule.cxx
 * @version $Revision: 176200 $
 * @date $Date: 2019-10-17 06:37:12 -0600 (Thu, 17 Oct 2019) $
 * @author Claudio Kopper
 */

#include <deque>
#include <boost/foreach.hpp>

#include "icetray/I3ConditionalModule.h"

#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "phys-services/I3RandomService.h"

#include "sim-services/I3PropagatorService.h"


/**
 * @brief Propagates all particles found in an MCTree that have
 * configured I3PropagatorServices. If one service returns
 * another particle that can be propagated by another service,
 * it will be passed to that one. The results, in turn, will
 * also be propagated until there are no more particles
 * to propagate.
 *
 * The ususal use case is to pass an MCTree to PROPOSAL,
 * the results to CMC and the resulting muons back to PROPOSAL.
 */
class I3PropagatorModule : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3PropagatorModule(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    ~I3PropagatorModule();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    void Configure();
    
    /**
     * The module needs to process Physics frames
     */
    void DAQ(I3FramePtr frame);
    
private:
    
    /// Parameter: Name of the I3MCTree frame object. 
    std::string inputMCTreeName_;

    /// Parameter: Name of the output I3MCTree frame object. If identical to the
    /// input or empty, the input object will be replaced.
    std::string outputMCTreeName_;
    
    /// Parameter: map of particle type to a propagator that should propagate this type.
    I3ParticleTypePropagatorServiceMapPtr particleToPropagatorServiceMap_;

    /// Parameter: a random number generator service
    I3RandomServicePtr random_;

    /// Parameter: set this to a non-empty string to save the RNG state to the frame
    ///  before each propagation (per event). If it already exist, it is loaded instead
    ///  and the RNG state is being reset before propagation. This allows to re-generate
    ///  MCTrees using just the output of ucr/corsika-reader or neutrino-generator and
    ///  the RNG state.
    std::string rngStateName_;

    /// Parameter: if true, skips nodes that appear to have been propagated already.
    bool skipAlreadyPropagated_;
  
private:
    // default, assignment, and copy constructor declared private
    I3PropagatorModule();
    I3PropagatorModule(const I3PropagatorModule&);
    I3PropagatorModule& operator=(const I3PropagatorModule&);
    
    SET_LOGGER("I3PropagatorModule");
};


I3_MODULE(I3PropagatorModule);

I3PropagatorModule::I3PropagatorModule(const I3Context& context) 
: I3ConditionalModule(context)
{
    
    AddParameter("InputMCTreeName",
                 "Name of the I3MCTree frame object.",
                 "I3MCTree");

    AddParameter("OutputMCTreeName",
                 "Name of the output I3MCTree frame object. If identical to the\n"
                 "input or empty, the input object will be replaced.",
                 "I3MCTree");

    AddParameter("PropagatorServices",
                 "A dict mapping I3Particle::ParticleType to I3PropagatorService.",
                 particleToPropagatorServiceMap_);

    AddParameter("RandomService",
                 "A random number generator service.",
                 random_);

    AddParameter("RNGStateName",
                 "Name under which to store the state of the supplied RNG",
                 "");

    AddParameter("SkipAlreadyPropagated",
                 "Don't re-propagate particles (experimental)",
                 false);

    // add an outbox
    AddOutBox("OutBox");

}

I3PropagatorModule::~I3PropagatorModule()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3PropagatorModule::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("InputMCTreeName", inputMCTreeName_);
    GetParameter("OutputMCTreeName", outputMCTreeName_);
    GetParameter("PropagatorServices", particleToPropagatorServiceMap_);
    GetParameter("RandomService", random_);
    GetParameter("RNGStateName", rngStateName_);
    GetParameter("SkipAlreadyPropagated", skipAlreadyPropagated_);

    if (!random_)
        log_fatal("No random number generator service was configured. Please set the \"RandomService\" parameter");

    if (!particleToPropagatorServiceMap_)
        log_fatal("Please configure a set of propagators and particle types using the \"PropagatorServices\" parameter");

    if (inputMCTreeName_=="")
        log_fatal("The \"InputMCTreeName\" parameter must not be empty.");

    // configure the services with our random number generator
    for (I3ParticleTypePropagatorServiceMap::const_iterator it =
        particleToPropagatorServiceMap_->begin(); it != particleToPropagatorServiceMap_->end(); ++it)
    {
        const I3PropagatorServicePtr &propagator = it->second;

        if (!propagator)
            log_fatal("(null) propagator in propagator map is invalid");

        propagator->SetRandomNumberGenerator(random_);
    }
}


void I3PropagatorModule::DAQ(I3FramePtr frame)
{
    log_debug("%s", __PRETTY_FUNCTION__);
    
    if (rngStateName_.size() > 0) {
        I3FrameObjectConstPtr rngState = frame->Get<I3FrameObjectConstPtr>(rngStateName_);

        if (rngState) {
            // there is a state in the frame. Set the RNG to that state.
            random_->RestoreState(rngState);
        } else {
            // there is no state in the frame. Get the current one and save it.
            rngState = random_->GetState();

            frame->Put(rngStateName_, rngState);
        }
    }

    I3MCTreeConstPtr inputMCTree = frame->Get<I3MCTreeConstPtr>(inputMCTreeName_);
    if (!inputMCTree) {
        log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                  inputMCTreeName_.c_str());
    }

    // allocate the output I3MCTree
    I3MCTreePtr outputMCTree(new I3MCTree(*inputMCTree));
    
    // extra stuff to be filled into the frame after propagation
    I3PropagatorService::DiagnosticMapPtr protoFrame(new I3PropagatorService::DiagnosticMap);

    // Extract a list of particles to work on
    std::deque<std::pair<I3MCTree::iterator, I3PropagatorServicePtr> > particlesToPropagate;
    
    // build a map of MCTree::iterator to I3PropagatorService
    for(I3MCTree::iterator t_iter = outputMCTree->begin();
        t_iter != outputMCTree->end(); t_iter++)
    {
        I3ParticleTypePropagatorServiceMap::const_iterator it =
            particleToPropagatorServiceMap_->find(t_iter->GetType());
        // don't propagate particle types that are not configures
        if (it == particleToPropagatorServiceMap_->end()) continue;
	
	if(skipAlreadyPropagated_ && !std::isnan(t_iter->GetLength()))
	   continue;
	
        // it's something we know how to propagate. Add it to the list.
        particlesToPropagate.push_back(std::make_pair(t_iter, it->second));
    }
    
    log_debug("Going to propagate %zu particles", particlesToPropagate.size());

    // The main propagation loop.
    while (!particlesToPropagate.empty())
    {
        // retrieve the first entry
        const std::pair<I3MCTree::iterator, I3PropagatorServicePtr> &currentItem =
            particlesToPropagate.front();

        const I3MCTree::iterator &currentParticle_it = currentItem.first;
        I3PropagatorServicePtr currentPropagator = currentItem.second;

        // propagate it!
        const std::vector<I3Particle> children =
        currentPropagator->Propagate(*currentParticle_it, protoFrame, frame);

        // Insert each of the children into the tree. While at it,
        // check to see if any of them are on the list and should be propagated.
        BOOST_FOREACH(const I3Particle& child, children)
        {
            const I3MCTree::iterator child_it =
                outputMCTree->append_child(currentParticle_it, child);
            
            I3ParticleTypePropagatorServiceMap::const_iterator it =
                particleToPropagatorServiceMap_->find(child.GetType());
            // In looping as in public health, don't consume your own output.
            if (it != particleToPropagatorServiceMap_->end() && (it->second != currentPropagator || std::isnan(child_it->GetLength())))
                particlesToPropagate.push_back(std::make_pair(child_it, it->second));
        }
        particlesToPropagate.pop_front();
    }

    // store the output I3MCTree
    if ((outputMCTreeName_=="") || (outputMCTreeName_==inputMCTreeName_)) {
        frame->Delete(inputMCTreeName_);
        frame->Put(inputMCTreeName_, outputMCTree);
    } else {
        frame->Put(outputMCTreeName_, outputMCTree);
    }
    
    // store auxiliary information from the propagators
    BOOST_FOREACH(const I3PropagatorService::DiagnosticMap::value_type &pair,
        std::make_pair(protoFrame->begin(), protoFrame->end()))
    {
        if (frame->Has(pair.first))
            log_fatal_stream("Frame already contains a key '"<<pair.first<<"'");
        frame->Put(pair.first, pair.second);
    }
    
    // that's it!
    PushFrame(frame);
}
