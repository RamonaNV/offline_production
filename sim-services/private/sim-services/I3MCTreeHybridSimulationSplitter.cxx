/**
 *  $Id: I3MCTreeHybridSimulationSplitter.cxx 108207 2013-07-13 15:19:02Z nwhitehorn $
 *  
 *  Copyright (C) 2012
 *  Claudio Kopper <ckopper@icecube.wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 */

#include <string>

#include <icetray/I3Module.h>
#include <icetray/I3Units.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTree.h>

#include <boost/foreach.hpp>

/**
 * @brief Splits an I3MCTree into two trees. One tree
 * will have all tracks set to shape "Dark", the other
 * one will have everything else set to "Dark".
 *
 * Can be used for hybrid simulations with two different
 * photon propagators, one for tracks and one for cascades.
 * (e.g. photonics/hit-maker and ppc).
 *
 * WARNING: Do not use this with current versions of PPC/IceTray!
 *          At the moment PPC will read *ALL* I3MCTrees from the
 *          frame and will simulate particles from all of them.
 *
 */
class I3MCTreeHybridSimulationSplitter : public I3Module
{
public:
    I3MCTreeHybridSimulationSplitter(const I3Context& context);
    void Configure();
    void DAQ(I3FramePtr frame);
    
private:
    bool ParticleTypeIsATrack(I3Particle::ParticleType type) const;
    bool ParticleShouldNotBeTouched(I3Particle::ParticleType type) const;

    std::string inputMCTreeName_;
    std::string outputMCTreeNameTracks_;
    std::string outputMCTreeNameCascades_;
    
	std::vector<I3Particle::ParticleType> lightProducingTracks_;
	std::vector<I3Particle::ParticleType> doNotTouchThese_;
};

I3_MODULE(I3MCTreeHybridSimulationSplitter);

I3MCTreeHybridSimulationSplitter::I3MCTreeHybridSimulationSplitter(const I3Context& context)
:
I3Module(context),
inputMCTreeName_("I3MCTree"),
outputMCTreeNameTracks_("I3MCTreeTracks"),
outputMCTreeNameCascades_("I3MCTreeCascades")
{
    AddParameter("InputMCTreeName",
                 "The input I3MCTree.",
                 inputMCTreeName_);

    AddParameter("OutputMCTreeNameTracks",
                 "A copy of the input MCTree with all cascades set to \"Dark\" will be stored here.",
                 outputMCTreeNameTracks_);

    AddParameter("OutputMCTreeNameCascades",
                 "A copy of the input MCTree with all tracks set to \"Dark\" will be stored here.",
                 outputMCTreeNameCascades_);

    // these will end up in the "Track" I3MCTree
    lightProducingTracks_.push_back(I3Particle::MuMinus);
    lightProducingTracks_.push_back(I3Particle::MuPlus);
    lightProducingTracks_.push_back(I3Particle::TauMinus);
    lightProducingTracks_.push_back(I3Particle::TauPlus);
    lightProducingTracks_.push_back(I3Particle::STauMinus);
    lightProducingTracks_.push_back(I3Particle::STauPlus);
    lightProducingTracks_.push_back(I3Particle::Monopole);

    // these will not be changed in either output tree
    doNotTouchThese_.push_back(I3Particle::NuE);
    doNotTouchThese_.push_back(I3Particle::NuEBar);
    doNotTouchThese_.push_back(I3Particle::NuMu);
    doNotTouchThese_.push_back(I3Particle::NuMuBar);
    doNotTouchThese_.push_back(I3Particle::NuTau);
    doNotTouchThese_.push_back(I3Particle::NuTauBar);

    AddOutBox("OutBox");
}

void
I3MCTreeHybridSimulationSplitter::Configure()
{
    GetParameter("InputMCTreeName", inputMCTreeName_);
    GetParameter("OutputMCTreeNameTracks", outputMCTreeNameTracks_);
    GetParameter("OutputMCTreeNameCascades", outputMCTreeNameCascades_);
}

bool I3MCTreeHybridSimulationSplitter::ParticleTypeIsATrack(I3Particle::ParticleType type) const
{
    BOOST_FOREACH(I3Particle::ParticleType t, lightProducingTracks_) {
        if (type == t) return true;
    }
    return false;
}

bool I3MCTreeHybridSimulationSplitter::ParticleShouldNotBeTouched(I3Particle::ParticleType type) const
{
    BOOST_FOREACH(I3Particle::ParticleType t, doNotTouchThese_) {
        if (type == t) return true;
    }
    return false;
}

void
I3MCTreeHybridSimulationSplitter::DAQ(I3FramePtr frame)
{
    I3MCTreeConstPtr inputMCTree = frame->Get<I3MCTreeConstPtr>(inputMCTreeName_);

    if (!inputMCTree)
        log_fatal("There is no I3MCTree named \"%s\" in this frame!",
                  inputMCTreeName_.c_str());

    // make two copies of the input tree
    I3MCTreePtr mcTreeTracks(new I3MCTree(*inputMCTree));
    I3MCTreePtr mcTreeCascades(new I3MCTree(*inputMCTree));
    
    // set all cascades to "Dark" in the "track" tree
    for (I3MCTree::iterator it = mcTreeTracks->begin();
         it != mcTreeTracks->end();
         ++it)
    {
        if (ParticleShouldNotBeTouched(it->GetType()))
            continue;
        
        if (!ParticleTypeIsATrack(it->GetType()))
            it->SetShape(I3Particle::Dark);
    }

    // set all tracks to "Dark" in the "cascade" tree
    for (I3MCTree::iterator it = mcTreeCascades->begin();
         it != mcTreeCascades->end();
         ++it)
    {
        if (ParticleShouldNotBeTouched(it->GetType()))
            continue;

        if (ParticleTypeIsATrack(it->GetType()))
            it->SetShape(I3Particle::Dark);
    }
    
    // store the results
    frame->Put(outputMCTreeNameTracks_, mcTreeTracks);
    frame->Put(outputMCTreeNameCascades_, mcTreeCascades);
    
    PushFrame(frame);
}
