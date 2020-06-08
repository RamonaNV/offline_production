
#include "sim-services/I3CosmicEventGenerator.h"
#include <boost/make_shared.hpp>

void
I3IncrementalEventGeneratorService::EndEvent(I3Frame &)
{}

I3CosmicEventGenerator::I3CosmicEventGenerator(I3IncrementalEventGeneratorServicePtr airShower,
    std::function<bool(const I3Particle &)> secondaryFilter,
    NeutrinoSelectorPtr neutrinoSelector,
    I3IncrementalEventGeneratorServicePtr neutrinoPropagator)
    : airShowerSimulator_(airShower), secondaryFilter_(secondaryFilter),
      neutrinoSelector_(neutrinoSelector), neutrinoPropagator_(neutrinoPropagator)
{}

void
I3CosmicEventGenerator::Generate(I3MCTree &mctree, I3Frame &frame, std::function<void(I3Particle&)> emitParticle)
{
    bool primary_neutrino = false;
    std::vector<I3Particle> secondary_neutrinos;
    boost::optional<I3ParticleID> biased_primary;
    
    bool biased = true;
    for (auto primary = I3MCTree::sibling_iterator(mctree.begin());
        primary != mctree.end_sibling(); primary++) {
        if (primary->IsNeutrino()) {
            if (primary_neutrino) {
                log_fatal("> 1 primary neutrino in the MCTree. Why?");
            } else {
                assert(neutrinoPropagator_);
                primary_neutrino = true;
                neutrinoPropagator_->StartShower(*primary, frame);
                while (true) {
                    I3Particle secondary;
                    if (!neutrinoPropagator_->NextParticle(secondary))
                        break;
                    if (secondaryFilter_(secondary)) {
                        mctree.append_child(primary, secondary);
                        emitParticle(*mctree.find(secondary));
                    }
                }
            }
        } else {
            // TODO: actually detect which primary is supposed to be biased
            
            airShowerSimulator_->StartShower(*primary, frame);
            
            if (biased)
                biased_primary = primary->GetID();
            while (true) {
                I3Particle secondary;
                if (!airShowerSimulator_->NextParticle(secondary))
                    break;
                if (secondary.IsNeutrino() && biased && !primary_neutrino) {
                    secondary_neutrinos.push_back(secondary);
                } else if (!secondaryFilter_ || secondaryFilter_(secondary)) {
                    mctree.append_child(primary, secondary);
                    emitParticle(*mctree.find(secondary));
                }
            }
            airShowerSimulator_->EndEvent(frame);
        }
        biased = false;
    }
    
    if (!secondary_neutrinos.empty() && neutrinoSelector_ && neutrinoPropagator_) {
        auto target = neutrinoSelector_->Select(secondary_neutrinos, frame);
        mctree.append_child(*biased_primary, *target);
        auto target_it = mctree.find(*target);
        neutrinoPropagator_->StartShower(*target_it, frame);
        while (true) {
            I3Particle secondary;
            if (!neutrinoPropagator_->NextParticle(secondary))
                break;
            if (secondaryFilter_(secondary)) {
                mctree.append_child(target_it, secondary);
                emitParticle(*mctree.find(secondary));
            }
        }
        neutrinoPropagator_->EndEvent(frame);
    }
}