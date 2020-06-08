
#include "clsim/I3CLSimLightSourcePropagatorFromI3PropagatorService.h"
#include "clsim/I3CLSimLightSource.h"

I3CLSimLightSourcePropagatorFromI3PropagatorService::I3CLSimLightSourcePropagatorFromI3PropagatorService(
    I3ParticleTypePropagatorServiceMapPtr particleToPropagatorServiceMap,
    bool trackParticleHistory,
    double cascadeBinWidth)
    : particleToPropagatorServiceMap_(particleToPropagatorServiceMap),
      initialized_(false),
      trackParticleHistory_(trackParticleHistory),
      cascadeBinWidth_(cascadeBinWidth/I3Constants::c)
{
    if (!particleToPropagatorServiceMap_)
        log_fatal("(null) propagator map is invalid");

    std::stringstream ss;
    for (auto &pair : *particleToPropagatorServiceMap_) {
        if (!pair.second)
            log_fatal("(null) propagator in map is invalid");
        I3Particle p;
        p.SetType(pair.first);

        std::string repr = boost::python::extract<const std::string>(boost::python::object(pair.second).attr("__repr__")());        
        ss << std::setw(20) << p.GetTypeString() << " => " << repr << std::endl;
    }
    log_info_stream("CLSim propagating particles with the following services:\n" << ss.str());
}

I3CLSimLightSourcePropagatorFromI3PropagatorService::~I3CLSimLightSourcePropagatorFromI3PropagatorService()
{}

void I3CLSimLightSourcePropagatorFromI3PropagatorService::SetRandomService(I3RandomServicePtr rng)
{
  for (auto &pair : *particleToPropagatorServiceMap_){
    pair.second->SetRandomNumberGenerator(rng);
  }
}

bool I3CLSimLightSourcePropagatorFromI3PropagatorService::IsValidForLightSource(const I3CLSimLightSource &source)
{
    // accept any particle in the type list that has not already been propagated
    return (source.GetType() == I3CLSimLightSource::Particle)
        && std::isnan(source.GetParticle().GetLength())
        && (particleToPropagatorServiceMap_->find(source.GetParticle().GetType()) != particleToPropagatorServiceMap_->end());
}

I3MCTreePtr I3CLSimLightSourcePropagatorFromI3PropagatorService::Convert(I3CLSimLightSourceConstPtr &lightSource, uint32_t identifier,
    secondary_callback emitSecondary, step_callback)
{
    std::deque<std::pair<I3Particle, I3PropagatorServicePtr>> queue;
    std::map<I3ParticleID, I3ParticleID> orphans;
    I3MCTreePtr history;
    if (trackParticleHistory_) {
        history.reset(new I3MCTree);
    }

    auto iter = particleToPropagatorServiceMap_->find(lightSource->GetParticle().GetType());
    assert( iter != particleToPropagatorServiceMap_->end());
    queue.emplace_back(lightSource->GetParticle(), iter->second);

    auto emitNonDark = [&](const I3Particle &particle) {
        if (particle.GetShape() != I3Particle::Dark) {
            I3CLSimLightSourceConstPtr secondary = boost::make_shared<I3CLSimLightSource>(particle);
            emitSecondary(secondary, identifier);
        }
    };
    
    while (!queue.empty()) {
        // retrieve the first entry
        auto &item = queue.front();

        // propagate it!
        auto secondaries = item.second->Propagate(item.first, NULL, NULL);
        // emit parent if the propagator didn't set its shape to dark
        emitNonDark(item.first);

        std::deque<I3Particle> coalesced_secondaries;
        for (auto &particle : secondaries) {
            
            auto iter = particleToPropagatorServiceMap_->find(particle.GetType());
            // return secondaries to propagation if they are not yet final
            if (iter != particleToPropagatorServiceMap_->end() &&
                (iter->second != item.second || std::isnan(particle.GetLength()))) {
                queue.emplace_back(particle, iter->second);
            } else {
                emitNonDark(particle);
            }
            if (trackParticleHistory_ &&
              !(
                // Omit particles that could be handled by the originating
                // propagator, but are marked as final, e.g. track slices
                (iter != particleToPropagatorServiceMap_->end()
                  && iter->second == item.second
                  && std::isfinite(particle.GetLength()))
                // Also omit neutrinos
                || particle.IsNeutrino()
              )) {
                if (coalesced_secondaries.empty() ||
                  !particle.IsCascade() ||
                  particle.GetType() == I3Particle::NuclInt ||
                  coalesced_secondaries.back().GetType() == I3Particle::NuclInt ||
                  particle.GetTime() - coalesced_secondaries.back().GetTime() > cascadeBinWidth_) {
                    coalesced_secondaries.push_back(particle);
                } else {
                    I3Particle &prev = coalesced_secondaries.back();
                    prev.SetEnergy(prev.GetEnergy()+particle.GetEnergy());
                    orphans.emplace(particle.GetID(),prev.GetID());
                }
            }
        }

        if (history) {
            auto parent = history->find(item.first);
            if (parent != history->end()) {
                *parent = item.first;
                for (auto &particle : coalesced_secondaries)
                    history->append_child(parent, particle);
            } else if (orphans.find(item.first.GetID()) == orphans.end()) {
                history->insert_after(item.first);
                parent = history->find(item.first);
                for (auto &particle : coalesced_secondaries)
                    history->append_child(parent, particle);
            }
            // NB: do not attach secondaries if parent was merged
        }

        queue.pop_front();
    }

    return history;
}

