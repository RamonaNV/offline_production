
#ifndef SIM_SERVICES_I3COSMICEVENTGENERATOR_H_INCLUDED
#define SIM_SERVICES_I3COSMICEVENTGENERATOR_H_INCLUDED

#include <sim-services/I3PropagatorService.h>
#include <dataclasses/physics/I3MCTree.h>

class I3IncrementalEventGeneratorService {
public:
    virtual void StartShower(I3Particle &primary, const I3Frame &frame) = 0;
    virtual bool NextParticle(I3Particle &particle) = 0;
    virtual void EndEvent(I3Frame &);
    virtual ~I3IncrementalEventGeneratorService(){;}
};

I3_POINTER_TYPEDEFS(I3IncrementalEventGeneratorService);

class NeutrinoSelector {
public:
    virtual std::vector<I3Particle>::const_iterator Select(const std::vector<I3Particle> &neutrinos, I3Frame &frame) = 0;
    virtual ~NeutrinoSelector(){;}
};

I3_POINTER_TYPEDEFS(NeutrinoSelector);

class I3CosmicEventGenerator {
public:
    I3CosmicEventGenerator(I3IncrementalEventGeneratorServicePtr airShower,
        std::function<bool(const I3Particle &)> secondaryFilter=std::function<bool(const I3Particle &)>(),
        NeutrinoSelectorPtr neutrinoSelector=nullptr,
        I3IncrementalEventGeneratorServicePtr neutrinoPropagator=nullptr);
    I3CosmicEventGenerator(const I3CosmicEventGenerator&) = delete;
    void Generate(I3MCTree &, I3Frame &frame, std::function<void(I3Particle&)> emitParticle);
    
protected:
    I3IncrementalEventGeneratorServicePtr airShowerSimulator_;
    std::function<bool(const I3Particle &)> secondaryFilter_;
    NeutrinoSelectorPtr neutrinoSelector_;
    I3IncrementalEventGeneratorServicePtr neutrinoPropagator_;
};

#endif // SIM_SERVICES_I3COSMICEVENTGENERATOR_H_INCLUDED
