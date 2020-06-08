/**
 * copyright  (C) 2020 The Icecube Collaboration, All Rights Reserved
 */

#include "simclasses/converter/I3CorsikaWeightConverter.h"

I3TableRowDescriptionPtr I3CorsikaWeightConverter::CreateDescription(const I3CorsikaWeight& weight) {
  I3TableRowDescriptionPtr desc(new I3TableRowDescription());
  desc->AddField<double>("zenith", "radian","zenith angle of particle origin");
  desc->AddField<double>("azimuth","radian","azimuthal angle of particle origin");
  desc->AddField<double>("energy", "GeV",   "energy of particle");
  MAKE_ENUM_VECTOR(particletype,I3Particle,I3Particle::ParticleType,I3PARTICLE_H_I3Particle_ParticleType);
  MAKE_ENUM_VECTOR(biastype,I3ShowerBias,I3ShowerBias::BiasParticleType,(Mu)(NuMu)(NuE));
  desc->AddEnumField<I3Particle::ParticleType> ("type",      particletype,"","");
  desc->AddEnumField<I3ShowerBias::BiasParticleType>("bias_type",biastype,"","");
  desc->AddField<double>("bias_target","","");
  desc->AddField<double>("weight","","");
  desc->AddField<double>("max_x","","");
  return desc;
};

size_t I3CorsikaWeightConverter::FillRows(const I3CorsikaWeight& w, I3TableRowPtr rows) {
  rows->Set<double>("zenith", w.primary.GetZenith());
  rows->Set<double>("azimuth",w.primary.GetAzimuth());
  rows->Set<double>("energy", w.primary.GetEnergy());
  rows->Set<I3Particle::ParticleType> ("type", w.primary.GetType());
  rows->Set<I3ShowerBias::BiasParticleType> ("bias_type", w.bias.type);  
  rows->Set<double>("bias_target", w.bias.target);
  rows->Set<double>("weight", w.weight);
  rows->Set<double>("max_x", w.max_x);
  return 1;
};
