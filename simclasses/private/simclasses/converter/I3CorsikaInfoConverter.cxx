/**
 * copyright  (C) 2020 The Icecube Collaboration
 *
 * $Id: $
 *
 * @version $Revision: $
 * @date $LastChangedDate:$
 * @author Kevin Meagher $LastChangedBy: $
 */

#include "simclasses/converter/I3CorsikaInfoConverter.h"

I3TableRowDescriptionPtr I3CorsikaInfoConverter::CreateDescription(const I3CorsikaInfo& info)
{
  I3TableRowDescriptionPtr desc = I3TableRowDescriptionPtr(new I3TableRowDescription() );
  desc->AddField<uint32_t>("run_id","","Run number for this neutrino generatior job");
  desc->AddField<uint32_t>("n_events","","Number of primary neutrino events produced by this job");
  MAKE_ENUM_VECTOR(primary_type,I3Particle,I3Particle::ParticleType,I3PARTICLE_H_I3Particle_ParticleType);
  desc->AddEnumField<I3Particle::ParticleType> ("primary_type",primary_type,"","PDG encoding of the primary neutrino");
  desc->AddField<uint32_t> ("atmosphere","","");
  desc->AddField<uint32_t> ("oversampling","","Number of times an air shower is repeated in the detector");  
  desc->AddField<double> ("cylinder_height","","Hight of cylinder particles are injected into");
  desc->AddField<double> ("cylinder_radius","","Radius of cylinder particles are injected into");
  desc->AddField<double> ("min_zenith","","Minimun zenith angle of angular distribution");
  desc->AddField<double> ("max_zenith","","Maximum zenigh angle of angular distribution");
  desc->AddField<double> ("min_energy","","Minimum energy of energy distribution");
  desc->AddField<double> ("max_energy","","Minimum energy of energy distribution");
  desc->AddField<double> ("power_law_index","","Exponent of energy distribution");
  return desc;
}

size_t I3CorsikaInfoConverter::FillRows(const I3CorsikaInfo& info, I3TableRowPtr rows)
{
  rows->Set<uint32_t>("run_id", info.run_id);
  rows->Set<uint32_t>("n_events",info.n_events);
  rows->Set<I3Particle::ParticleType> ("primary_type",info.primary_type);
  rows->Set<uint32_t> ("atmosphere",info.atmosphere);
  rows->Set<uint32_t> ("oversampling",info.oversampling);  
  rows->Set<double> ("cylinder_height",info.cylinder_height);
  rows->Set<double> ("cylinder_radius",info.cylinder_radius);
  rows->Set<double> ("min_zenith",info.min_zenith);
  rows->Set<double> ("max_zenith",info.max_zenith);
  rows->Set<double> ("min_energy",info.min_energy);
  rows->Set<double> ("max_energy",info.max_energy);
  rows->Set<double> ("power_law_index",info.power_law_index);
  return 1;
}

