/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3CorsikaShowerInfoConverter.cxx 179791 2020-04-16 16:29:12Z kath $
 *
 * @version $Revision: 179791 $
 * @date $LastChangedDate: 2020-04-16 10:29:12 -0600 (Thu, 16 Apr 2020) $
 * @author Fabian Kislat <fabian.kislat@desy.de> $LastChangedBy: kath $
 */

#include "simclasses/converter/I3CorsikaShowerInfoConverter.h"
#include <icetray/I3Units.h>

static const double gcm2 = I3Units::g/I3Units::cm2;

I3CorsikaShowerInfoConverter::I3CorsikaShowerInfoConverter()
  : I3ConverterImplementation<I3CorsikaShowerInfo>(),
    nLongSteps_(0)
{}

I3TableRowDescriptionPtr I3CorsikaShowerInfoConverter::CreateDescription(const I3CorsikaShowerInfo& info)
{
  I3TableRowDescriptionPtr desc = I3TableRowDescriptionPtr(new I3TableRowDescription() );

  /*
   * This implementation assumes that all longitudinal profiles have the same 
   * number of steps. This is the case for files based on a single CORSIKA binary
   * input file. Restriction: I do not know if this is still true in case of the
   * SLANT or CURVED options.
   */
  nLongSteps_ = info.longProfile.size();
  
  desc->AddField<int32_t>("crsRunID", "", "Corsika run number");
  desc->AddField<int32_t>("crsEventID", "", "Corsika event number");
  desc->AddField<int32_t>("crsSampleID", "", "Resampling index of this shower");
  desc->AddField<double>("firstIntHeight", "m", "Height of the first interaction");
  desc->AddField<double>("firstIntDepth", "g/cm^2", "Atmospheric depth of the first interaction");
  desc->AddField<double>("obsLevelHeight", "m", "Observation level (in case of multiple observation levels, only the lowest one is used)");
  desc->AddField<double>("ghMaxNum", "", "Gaisser-Hillas fit parameter: maximum number of particles");
  desc->AddField<double>("ghStartDepth", "g/cm^2", "Gaisser-Hillas fit parameter: starting depth of the shower");
  desc->AddField<double>("ghMaxDepth", "g/cm^2", "Gaisser-Hillas fit parameter: depth of the shower maximum");
  desc->AddField<double>("ghLambdaa", "g/cm^2", "Gaisser-Hillas fit parameter: \"a\" part of lambda = a + bt + ct^2");
  desc->AddField<double>("ghLambdab", "", "Gaisser-Hillas fit parameter: \"b\" part of lambda = a + bt + ct^2");
  desc->AddField<double>("ghLambdac", "1/(g/cm^2)", "Gaisser-Hillas fit parameter: \"c\" part of lambda = a + bt + ct^2");
  desc->AddField<double>("ghRedChiSqr", "", "Reduced chi-square of the Gaisser-Hillas fit");
  desc->AddField<double>("resampleRadius", "m", "Radius of the resampling area");
  desc->AddField<double>("weight", "", "Weighting factor due to resampling");
  desc->AddField<uint16_t>("nResample", "", "Number of samples created of the shower");
  desc->AddField<uint16_t>("nResampleNominal", "", "Nominal number of samples of the shower (what was requested).");

  if (nLongSteps_) {
    desc->AddField<tableio_size_t>("nLongSteps", "", "Number of steps in the longitudinal profile");
    desc->AddField<double>("longDepth", "g/cm^2", "Longitudinal profile: atmospheric depth", nLongSteps_);
    desc->AddField<uint64_t>("longNumGamma", "", "Longitudinal profile: number of gammas", nLongSteps_);
    desc->AddField<uint64_t>("longNumEMinus", "", "Longitudinal profile: number of electrons", nLongSteps_);
    desc->AddField<uint64_t>("longNumEPlus", "", "Longitudinal profile: number of positrons", nLongSteps_);
    desc->AddField<uint64_t>("longNumMuMinus", "", "Longitudinal profile: number of muons", nLongSteps_);
    desc->AddField<uint64_t>("longNumMuPlus", "", "Longitudinal profile: number of anti-muons", nLongSteps_);
    desc->AddField<uint64_t>("longNumHadron", "", "Longitudinal profile: number of hadrons", nLongSteps_);
    desc->AddField<uint64_t>("longNumCharged", "", "Longitudinal profile: number of charged particles", nLongSteps_);
    desc->AddField<uint64_t>("longNumNuclei", "", "Longitudinal profile: number of nuclei", nLongSteps_);
    desc->AddField<uint64_t>("longNumCherenkov", "", "Longitudinal profile: number of Cherenkov photons", nLongSteps_);
  }

  return desc;
}

size_t I3CorsikaShowerInfoConverter::FillRows(const I3CorsikaShowerInfo& info, I3TableRowPtr rows)
{
  // A temporary hack to see if we're reading files in which topsimulator stored depths in g/cm^2 (without converting).
  // We'll try to predict it and act accordingly.
  // This bit should go away in the future when all topsimulator files are storing Corsika-style units.
  // This routine will serve up numbers in "original Corsika-style units" such as g/cm^2 etc, and NOT I3Units.
  bool depths_not_stored_using_I3Units = 0;
  double big = 0.001 * gcm2;
  if (info.ghMaxDepth < big && info.ghStartDepth < big && info.firstIntDepth < big) {
    log_debug("It looks like depths were stored as g/cm^2 instead of I3Units native units. I will NOT convert them.");
    depths_not_stored_using_I3Units = 1;
  } else {
    log_info("This looks like an old file where depths were stored in I3Units native units. I will convert them to Corsika units.");
  }

  rows->Set<int32_t>("crsRunID", info.crsRunID);
  rows->Set<int32_t>("crsEventID", info.crsEventID);
  rows->Set<int32_t>("crsSampleID", info.crsSampleID);
  rows->Set<double>("firstIntHeight", info.firstIntHeight/I3Units::m);
  if (depths_not_stored_using_I3Units) 
    rows->Set<double>("firstIntDepth", info.firstIntDepth);
  else 
    rows->Set<double>("firstIntDepth", info.firstIntDepth/gcm2);
  rows->Set<double>("obsLevelHeight", info.obsLevelHeight/I3Units::m);
  rows->Set<double>("ghMaxNum", info.ghMaxNum);
  if (depths_not_stored_using_I3Units) {
    rows->Set<double>("ghStartDepth", info.ghStartDepth);
    rows->Set<double>("ghMaxDepth", info.ghMaxDepth);
  } else {
    rows->Set<double>("ghStartDepth", info.ghStartDepth/gcm2);
    rows->Set<double>("ghMaxDepth", info.ghMaxDepth/gcm2);
  }
  rows->Set<double>("ghLambdaa", info.ghLambdaa); // Corsika units (g/cm^2)
  rows->Set<double>("ghLambdab", info.ghLambdab); // Corsika units (unitless)
  rows->Set<double>("ghLambdac", info.ghLambdac); // Corsika units (1/(g/cm^2))
  rows->Set<double>("ghRedChiSqr", info.ghRedChiSqr);
  rows->Set<double>("resampleRadius", info.resampleRadius);
  rows->Set<double>("weight", info.weight);
  rows->Set<uint16_t>("nResample", info.nResample);
  rows->Set<uint16_t>("nResampleNominal", info.nResampleNominal);


  if (nLongSteps_) {
    if (info.longProfile.size() != nLongSteps_) {
      log_info("%s: This event has a longitudinal profile with %zu steps, while "
	       "the converter is configured to hold %zu longitudinal steps.",
	       __PRETTY_FUNCTION__, info.longProfile.size(), nLongSteps_);
    }
    rows->Set<tableio_size_t>("nLongSteps", nLongSteps_);
    double *longDepthBuffer = rows->GetPointer<double>("longDepth");
    uint64_t *longNumGammaBuffer = rows->GetPointer<uint64_t>("longNumGamma");
    uint64_t *longNumEMinusBuffer = rows->GetPointer<uint64_t>("longNumEMinus");
    uint64_t *longNumEPlusBuffer = rows->GetPointer<uint64_t>("longNumEPlus");
    uint64_t *longNumMuMinusBuffer = rows->GetPointer<uint64_t>("longNumMuMinus");
    uint64_t *longNumMuPlusBuffer = rows->GetPointer<uint64_t>("longNumMuPlus");
    uint64_t *longNumHadronBuffer = rows->GetPointer<uint64_t>("longNumHadron");
    uint64_t *longNumChargedBuffer = rows->GetPointer<uint64_t>("longNumCharged");
    uint64_t *longNumNucleiBuffer = rows->GetPointer<uint64_t>("longNumNuclei");
    uint64_t *longNumCherenkovBuffer = rows->GetPointer<uint64_t>("longNumCherenkov");

    for (size_t i = 0; i < std::min(info.longProfile.size(), nLongSteps_); ++i) {
      if (depths_not_stored_using_I3Units) 
	longDepthBuffer[i] = info.longProfile[i].depth;
      else 
	longDepthBuffer[i] = info.longProfile[i].depth/gcm2;
      longNumGammaBuffer[i] = info.longProfile[i].numGamma;
      longNumEMinusBuffer[i] = info.longProfile[i].numEMinus;
      longNumEPlusBuffer[i] = info.longProfile[i].numEPlus;
      longNumMuMinusBuffer[i] = info.longProfile[i].numMuMinus;
      longNumMuPlusBuffer[i] = info.longProfile[i].numMuPlus;
      longNumHadronBuffer[i] = info.longProfile[i].numHadron;
      longNumChargedBuffer[i] = info.longProfile[i].numCharged;
      longNumNucleiBuffer[i] = info.longProfile[i].numNuclei;
      longNumCherenkovBuffer[i] = info.longProfile[i].numCherenkov;
    }
  }

  return 1;
}
