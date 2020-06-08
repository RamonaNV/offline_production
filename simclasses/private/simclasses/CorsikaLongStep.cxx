/**
 * Copyright (C) 2010
 * The IceCube collaboration
 * ID: $Id: $
 *
 * @file CorsikaLongStep.cxx
 * @version $Rev: $
 * @date $Date: $
 * @author Tilo Waldenmaier
 */

#include "simclasses/CorsikaLongStep.h"
#include <icetray/serialization.h>
#include <iomanip>

CorsikaLongStep::CorsikaLongStep():
  depth(NAN),
  numGamma(0),
  numEMinus(0),
  numEPlus(0),
  numMuMinus(0),
  numMuPlus(0),
  numHadron(0),
  numCharged(0),
  numNuclei(0),
  numCherenkov(0)
{
  
}

CorsikaLongStep::~CorsikaLongStep()
{
  
}

template <class Archive>
void CorsikaLongStep::serialize(Archive& ar, unsigned version)
{
  if(version>corsikalongstep_version_)
  {
    log_fatal("Attempting to read version %u from file but running version %u of CorsikaLongStep class.",
	      version,
	      corsikalongstep_version_);
  }
  
  ar & make_nvp("depth",        depth);
  ar & make_nvp("numGamma",     numGamma);
  ar & make_nvp("numEMinus",    numEMinus);
  ar & make_nvp("numEPlus",     numEPlus);
  ar & make_nvp("numMuMinus",   numMuMinus);
  ar & make_nvp("numMuPlus",    numMuPlus);
  ar & make_nvp("numHadron",    numHadron);
  ar & make_nvp("numCharged",   numCharged);
  ar & make_nvp("numNuclei",    numNuclei);
  ar & make_nvp("numCherenkov", numCherenkov);
}

bool CorsikaLongStep::operator==(const CorsikaLongStep& rhs) const
{
  return 
    (depth == rhs.depth) && 
    (numGamma == rhs.numGamma) &&
    (numEMinus == rhs.numEMinus) &&
    (numEPlus == rhs.numEPlus) &&
    (numMuMinus == rhs.numMuMinus) &&
    (numMuPlus == rhs.numMuPlus) &&
    (numHadron == rhs.numHadron) &&
    (numCharged == rhs.numCharged) &&
    (numNuclei == rhs.numNuclei) &&
    (numCherenkov == rhs.numCherenkov)
    ;
}

std::ostream& operator<<(std::ostream& os, const CorsikaLongStep& s) {
  os << "[ CorsikaLongStep::"
     << "\n  depth       :" << s.depth
     << "\n  numGamma    :" << s.numGamma
     << "\n  numEMinus   :" << s.numEMinus
     << "\n  numEPlus    :" << s.numEPlus
     << "\n  numMuMinus  :" << s.numMuMinus
     << "\n  numMuPlus   :" << s.numMuPlus
     << "\n  numHadron   :" << s.numHadron
     << "\n  numCharged  :" << s.numCharged
     << "\n  numNuclei   :" << s.numNuclei
     << "\n  numCherenkov:" << s.numCherenkov
     << " ]";
  return os;
}

std::ostream& operator<<(std::ostream& os, const CorsikaLongProfile& p) {
  using std::setw;
  os << "[ CorsikaLongProfile::";
  os << "\n  " 
     << setw(6) << "depth" 
     << setw(13) << "numGamma" 
     << setw(13) << "numEMinus"
     << setw(13) << "numEPlus"
     << setw(13) << "numMuMinus"
     << setw(13) << "numMuPlus"
     << setw(13) << "numHadron"
     << setw(13) << "numCharged"
     << setw(13) << "numNuclei"
     << setw(13) << "numCherenkov";
  for (CorsikaLongProfile::const_iterator
         iter = p.begin(), end = p.end();
       iter != end; ++iter)
    os << "\n  " 
       << setw(6) << iter->depth 
       << setw(13) << iter->numGamma
       << setw(13) << iter->numEMinus
       << setw(13) << iter->numEPlus
       << setw(13) << iter->numMuMinus
       << setw(13) << iter->numMuPlus
       << setw(13) << iter->numHadron
       << setw(13) << iter->numCharged
       << setw(13) << iter->numNuclei
       << setw(13) << iter->numCherenkov;
  os << " ]";
  return os;
}

I3_SERIALIZABLE(CorsikaLongStep);
