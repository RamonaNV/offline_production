/**
 * Copyright (C) 2010
 * The IceCube collaboration
 * ID: $Id: $
 *
 * @file CorsikaLongStep.h
 * @version $Rev: $
 * @date $Date: $
 * @author Tilo Waldenmaier
 */

#ifndef _TOPSIMULATOR_CORSIKALONGSTEP_H_
#define _TOPSIMULATOR_CORSIKALONGSTEP_H_

#include <vector>
#include <dataclasses/Utility.h>

static const unsigned corsikalongstep_version_ = 0;

struct CorsikaLongStep
{
  double depth;
  uint64_t numGamma;
  uint64_t numEMinus;
  uint64_t numEPlus;
  uint64_t numMuMinus;
  uint64_t numMuPlus;
  uint64_t numHadron;
  uint64_t numCharged;
  uint64_t numNuclei;
  uint64_t numCherenkov;
  
  CorsikaLongStep();
  
  virtual ~CorsikaLongStep();

  bool operator==(const CorsikaLongStep& rhs) const;
  
  template <class Archive> void serialize(Archive & ar, unsigned version);
};

typedef std::vector<CorsikaLongStep> CorsikaLongProfile;

std::ostream& operator<<(std::ostream&, const CorsikaLongStep&);
std::ostream& operator<<(std::ostream&, const CorsikaLongProfile&);

I3_CLASS_VERSION(CorsikaLongStep, corsikalongstep_version_);
I3_POINTER_TYPEDEFS(CorsikaLongStep);
I3_POINTER_TYPEDEFS(CorsikaLongProfile);

#endif
