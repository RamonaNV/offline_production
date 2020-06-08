/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: CorsikaLongStep.cxx 126444 2014-12-02 00:19:10Z hdembinski $
 *
 * @version $Revision: 126444 $
 * @date $LastChangedDate: 2014-12-01 17:19:10 -0700 (Mon, 01 Dec 2014) $
 * @author Fabian Kislat <fabian.kislat@desy.de> $LastChangedBy: hdembinski $
 */

#include "simclasses/CorsikaLongStep.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <icetray/python/dataclass_suite.hpp>

namespace bp = boost::python;

void register_CorsikaLongStep()
{
  bp::class_<CorsikaLongStep, CorsikaLongStepPtr>("CorsikaLongStep")
    .def(bp::dataclass_suite<CorsikaLongStep>())
    .def_readwrite("depth", &CorsikaLongStep::depth)
    .def_readwrite("numGamma", &CorsikaLongStep::numGamma)
    .def_readwrite("numEMinus", &CorsikaLongStep::numEMinus)
    .def_readwrite("numEPlus", &CorsikaLongStep::numEPlus)
    .def_readwrite("numMuMinus", &CorsikaLongStep::numMuMinus)
    .def_readwrite("numMuPlus", &CorsikaLongStep::numMuPlus)
    .def_readwrite("numHadron", &CorsikaLongStep::numHadron)
    .def_readwrite("numCharged", &CorsikaLongStep::numCharged)
    .def_readwrite("numNuclei", &CorsikaLongStep::numNuclei)
    .def_readwrite("numCherenkov", &CorsikaLongStep::numCherenkov)
    ;

  bp::class_< std::vector<CorsikaLongStep> >("CorsikaLongProfile")
    .def(bp::dataclass_suite<std::vector<CorsikaLongStep> >())
    ;
}
