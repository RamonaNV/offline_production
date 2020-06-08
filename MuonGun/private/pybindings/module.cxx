/** $Id: module.cxx 162249 2018-04-19 19:01:37Z olivas $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 162249 $
 * $Date: 2018-04-19 13:01:37 -0600 (Thu, 19 Apr 2018) $
 */

#include <icetray/I3FrameObject.h>
#include <icetray/load_project.h>

#ifdef USE_NUMPY
#if BOOST_VERSION < 106300
#include <boost/numpy.hpp>
#define BOOST_NUMPY boost::numpy
#else
#include <boost/python/numpy.hpp>
#define BOOST_NUMPY boost::python::numpy
#endif // BOOST_VERSION < 106300
#endif // USE_NUMPY

namespace bp = boost::python;

#define REGISTER_THESE_THINGS                                           \
  (I3MuonGun)(RadialDistribution)(EnergyDistribution)                   \
  (Flux)(Generator)(Surface)(WeightCalculator)(CanCan)                  \
  (CORSIKAGenerationProbability)(Track)

#define REGISTER_EXTRA (histogram)(TrackBinner)(MuonPropagator)(CompactTrack)

#define I3_REGISTRATION_FN_DECL(r, data, t) void BOOST_PP_CAT(register_,t)();
#define I3_REGISTER(r, data, t) BOOST_PP_CAT(register_,t)();
BOOST_PP_SEQ_FOR_EACH(I3_REGISTRATION_FN_DECL, ~, REGISTER_THESE_THINGS)
#ifdef USE_PROPOSAL
BOOST_PP_SEQ_FOR_EACH(I3_REGISTRATION_FN_DECL, ~, REGISTER_EXTRA);
#endif

I3_PYTHON_MODULE(MuonGun)
{
#ifdef USE_NUMPY
	BOOST_NUMPY::initialize();
#endif
	load_project("MuonGun", false);
	BOOST_PP_SEQ_FOR_EACH(I3_REGISTER, ~, REGISTER_THESE_THINGS);
	#ifdef USE_PROPOSAL
	BOOST_PP_SEQ_FOR_EACH(I3_REGISTER, ~, REGISTER_EXTRA);
	#endif
}

