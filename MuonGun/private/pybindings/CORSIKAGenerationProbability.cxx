/** $Id: CORSIKAGenerationProbability.cxx 108206 2013-07-13 15:18:23Z nwhitehorn $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 108206 $
 * $Date: 2013-07-13 09:18:23 -0600 (Sat, 13 Jul 2013) $
 */

#include <MuonGun/CORSIKAGenerationProbability.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

void register_CORSIKAGenerationProbability()
{
	using namespace I3MuonGun;
	using namespace boost::python;
	
	class_<CORSIKAGenerationProbability, boost::shared_ptr<CORSIKAGenerationProbability>, bases<GenerationProbability> >
	    ("CORSIKAGenerationProbability", init<SamplingSurfacePtr, FluxPtr, RadialDistributionPtr, EnergyDistributionPtr>())
	    #define RO_PROPS (Flux)(RadialDistribution)(EnergyDistribution)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP_RO, CORSIKAGenerationProbability, RO_PROPS)
	    #undef RO_PROPS
	;
}
