/** $Id: WeightCalculator.cxx 161217 2018-02-26 14:28:37Z thomas.kintscher $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 161217 $
 * $Date: 2018-02-26 07:28:37 -0700 (Mon, 26 Feb 2018) $
 */

#include <MuonGun/WeightCalculator.h>
#include <MuonGun/SamplingSurface.h>
#include <MuonGun/Cylinder.h>
#include <MuonGun/Flux.h>
#include <MuonGun/RadialDistribution.h>
#include <MuonGun/EnergyDistribution.h>

#include <tableio/converter/pybindings.h>

using namespace I3MuonGun;
using namespace boost::python;
 
#ifdef USE_NUMPY
#if BOOST_VERSION < 106300
#include <boost/numpy.hpp>
#define BOOST_NUMPY boost::numpy
#else
#include <boost/python/numpy.hpp>
#define BOOST_NUMPY boost::python::numpy
#endif // BOOST_VERSION < 106300

namespace {

template <typename T>
inline T
get(const BOOST_NUMPY::ndarray &a, int i0)
{
	return *reinterpret_cast<T*>(a.get_data() + i0*a.strides(0));
}

template <typename T>
inline T
get(const BOOST_NUMPY::ndarray &a, int i0, int i1)
{
	return *reinterpret_cast<T*>(a.get_data() + i0*a.strides(0) + i1*a.strides(1));
}

}

// An adapter function to use the standard WeightCalculator with Numpy arrays,
// most likely from tableio/hdfwriter
object
GetWeight(const WeightCalculator& weighter, object &xo, object &yo, object &zo,
    object &zeno, object &azio, object &mo, object &eno, object &rado)
{
	using namespace BOOST_NUMPY;
	
	// For efficiency's sake, match the types emitted by the default tableio converters
	ndarray x   = from_object(xo,   dtype::get_builtin<double>(), 0, 1);
	ndarray y   = from_object(yo,   dtype::get_builtin<double>(), 0, 1);
	ndarray z   = from_object(zo,   dtype::get_builtin<double>(), 0, 1);
	ndarray zen = from_object(zeno, dtype::get_builtin<double>(), 0, 1);
	ndarray azi = from_object(azio, dtype::get_builtin<double>(), 0, 1);
	ndarray mult     = from_object(mo,   dtype::get_builtin<uint32_t>(), 0, 1);
	ndarray energies = from_object(eno,  dtype::get_builtin<float>(),    1, 2);
	ndarray radii    = from_object(rado, dtype::get_builtin<float>(),    1, 2);
	
	// Reshape if we get scalar arguments
	if (x.get_nd() == 0) {
		x   = x.reshape(make_tuple(1));
		y   = y.reshape(make_tuple(1));
		z   = z.reshape(make_tuple(1));
		zen = zen.reshape(make_tuple(1));
		azi = azi.reshape(make_tuple(1));
		mult = mult.reshape(make_tuple(1));
		energies = energies.reshape(make_tuple(1, energies.shape(0)));
		radii = radii.reshape(make_tuple(1, radii.shape(0)));
	}
	
	int nrows = x.shape(0);
	int ncols = std::min(radii.shape(1), energies.shape(1));
	if (y.shape(0) != nrows || z.shape(0) != nrows
	    || zen.shape(0) != nrows || azi.shape(0) != nrows
	    || energies.shape(0) != nrows || radii.shape(0) != nrows)
	    throw(std::runtime_error("shape mismatch!"));
	
	ndarray weights = zeros(x.get_nd(), x.get_shape(), dtype::get_builtin<double>());	
	I3Particle axis;
	for (int i=0; i < nrows; i++) {
		
		int m = get<uint32_t>(mult, i);
		BundleConfiguration spec;
		for (int j=0; j < std::min(ncols, m); j++) {
			spec.push_back(BundleEntry(
			    get<float>(radii, i, j), get<float>(energies, i, j)));
		}
		axis.SetPos(get<double>(x, i), get<double>(y, i), get<double>(z, i));
		axis.SetDir(get<double>(zen, i), get<double>(azi, i));
		
		reinterpret_cast<double*>(weights.get_data())[i] = weighter.GetWeight(axis, spec);
	}
	
	return weights.scalarize();
}

#endif

void register_WeightCalculator()
{
	def("muons_at_surface", &GetMuonsAtSurface);

	class_<BundleModel>("BundleModel", init<FluxPtr, RadialDistributionPtr, EnergyDistributionPtr>(
	    (arg("flux"), "radius", "energy")))
	    .def_readwrite("flux", &BundleModel::flux)
	    .def_readwrite("radius", &BundleModel::radius)
	    .def_readwrite("energy", &BundleModel::energy)
	;

	class_<WeightCalculator>("WeightCalculator", init<const BundleModel&, GenerationProbabilityPtr>(
	    (arg("model"), "generator")))
	    .def("__call__", &WeightCalculator::GetWeight)
#ifdef USE_NUMPY
	    .def("__call__", &GetWeight, (bp::arg("x"), "y", "z", "zenith", "azimuth",
	        "multiplicity", "energies", "radii"))
#endif
	    #define PROPS (Surface)
	    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, WeightCalculator, PROPS)
	    #undef PROPS
	;
	
	I3CONVERTER_NAMESPACE(MuonGun);
	I3CONVERTER_EXPORT(MuonBundleConverter, "foo")
	    .def(init<uint32_t, SamplingSurfaceConstPtr>((
	    arg("maxMultiplicity")=25,
	    arg("surface")=boost::make_shared<Cylinder>(1600, 800))))
	;
}
