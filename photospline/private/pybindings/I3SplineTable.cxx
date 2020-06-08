
#define NO_IMPORT_ARRAY /* Just use the headers */
#define PY_ARRAY_UNIQUE_SYMBOL photospline_PyArray_API
#include <numpy/numpyconfig.h>
#ifdef NPY_1_7_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif
#include <numpy/ndarrayobject.h>

#include "photospline/I3SplineTable.h"

namespace bp = boost::python;

#define PY_TYPESTRING(pyobj) \
	pyobj.ptr()->ob_type->tp_name

static double
splinetableeval(I3SplineTable &self, bp::object coordinates, bp::object derivatives)
{
	double retvalue(NAN);
	double *coord_ptr;
	unsigned *deriv_ptr=NULL;

	PyObject *coords;	
	coords = PyArray_ContiguousFromObject(coordinates.ptr(), NPY_DOUBLE, 0, 0);
	if (!coords) {
		PyErr_Format(PyExc_ValueError, "Can't convert object of type"
		    "'%s' to an array of doubles!", PY_TYPESTRING(coordinates));
		bp::throw_error_already_set();
	}
	coord_ptr = (double*)PyArray_DATA((PyArrayObject *)coords);
	
	PyObject *derivs=NULL;
	if (derivatives.ptr() != Py_None) {
		// This conversion will wrap around if the contents of derivatives are
		// negative, resulting in large numbers of differentiations. As long
		// as the splines are of reasonably small order, though, the result of
		// the spline evaluation will end up being zero.
		derivs = PyArray_ContiguousFromObject(derivatives.ptr(), NPY_UINT, 0, 0);
		if (!derivs) {
			PyErr_Format(PyExc_ValueError, "Can't convert object of type"
			    "'%s' to an array of integers!", PY_TYPESTRING(derivatives));
			bp::throw_error_already_set();
		}
		deriv_ptr = (unsigned*)PyArray_DATA((PyArrayObject *)derivs);
	}
	
	self.Eval(coord_ptr, &retvalue, deriv_ptr);
	Py_XDECREF(coords);
	if (derivs != NULL)
		Py_XDECREF(derivs);

	return retvalue;
}

static void
splinetableconvolve(I3SplineTable &self, int dim, bp::object knots)
{
	double *knot_ptr;
	double retvalue(NAN);
	double *coord_ptr;
	unsigned *deriv_ptr=NULL;

	PyObject *knot_array = PyArray_ContiguousFromObject(knots.ptr(), NPY_DOUBLE, 0, 0);
	if (!knot_array) {
		PyErr_Format(PyExc_ValueError, "Can't convert object of type"
		    "'%s' to an array of doubles!", PY_TYPESTRING(knots));
		bp::throw_error_already_set();
	}
	knot_ptr = (double*)PyArray_DATA((PyArrayObject *)knot_array);
	
	self.Convolve(dim, knot_ptr, PyArray_SIZE((PyArrayObject *)knot_array));
	
	Py_XDECREF(knot_array);
}


static bp::list
GetExtents(const I3SplineTable &self)
{
	bp::list extents;
	for (unsigned i=0; i < self.GetNDim(); i++) {
		std::pair<double, double> ext = self.GetExtents(i);
		extents.append(bp::make_tuple(ext.first, ext.second));
	}
	
	return extents;
}

void register_I3SplineTable() {
	bp::class_<I3SplineTable, boost::shared_ptr<I3SplineTable>, boost::noncopyable>
	    ("I3SplineTable", bp::init<const std::string&>(bp::arg("path")))
	    .def("convolve", splinetableconvolve, (bp::args("dimension"), "knots"))
	    .def("eval", splinetableeval, (bp::args("coordinates"), bp::arg("derivatives")=bp::object()),
	        "Evaluate the spline surface at the given coordinates.\n\n"
	        ":param coordinates: N-dimensonal coordinates at which to evaluate\n"
	        ":param derivatives: An array indicating the type of basis to use"
	                          "in each dimension. If an entry corresponding to"
	                          "a dimension is N>0, the basis in that"
	                          "dimension will consist of the Nth derivatives of"
	                          "the usual B-spline basis, and result"
	                          "will be the Nth-order gradient of the surface in"
	                          "that dimension. If NULL, 0 will be assumed.")
	    .add_property("ndim", &I3SplineTable::GetNDim)
	    .add_property("extents", &GetExtents)
	;
}

