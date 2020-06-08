/*
 * copyright  (C) 2012
 * The Icecube Collaboration
 *
 * $Id: I3ConverterMill.cxx 125933 2014-11-19 20:53:49Z jvansanten $
 *
 * @version $Revision: 125933 $
 * @date $LastChangedDate: 2014-11-19 13:53:49 -0700 (Wed, 19 Nov 2014) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy: jvansanten $
 */

#include "tableio/detail/I3ConverterMill.h"
#include <boost/python/extract.hpp>

/// @cond
I3ConverterMill::I3ConverterMill(boost::python::object callable)
    : callable_(callable)
{
	thneed_ = (*this)();
	if (!thneed_)
		throw std::runtime_error("Couldn't instantiate converter!");
}
	
bool
I3ConverterMill::CanConvert(I3FrameObjectConstPtr object)
{
	return thneed_->CanConvert(object);
}
	
I3ConverterPtr
I3ConverterMill::operator()()
{
	return boost::python::extract<I3ConverterPtr>(callable_());
}
/// @endcond
