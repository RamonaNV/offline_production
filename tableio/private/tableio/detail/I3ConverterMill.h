/*
 * copyright  (C) 2012
 * The Icecube Collaboration
 *
 * $Id: I3ConverterMill.h 125933 2014-11-19 20:53:49Z jvansanten $
 *
 * @version $Revision: 125933 $
 * @date $LastChangedDate: 2014-11-19 13:53:49 -0700 (Wed, 19 Nov 2014) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy: jvansanten $
 */

#ifndef TABLEIO_DETAIL_I3CONVERTERMILL_H_INCLUDED
#define TABLEIO_DETAIL_I3CONVERTERMILL_H_INCLUDED

#include "tableio/I3Converter.h"
#include <boost/python/object.hpp>

/// @cond 
class I3ConverterMill {
public:
	I3ConverterMill(boost::python::object);
	bool CanConvert(I3FrameObjectConstPtr);
	I3ConverterPtr operator()();
private:
	boost::python::object callable_;
	I3ConverterPtr thneed_;
};

I3_POINTER_TYPEDEFS(I3ConverterMill);

/// @endcond

#endif // TABLEIO_DETAIL_I3CONVERTERMILL_H_INCLUDED